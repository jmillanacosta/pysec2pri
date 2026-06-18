"""Build/read the per-mapping first-seen-date index for any versioned datasource.

Some datasources (HGNC, NCBI) carry a real per-row change date in their flat
files (``date_symbol_changed``, ``Discontinue_Date``) -- those are already
wired directly into ``mapping_date`` during normal parsing and don't need
this module at all. Others (ChEBI, UniProt) carry no per-row date anywhere:
a given ``(primary, secondary)`` pair's deprecation date isn't recorded in
any single release. Since :meth:`~pysec2pri.parsers.base.BaseParser._record_id`
hashes a mapping's ``record_id`` without baking in the release version, the
same pair hashes identically across every release of a given datasource --
so walking every historical release once and recording the first release a
``record_id`` appears in recovers its true historical date.

This module supports two consolidation strategies, selected via *mode*:

- ``"release"``: walk every historical release (oldest first) and record
  the version/date each mapping first/last appeared in. Requires the
  datasource to publish a versioned archive (ChEBI, HGNC, UniProt).
- ``"date"``: take whatever real per-row ``mapping_date`` the datasource's
  parser already produces from a single (current) parse -- no historical
  walk needed. If the datasource never produces a real per-row date (ChEBI,
  UniProt), this falls back to ``"release"`` mode with a warning.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import warnings
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any

from pysec2pri.logging import logger
from pysec2pri.parsers.base import _cmp_versions

__all__ = [
    "SUPPORTED_DATASOURCES",
    "consolidate_mapping_dates",
    "default_cache_dir",
    "load_mapping_dates",
]

# Datasources this module knows how to download+parse per release.
# NCBI has no versioned archive, so only mode="date" applies to it.
SUPPORTED_DATASOURCES = ("chebi", "hgnc", "ncbi", "uniprot")

# mapping_sets kinds each datasource's parser actually supports.
_SUPPORTED_MAPPING_SETS: dict[str, tuple[str, ...]] = {
    "chebi": ("ids", "labels"),
    "hgnc": ("ids", "labels"),
    "ncbi": ("ids", "labels"),
    "uniprot": ("ids",),
}


def _chebi_versions(subset: str) -> list[str]:
    from pysec2pri.parsers.chebi import ChEBIDownloader

    return ChEBIDownloader(subset=subset).list_versions()


def _hgnc_versions(subset: str) -> list[str]:
    from pysec2pri.parsers.hgnc import HGNCDownloader

    return HGNCDownloader().list_versions()


def _uniprot_versions(subset: str) -> list[str]:
    from pysec2pri.parsers.uniprot import UniProtDownloader

    return UniProtDownloader().list_versions()


# Datasources with a versioned archive that can be walked for mode="release".
_LIST_VERSIONS_FNS: dict[str, Callable[[str], list[str]]] = {
    "chebi": _chebi_versions,
    "hgnc": _hgnc_versions,
    "uniprot": _uniprot_versions,
}


def default_cache_dir() -> Path:
    """Return the default cache directory for consolidated mapping-date indexes.

    Returns:
        ``$PYSEC2PRI_CACHE_DIR`` when set, otherwise ``~/.cache/pysec2pri``.
    """
    env = os.environ.get("PYSEC2PRI_CACHE_DIR")
    return Path(env) if env else Path.home() / ".cache" / "pysec2pri"


def _cache_path(cache_dir: Path, datasource: str, subset: str, mapping_sets: str) -> Path:
    return cache_dir / f"{datasource}_{subset}_{mapping_sets}_mapping_dates.tsv"


def _meta_path(cache_dir: Path, datasource: str, subset: str, mapping_sets: str) -> Path:
    return cache_dir / f"{datasource}_{subset}_{mapping_sets}_mapping_dates.meta.json"


_CACHE_COLUMNS = (
    "record_id",
    "first_seen_version",
    "first_seen_date",
    "last_seen_version",
    "last_seen_date",
)


def _read_cache(cache_path: Path) -> dict[str, dict[str, str]]:
    """Read a consolidated mapping-date cache TSV into a dict keyed by record_id."""
    if not cache_path.exists():
        return {}

    import polars as pl

    df = pl.read_csv(cache_path, separator="\t", schema_overrides={"record_id": pl.Utf8})
    return {
        str(row["record_id"]): {col: str(row[col]) for col in _CACHE_COLUMNS[1:]}
        for row in df.iter_rows(named=True)
    }


def _write_cache(cache_path: Path, records: dict[str, dict[str, str]]) -> None:
    """Write the merged ``record_id -> first/last seen`` dict to a TSV file."""
    import polars as pl

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    rows = [{"record_id": rid, **fields} for rid, fields in records.items()]
    schema = list(_CACHE_COLUMNS)
    df = pl.DataFrame(rows, schema=schema) if rows else pl.DataFrame(schema=schema)
    df.write_csv(cache_path, separator="\t")


def _read_meta(meta_path: Path) -> str | None:
    """Read the ``last_version`` sidecar, or ``None`` if absent/unreadable."""
    if not meta_path.exists():
        return None
    try:
        data: dict[str, Any] = json.loads(meta_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    last_version = data.get("last_version")
    return str(last_version) if last_version is not None else None


def _write_meta(meta_path: Path, last_version: str) -> None:
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    meta_path.write_text(json.dumps({"last_version": last_version}), encoding="utf-8")


def load_mapping_dates(
    datasource: str,
    cache_dir: Path | None = None,
    *,
    subset: str = "3star",
    mapping_sets: str = "ids",
) -> dict[str, str]:
    """Load the consolidated ``record_id -> first_seen_date`` index for *datasource*.

    Safe to call even when the index hasn't been built yet: returns ``{}``,
    which leaves ``Mapping.mapping_date`` unset for every row (the caller
    falls back to the set-level pinned release date, exactly as before this
    feature existed).

    Args:
        datasource: Datasource name (e.g. ``"chebi"``, ``"uniprot"``).
        cache_dir: Directory holding the cache file. Defaults to
            :func:`default_cache_dir`.
        subset: ``"3star"`` or ``"complete"`` (only meaningful for ChEBI) --
            must match the subset used when the index was built.
        mapping_sets: ``"ids"`` or ``"labels"`` -- must match the
            mapping-set kind used when the index was built.

    Returns:
        Dict mapping each ``record_id`` to its first-seen ISO date string.
    """
    cache_path = _cache_path(cache_dir or default_cache_dir(), datasource, subset, mapping_sets)
    records = _read_cache(cache_path)
    return {rid: fields["first_seen_date"] for rid, fields in records.items()}


def _parse_mapping_set(
    datasource: str,
    files: dict[str, Path],
    version: str | None,
    subset: str,
    mapping_sets: str,
    tax_id: str,
) -> Any:
    """Parse one downloaded release into a mapping set, dispatched by datasource.

    Calls each parser with explicit paths from *files* rather than relying
    on directory auto-discovery, so consolidation never depends on the
    tmpdir-guessing logic used by the ``generate_*`` convenience functions.
    """
    if datasource == "chebi":
        return _parse_chebi_version(files, version, subset, mapping_sets)
    if datasource == "hgnc":
        return _parse_hgnc_version(files, version, mapping_sets)
    if datasource == "ncbi":
        return _parse_ncbi_version(files, version, mapping_sets, tax_id)
    if datasource == "uniprot":
        return _parse_uniprot_version(files, version, mapping_sets)
    raise ValueError(f"Unsupported datasource for consolidation: {datasource!r}")


def _parse_chebi_version(
    files: dict[str, Path], version: str | None, subset: str, mapping_sets: str
) -> Any:
    from pysec2pri.parsers.chebi import ChEBIParser

    parser = ChEBIParser(version=version, show_progress=False, subset=subset)
    if "sdf" in files:
        return (
            parser.parse(files["sdf"])
            if mapping_sets == "ids"
            else parser.parse_synonyms(files["sdf"])
        )
    if mapping_sets == "ids":
        return parser.parse(
            secondary_ids_path=files["secondary_ids"],
            compounds_path=files.get("compounds"),
        )
    return parser.parse_synonyms(
        names_path=files["names"],
        compounds_path=files.get("compounds"),
    )


def _parse_hgnc_version(files: dict[str, Path], version: str | None, mapping_sets: str) -> Any:
    from pysec2pri.parsers.hgnc import HGNCParser

    parser = HGNCParser(version=version, show_progress=False)
    if mapping_sets == "ids":
        return parser.parse(files["withdrawn"], complete_set_path=files["complete"])
    return parser.parse_labels(files["complete"])


def _parse_ncbi_version(
    files: dict[str, Path], version: str | None, mapping_sets: str, tax_id: str
) -> Any:
    from pysec2pri.parsers.ncbi import NCBIParser

    parser = NCBIParser(version=version, show_progress=False)
    if mapping_sets == "ids":
        return parser.parse(files["gene_history"], tax_id=tax_id, gene_info_path=files["gene_info"])
    return parser.parse_labels(files["gene_info"], tax_id=tax_id)


def _parse_uniprot_version(files: dict[str, Path], version: str | None, mapping_sets: str) -> Any:
    from pysec2pri.parsers.uniprot import UniProtParser

    if mapping_sets != "ids":
        raise ValueError("uniprot only supports mapping_sets='ids'")
    parser = UniProtParser(version=version, show_progress=False)
    return parser.parse(files.get("sec_ac"), delac_path=files.get("delac_sp"))


def _run_one_version(
    datasource: str,
    version: str | None,
    subset: str,
    mapping_sets: str,
    tax_id: str,
) -> Any:
    """Download one release into a scratch tmpdir, parse it, then clean up.

    The tmpdir is always removed afterwards: walking ~250 releases must not
    accumulate downloaded files on disk.
    """
    from pysec2pri.download import download_datasource_with_release

    tmpdir = Path(tempfile.mkdtemp(prefix=f"pysec2pri_consolidate_{datasource}_"))
    try:
        files, _ = download_datasource_with_release(
            datasource, tmpdir, version=version, subset=subset
        )
        return _parse_mapping_set(datasource, files, version, subset, mapping_sets, tax_id)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _consolidate_by_date(
    datasource: str,
    cache_path: Path,
    subset: str,
    mapping_sets: str,
    tax_id: str,
) -> bool:
    """Single-pass "date" mode: capture each row's own real ``mapping_date``.

    Returns ``True`` and writes the cache when *datasource* produced at
    least one real per-row date; returns ``False`` (cache untouched) when it
    produced none, so the caller can fall back to ``"release"`` mode.
    """
    mapping_set = _run_one_version(datasource, None, subset, mapping_sets, tax_id)
    dated = [m for m in (mapping_set.mappings or []) if getattr(m, "mapping_date", None)]
    if not dated:
        return False

    version_label = str(getattr(mapping_set, "mapping_set_version", None) or "current")
    records: dict[str, dict[str, str]] = {}
    for m in dated:
        record_id = str(getattr(m, "record_id", None) or "")
        if not record_id:
            continue
        date_str = str(m.mapping_date)
        records[record_id] = {
            "first_seen_version": version_label,
            "first_seen_date": date_str,
            "last_seen_version": version_label,
            "last_seen_date": date_str,
        }
    _write_cache(cache_path, records)
    return True


def _consolidate_by_release(
    datasource: str,
    cache_path: Path,
    meta_path: Path,
    subset: str,
    mapping_sets: str,
    tax_id: str,
    show_progress: bool,
    force: bool,
) -> None:
    """Historical-walk "release" mode: track first/last-seen release per mapping."""
    from pysec2pri.download import resolve_release_date

    list_versions_fn = _LIST_VERSIONS_FNS.get(datasource)
    if list_versions_fn is None:
        raise ValueError(f"{datasource!r} has no versioned archive; only mode='date' is supported.")

    records: dict[str, dict[str, str]] = {} if force else _read_cache(cache_path)
    last_version = None if force else _read_meta(meta_path)

    versions = list_versions_fn(subset)
    if last_version is not None:
        versions = [v for v in versions if _cmp_versions(v, last_version) > 0]

    iterator: Iterable[str] = versions
    if show_progress:
        from tqdm import tqdm

        iterator = tqdm(versions, desc=f"Consolidating {datasource.upper()} mapping dates")

    for v in iterator:
        try:
            mapping_set = _run_one_version(datasource, v, subset, mapping_sets, tax_id)
            release_date = resolve_release_date(datasource, v, subset=subset)
            date_str = release_date.date().isoformat() if release_date else v

            for m in mapping_set.mappings or []:
                record_id = str(getattr(m, "record_id", None) or "")
                if not record_id:
                    continue
                entry = records.setdefault(
                    record_id,
                    {"first_seen_version": v, "first_seen_date": date_str},
                )
                entry["last_seen_version"] = v
                entry["last_seen_date"] = date_str
        except Exception:
            logger.warning(
                "Skipping %s version %s during consolidation", datasource, v, exc_info=True
            )
            continue

        last_version = v
        _write_cache(cache_path, records)
        _write_meta(meta_path, last_version)


def consolidate_mapping_dates(
    datasource: str,
    *,
    mode: str = "release",
    cache_dir: Path | None = None,
    subset: str = "3star",
    mapping_sets: str = "ids",
    tax_id: str = "9606",
    show_progress: bool = True,
    force: bool = False,
) -> Path:
    """Build/update the first-seen-date index for *datasource*.

    Two strategies, selected via *mode*:

    - ``"release"``: walk every historical release of *datasource* once
      (oldest first, resuming from the last completed version unless
      *force* is set). For every mapping seen, records the version/date it
      first appeared and keeps bumping the version/date it was last seen.
      This is a slow, network-heavy operation (~250 releases for ChEBI) and
      is meant to be run manually/as a one-off, not as part of normal
      mapping generation. Requires a versioned archive (ChEBI, HGNC,
      UniProt); raises :class:`ValueError` for datasources without one
      (e.g. NCBI).
    - ``"date"``: a single (current) parse, capturing whatever real
      per-row ``mapping_date`` the datasource's parser already produces
      (e.g. HGNC's ``date_symbol_changed``, NCBI's ``Discontinue_Date``).
      Fast -- no historical walk. If *datasource* never produces a real
      per-row date (ChEBI, UniProt), falls back to ``"release"`` mode and
      warns the user.

    Args:
        datasource: One of :data:`SUPPORTED_DATASOURCES`.
        mode: ``"release"`` or ``"date"``.
        cache_dir: Directory to write the cache file. Defaults to
            :func:`default_cache_dir`.
        subset: ``"3star"`` or ``"complete"`` (only meaningful for ChEBI).
        mapping_sets: ``"ids"`` or ``"labels"``.
        tax_id: NCBI taxonomy ID to filter (only meaningful for NCBI).
        show_progress: Whether to show a progress bar over releases
            (``"release"`` mode only).
        force: Re-scan every release from scratch, ignoring any existing
            cache/resume state (``"release"`` mode only).

    Returns:
        Path to the written cache TSV.

    Raises:
        ValueError: For an unknown *mode*, an unsupported *datasource*, an
            unsupported *mapping_sets* for *datasource*, or
            ``mode="release"`` against a datasource with no versioned
            archive.
    """
    if mode not in ("release", "date"):
        raise ValueError(f"mode must be 'release' or 'date', got {mode!r}")
    if datasource not in SUPPORTED_DATASOURCES:
        raise ValueError(
            f"Unsupported datasource: {datasource!r}. Supported: {SUPPORTED_DATASOURCES}"
        )
    if mapping_sets not in _SUPPORTED_MAPPING_SETS[datasource]:
        raise ValueError(
            f"{datasource!r} does not support mapping_sets={mapping_sets!r}. "
            f"Supported: {_SUPPORTED_MAPPING_SETS[datasource]}"
        )

    cache_dir = cache_dir or default_cache_dir()
    cache_path = _cache_path(cache_dir, datasource, subset, mapping_sets)
    meta_path = _meta_path(cache_dir, datasource, subset, mapping_sets)

    if mode == "date":
        if _consolidate_by_date(datasource, cache_path, subset, mapping_sets, tax_id):
            return cache_path
        msg = (
            f"{datasource!r} has no parseable per-row mapping dates; "
            "falling back to mode='release'."
        )
        logger.warning(msg)
        warnings.warn(msg, UserWarning, stacklevel=2)

    _consolidate_by_release(
        datasource, cache_path, meta_path, subset, mapping_sets, tax_id, show_progress, force
    )
    return cache_path
