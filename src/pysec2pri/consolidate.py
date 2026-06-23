"""Build/read the per-mapping first-seen-date index for any versioned datasource.

Some datasources (HGNC, NCBI) carry a real per-row change date in their flat
files (``date_symbol_changed``, ``Discontinue_Date``) -- those are already
wired directly into ``mapping_date`` during normal parsing and don't need
this module at all. Others (ChEBI, UniProt) carry no per-row date anywhere:
a given ``(primary, secondary)`` pair's deprecation date isn't recorded in
any single release.

``Mapping.record_id`` itself is release-scoped (it's the row's OWL Axiom IRI
in SSSOM's RDF/OWL output, and is asserted as part of a specific, versioned
mapping set -- see :meth:`~pysec2pri.parsers.base.BaseParser._record_id`), so
it does *not* match across releases. The version-independent join key this
module needs is :meth:`~pysec2pri.parsers.base.BaseParser._pair_hash` --
the same ``(pri, sec)`` pair always hashes identically regardless of
release, and conveniently is always the trailing 16 hex characters of
``record_id`` (so this module can read it straight off a parsed
``Mapping`` without needing per-parser knowledge of *pri*/*sec*). Walking
every historical release once and recording the first release a pair hash
appears in recovers its true historical date.

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

from pysec2pri.constants import ALL_DATASOURCES
from pysec2pri.logging import logger
from pysec2pri.parsers.base import Sec2PriMappingSet, _cmp_versions

__all__ = [
    "SUPPORTED_DATASOURCES",
    "build_label_history",
    "consolidate_mapping_dates",
    "default_cache_dir",
    "load_mapping_dates",
]

# Datasources this module knows how to download+parse per release.
# NCBI and VGNC have no versioned archive, so only mode="date" applies to them.
SUPPORTED_DATASOURCES = ("chebi", "ensembl", "hgnc", "ncbi", "uniprot", "vgnc")

# mapping_sets kinds each datasource's parser actually supports.
_SUPPORTED_MAPPING_SETS: dict[str, tuple[str, ...]] = {
    "chebi": ("ids", "labels"),
    "ensembl": ("ids", "labels"),
    "hgnc": ("ids", "labels"),
    "ncbi": ("ids", "labels"),
    "uniprot": ("ids",),
    "vgnc": ("ids", "labels"),
}


def _chebi_versions(**kwargs: Any) -> list[str]:
    from pysec2pri.parsers.chebi import ChEBIDownloader

    config = ALL_DATASOURCES["chebi"]
    subset = kwargs.get("subset") or config.default_subset() or "3star"
    return ChEBIDownloader(subset=subset).list_versions()


def _ensembl_versions(**kwargs: Any) -> list[str]:
    from pysec2pri.parsers.ensembl import EnsemblDownloader

    return EnsemblDownloader(show_progress=False).list_versions()


def _hgnc_versions(**kwargs: Any) -> list[str]:
    from pysec2pri.parsers.hgnc import HGNCDownloader

    return HGNCDownloader().list_versions()


def _uniprot_versions(**kwargs: Any) -> list[str]:
    from pysec2pri.parsers.uniprot import UniProtDownloader

    return UniProtDownloader().list_versions()


# Datasources with a versioned archive that can be walked for mode="release".
_LIST_VERSIONS_FNS: dict[str, Callable[..., list[str]]] = {
    "chebi": _chebi_versions,
    "ensembl": _ensembl_versions,
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


def _product_slugs(datasource: str, **kwargs: Any) -> tuple[str, ...]:
    """Return the IRI/path slug(s) disambiguating *this* datasource's product.

    Driven entirely by the datasource's own config -- a ``subset`` block
    (ChEBI) or a ``species`` block (NCBI/Ensembl) -- rather than a hardcoded
    per-datasource knob list, mirroring
    :meth:`~pysec2pri.parsers.base.BaseParser._product_slug`. Datasources
    with neither block get no slug.
    """
    config = ALL_DATASOURCES.get(datasource)
    if config is None:
        return ()
    if config.subset:
        return (str(kwargs.get("subset") or config.default_subset()),)
    if config.species:
        return (str(kwargs.get("species", config.default_species())),)
    return ()


def _cache_dir_for(cache_dir: Path, datasource: str, **kwargs: Any) -> Path:
    """Return ``{cache_dir}/{datasource}/{product_slugs...}/consolidated``."""
    slugs = _product_slugs(datasource, **kwargs)
    return cache_dir.joinpath(datasource, *slugs, "consolidated")


def _cache_path(cache_dir: Path, datasource: str, mapping_sets: str, **kwargs: Any) -> Path:
    base = _cache_dir_for(cache_dir, datasource, **kwargs)
    return base / f"{mapping_sets}_mapping_dates.tsv"


def _meta_path(cache_dir: Path, datasource: str, mapping_sets: str, **kwargs: Any) -> Path:
    base = _cache_dir_for(cache_dir, datasource, **kwargs)
    return base / f"{mapping_sets}_mapping_dates.meta.json"


def _sssom_output_path(cache_path: Path) -> Path:
    """Return the companion SSSOM mapping-set path for a consolidated cache file."""
    return cache_path.with_name(cache_path.stem + "_sssom.tsv")


_CACHE_COLUMNS = (
    "record_id",
    "first_seen_version",
    "first_seen_date",
    "last_seen_version",
    "last_seen_date",
    # JSON-encoded snapshot of the mapping's own fields (subject_id,
    # object_id, predicate_id, ...) as last seen -- used to materialize the
    # consolidated index as a real SSSOM mapping set (see
    # _build_consolidated_mapping_set). Empty for legacy/hand-built rows.
    "fields_json",
)


def _read_cache(cache_path: Path) -> dict[str, dict[str, str]]:
    """Read a consolidated mapping-date cache TSV into a dict keyed by record_id."""
    if not cache_path.exists():
        return {}

    import polars as pl

    df = pl.read_csv(cache_path, separator="\t", schema_overrides={"record_id": pl.Utf8})
    cols = [c for c in _CACHE_COLUMNS[1:] if c in df.columns]
    return {
        str(row["record_id"]): {col: str(row[col]) for col in cols}
        for row in df.iter_rows(named=True)
    }


def _write_cache(cache_path: Path, records: dict[str, dict[str, str]]) -> None:
    """Write the merged ``record_id -> first/last seen`` dict to a TSV file."""
    import polars as pl

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {"record_id": rid, **{c: fields.get(c, "") for c in _CACHE_COLUMNS[1:]}}
        for rid, fields in records.items()
    ]
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
    mapping_sets: str = "ids",
    **kwargs: Any,
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
        mapping_sets: ``"ids"`` or ``"labels"`` -- must match the
            mapping-set kind used when the index was built.
        **kwargs: Datasource-specific knobs (``subset`` for ChEBI, ``species``
            for NCBI/Ensembl) -- must match what was used when the index was
            built; ignored for datasources with no such config block.

    Returns:
        Dict mapping each ``record_id`` to its first-seen ISO date string.
    """
    cache_path = _cache_path(cache_dir or default_cache_dir(), datasource, mapping_sets, **kwargs)
    records = _read_cache(cache_path)
    return {rid: fields["first_seen_date"] for rid, fields in records.items()}


def _parse_mapping_set(
    datasource: str,
    files: dict[str, Path],
    version: str | None,
    mapping_sets: str,
    **kwargs: Any,
) -> Any:
    """Parse one downloaded release into a mapping set, dispatched by datasource.

    Calls each parser with explicit paths from *files* rather than relying
    on directory auto-discovery, so consolidation never depends on the
    tmpdir-guessing logic used by the ``generate_*`` convenience functions.
    Each branch pulls only the kwarg it needs (falling back to the
    datasource's own config default), so unsupported kwargs are just ignored.
    """
    config = ALL_DATASOURCES.get(datasource)
    if datasource == "chebi":
        subset = kwargs.get("subset") or (config.default_subset() if config else None) or "3star"
        return _parse_chebi_version(files, version, subset, mapping_sets)
    if datasource == "hgnc":
        return _parse_hgnc_version(files, version, mapping_sets)
    if datasource == "ncbi":
        species = kwargs.get("species", config.default_species() if config else "9606")
        return _parse_ncbi_version(files, version, mapping_sets, species)
    if datasource == "ensembl":
        species = kwargs.get("species", config.default_species() if config else "9606")
        return _parse_ensembl_version(files, version, mapping_sets, species)
    if datasource == "uniprot":
        return _parse_uniprot_version(files, version, mapping_sets)
    if datasource == "vgnc":
        species = kwargs.get("species", config.default_species() if config else "9598")
        return _parse_vgnc_version(files, version, mapping_sets, species)
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
    files: dict[str, Path], version: str | None, mapping_sets: str, species: str
) -> Any:
    from pysec2pri.parsers.ncbi import NCBIParser

    parser = NCBIParser(version=version, show_progress=False)
    if mapping_sets == "ids":
        return parser.parse(
            files["gene_history"], species=species, gene_info_path=files["gene_info"]
        )
    return parser.parse_labels(files["gene_info"], species=species)


def _parse_ensembl_version(
    files: dict[str, Path], version: str | None, mapping_sets: str, species: str
) -> Any:
    from pysec2pri.parsers.ensembl import EnsemblParser

    parser = EnsemblParser(version=version, show_progress=False, species=species)
    if mapping_sets == "ids":
        return parser.parse(
            files["stable_id_event"],
            mapping_session_path=files.get("mapping_session"),
            gene_path=files.get("gene"),
        )
    return parser.parse_labels(files.get("gene"), files.get("xref"), files.get("external_synonym"))


def _parse_uniprot_version(files: dict[str, Path], version: str | None, mapping_sets: str) -> Any:
    from pysec2pri.parsers.uniprot import UniProtParser

    if mapping_sets != "ids":
        raise ValueError("uniprot only supports mapping_sets='ids'")
    parser = UniProtParser(version=version, show_progress=False)
    return parser.parse(files.get("sec_ac"), delac_path=files.get("delac_sp"))


def _parse_vgnc_version(
    files: dict[str, Path], version: str | None, mapping_sets: str, species: str
) -> Any:
    from pysec2pri.parsers.vgnc import VGNCParser

    parser = VGNCParser(version=version, show_progress=False)
    if mapping_sets == "ids":
        return parser.parse(files["withdrawn"], complete_set_path=files.get("complete"))
    return parser.parse_labels(files["complete"], species=species)


def _run_one_version(
    datasource: str,
    version: str | None,
    mapping_sets: str,
    **kwargs: Any,
) -> Any:
    """Download one release into a scratch tmpdir, parse it, then clean up.

    The tmpdir is always removed afterwards: walking ~250 releases must not
    accumulate downloaded files on disk.
    """
    from pysec2pri.download import download_datasource_with_release

    tmpdir = Path(tempfile.mkdtemp(prefix=f"pysec2pri_consolidate_{datasource}_"))
    try:
        files, _ = download_datasource_with_release(datasource, tmpdir, version=version, **kwargs)
        return _parse_mapping_set(datasource, files, version, mapping_sets, **kwargs)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _mapping_fields_json(m: Any) -> str:
    """JSON-encode a mapping's own fields (excluding mapping_date/record_id).

    Snapshots a mapping's current shape (subject_id, object_id,
    predicate_id, ...) so the consolidated index can later be materialized
    as a real SSSOM mapping set (see :func:`_build_consolidated_mapping_set`).
    ``default=str`` handles enum-valued fields like ``mapping_cardinality``,
    which round-trip cleanly back through ``Mapping(**fields)``. Returns
    ``"{}"`` for non-dataclass stand-ins (e.g. test doubles) that don't carry
    a full field set -- callers skip materializing those rows.
    """
    from dataclasses import fields as dataclass_fields
    from dataclasses import is_dataclass

    if not is_dataclass(m):
        return "{}"
    fields = {
        f.name: getattr(m, f.name)
        for f in dataclass_fields(m)
        if f.name not in ("mapping_date", "record_id") and getattr(m, f.name, None) is not None
    }
    return json.dumps(fields, default=str)


def _consolidate_by_date(
    datasource: str,
    cache_path: Path,
    meta_path: Path,
    mapping_sets: str,
    **kwargs: Any,
) -> bool:
    """Single-pass "date" mode: capture each row's own real ``mapping_date``.

    Returns ``True`` and writes the cache when *datasource* produced at
    least one real per-row date; returns ``False`` (cache untouched) when it
    produced none, so the caller can fall back to ``"release"`` mode.
    """
    mapping_set = _run_one_version(datasource, None, mapping_sets, **kwargs)
    dated = [m for m in (mapping_set.mappings or []) if getattr(m, "mapping_date", None)]
    if not dated:
        return False

    version_label = str(getattr(mapping_set, "mapping_set_version", None) or "current")
    records: dict[str, dict[str, str]] = {}
    for m in dated:
        # record_id is now release-scoped (see BaseParser._record_id); its
        # trailing 16 hex chars are always the version-independent pair hash
        # (BaseParser._pair_hash), which is what must match across releases.
        pair_key = str(getattr(m, "record_id", None) or "")[-16:]
        if not pair_key:
            continue
        date_str = str(m.mapping_date)
        records[pair_key] = {
            "first_seen_version": version_label,
            "first_seen_date": date_str,
            "last_seen_version": version_label,
            "last_seen_date": date_str,
            "fields_json": _mapping_fields_json(m),
        }
    _write_cache(cache_path, records)
    _write_meta(meta_path, version_label)
    return True


def _build_consolidated_mapping_set(
    datasource: str,
    mapping_sets: str,
    records: dict[str, dict[str, str]],
    last_version: str | None,
) -> Any:
    """Materialize the consolidated index as a real SSSOM mapping set.

    Each record's stored field snapshot (from the most recent release where
    it was seen) is rebuilt into a ``Mapping``, with ``mapping_date``
    overridden to its ``first_seen_date`` -- the date of appearance, rather
    than whatever date the snapshot's own release happened to carry -- and
    ``record_id`` re-scoped with a trailing ``consolidate`` segment (mirrors
    :meth:`~pysec2pri.parsers.ensembl.EnsemblParser.parse_label_history`),
    marking the IRI as a derived, cross-release product. *last_version* is
    the most recent release the walk has processed.
    """
    from sssom_schema import Mapping

    from pysec2pri.parsers.base import IdMappingSet, LabelMappingSet

    config = ALL_DATASOURCES[datasource]
    ms_meta = config.mappingset_metadata
    base_ns = str(config.mapping_metadata.get("record_id") or "")
    version_label = str(last_version) if last_version else "current"

    mappings = []
    for pair_key, fields in records.items():
        fields_json = fields.get("fields_json") or ""
        if not fields_json:
            continue
        try:
            row_fields = json.loads(fields_json)
        except json.JSONDecodeError:
            continue
        if not row_fields:
            # Non-dataclass stand-ins (see _mapping_fields_json) have
            # nothing to materialize into a real Mapping; date bookkeeping
            # for them still lives in the TSV cache, just not in the SSSOM.
            continue
        row_fields["mapping_date"] = fields.get("first_seen_date")
        row_fields["record_id"] = f"{base_ns}{version_label}/consolidate/{pair_key}"
        mappings.append(Mapping(**row_fields))

    cls = LabelMappingSet if mapping_sets == "labels" else IdMappingSet
    base_ms_id = str(ms_meta.get("mapping_set_id") or "")
    mapping_set = cls(
        mappings=mappings,
        curie_map=ms_meta.get("curie_map") or {},
        mapping_set_id=f"{base_ms_id}/{version_label}/consolidate",
        mapping_set_version=version_label,
        mapping_set_title=ms_meta.get("mapping_set_title"),
        mapping_set_description=ms_meta.get("mapping_set_description"),
        license=ms_meta.get("license"),
    )
    mapping_set.compute_cardinalities()
    return mapping_set


def _write_consolidated_sssom(
    datasource: str, mapping_sets: str, cache_path: Path, meta_path: Path
) -> Path:
    """Build and save the companion SSSOM mapping set next to the cache file."""
    records = _read_cache(cache_path)
    last_version = _read_meta(meta_path)
    mapping_set = _build_consolidated_mapping_set(datasource, mapping_sets, records, last_version)
    output_path = _sssom_output_path(cache_path)
    mapping_set.save("sssom", output_path)
    return output_path


def _consolidate_by_release(
    datasource: str,
    cache_path: Path,
    meta_path: Path,
    mapping_sets: str,
    show_progress: bool,
    force: bool,
    **kwargs: Any,
) -> None:
    """Historical-walk "release" mode: track first/last-seen release per mapping."""
    from pysec2pri.download import resolve_release_date

    list_versions_fn = _LIST_VERSIONS_FNS.get(datasource)
    if list_versions_fn is None:
        raise ValueError(f"{datasource!r} has no versioned archive; only mode='date' is supported.")

    records: dict[str, dict[str, str]] = {} if force else _read_cache(cache_path)
    last_version = None if force else _read_meta(meta_path)

    versions = list_versions_fn(**kwargs)
    if last_version is not None:
        versions = [v for v in versions if _cmp_versions(v, last_version) > 0]

    iterator: Iterable[str] = versions
    if show_progress:
        from tqdm import tqdm

        iterator = tqdm(versions, desc=f"Consolidating {datasource.upper()} mapping dates")

    for v in iterator:
        try:
            mapping_set = _run_one_version(datasource, v, mapping_sets, **kwargs)
            release_date = resolve_release_date(datasource, v, **kwargs)
            date_str = release_date.date().isoformat() if release_date else v

            for m in mapping_set.mappings or []:
                # record_id is release-scoped; match across releases on its
                # trailing pair hash instead (see BaseParser._pair_hash).
                pair_key = str(getattr(m, "record_id", None) or "")[-16:]
                if not pair_key:
                    continue
                entry = records.setdefault(
                    pair_key,
                    {"first_seen_version": v, "first_seen_date": date_str},
                )
                entry["last_seen_version"] = v
                entry["last_seen_date"] = date_str
                entry["fields_json"] = _mapping_fields_json(m)
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
    mapping_sets: str = "ids",
    show_progress: bool = True,
    force: bool = False,
    **kwargs: Any,
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

    Either way, alongside the internal cache file this also writes a
    companion real SSSOM mapping set (see :func:`_sssom_output_path`) where
    every row's ``mapping_date`` is its true first-seen date, rather than
    whichever release happened to be parsed last.

    Args:
        datasource: One of :data:`SUPPORTED_DATASOURCES`.
        mode: ``"release"`` or ``"date"``.
        cache_dir: Directory to write the cache file. Defaults to
            :func:`default_cache_dir`.
        mapping_sets: ``"ids"`` or ``"labels"``.
        show_progress: Whether to show a progress bar over releases
            (``"release"`` mode only).
        force: Re-scan every release from scratch, ignoring any existing
            cache/resume state (``"release"`` mode only).
        **kwargs: Datasource-specific knobs (``subset`` for ChEBI, ``species``
            for NCBI/Ensembl); ignored for datasources with no such config
            block.

    Returns:
        Path to the written cache TSV (see :func:`_sssom_output_path` for
        the companion SSSOM mapping set written alongside it).

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
    cache_path = _cache_path(cache_dir, datasource, mapping_sets, **kwargs)
    meta_path = _meta_path(cache_dir, datasource, mapping_sets, **kwargs)

    if mode == "date":
        if _consolidate_by_date(datasource, cache_path, meta_path, mapping_sets, **kwargs):
            _write_consolidated_sssom(datasource, mapping_sets, cache_path, meta_path)
            return cache_path
        msg = (
            f"{datasource!r} has no parseable per-row mapping dates; "
            "falling back to mode='release'."
        )
        logger.warning(msg)
        warnings.warn(msg, UserWarning, stacklevel=2)

    _consolidate_by_release(
        datasource, cache_path, meta_path, mapping_sets, show_progress, force, **kwargs
    )
    _write_consolidated_sssom(datasource, mapping_sets, cache_path, meta_path)
    return cache_path


# Cross-release label history (e.g. Ensembl, whose core schema has no
# previous-gene-symbol table -- genuine previous->current symbol transitions
# are recovered by diffing each release's current-label snapshot).


def _label_transitions(
    prev_map: dict[str, str], curr_map: dict[str, str]
) -> list[tuple[str, str, str]]:
    """Return ``(stable_id, prev_label, curr_label)`` for every changed gene.

    Compares two ``{stable_id -> label}`` snapshots (see
    :meth:`~pysec2pri.parsers.ensembl.EnsemblParser.current_label_snapshot`)
    and yields one entry per gene whose label differs between them. Genes
    present in only one snapshot, or with an unchanged label, are skipped --
    a pure function so the diff logic is testable without any I/O.

    Args:
        prev_map: ``{stable_id -> label}`` snapshot from the earlier release.
        curr_map: ``{stable_id -> label}`` snapshot from the later release.

    Returns:
        List of ``(stable_id, prev_label, curr_label)`` tuples.
    """
    transitions: list[tuple[str, str, str]] = []
    for stable_id, curr_label in curr_map.items():
        prev_label = prev_map.get(stable_id)
        if prev_label is not None and prev_label != curr_label:
            transitions.append((stable_id, prev_label, curr_label))
    return transitions


def _label_history_dir(cache_dir: Path, datasource: str, species: str) -> Path:
    return cache_dir.joinpath(datasource, str(species), "consolidated")


def _read_label_history_state(state_path: Path) -> dict[str, Any]:
    if not state_path.exists():
        return {}
    try:
        data: dict[str, Any] = json.loads(state_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return data


def _write_label_history_state(state_path: Path, state: dict[str, Any]) -> None:
    state_path.parent.mkdir(parents=True, exist_ok=True)
    state_path.write_text(json.dumps(state), encoding="utf-8")


def build_label_history(
    datasource: str = "ensembl",
    *,
    species: str | int = 9606,
    from_version: str | None = None,
    to_version: str | None = None,
    cache_dir: Path | None = None,
    show_progress: bool = True,
    force: bool = False,
) -> Sec2PriMappingSet:
    """Derive previous-label -> current-label transitions across releases.

    Ensembl's core schema has no previous-gene-symbol table, so genuine
    previous -> current symbol transitions are recovered by diffing each
    release's current-label snapshot (``gene``+``xref`` only) against the
    previous release's, oldest to newest. Network-heavy and resumable --
    meant to be run on demand, not as part of normal mapping generation
    (mirrors the release-walk pattern in :func:`_consolidate_by_release`,
    but tracks a carried label map instead of first/last-seen dates).

    Args:
        datasource: Currently only ``"ensembl"`` is supported.
        species: Canonical NCBI taxon ID.
        from_version: Optional lower bound (inclusive) on the release walk.
        to_version: Optional upper bound (inclusive) on the release walk.
        cache_dir: Directory for the resumable state + output. Defaults to
            :func:`default_cache_dir`.
        show_progress: Whether to show a progress bar over releases.
        force: Re-walk every release from scratch, ignoring resume state.

    Returns:
        :class:`~pysec2pri.parsers.base.LabelMappingSet` of
        previous -> current symbol transitions.

    Raises:
        ValueError: If *datasource* has no configured version-list function.
    """
    list_versions_fn = _LIST_VERSIONS_FNS.get(datasource)
    if list_versions_fn is None:
        raise ValueError(f"build_label_history does not support {datasource!r}")

    from pysec2pri.download import download_datasource_with_release, resolve_release_date
    from pysec2pri.parsers.ensembl import EnsemblParser

    cache_dir = cache_dir or default_cache_dir()
    out_dir = _label_history_dir(cache_dir, datasource, str(species))
    state_path = out_dir / "label_history_state.json"
    output_path = out_dir / "label_history_sssom.tsv"

    versions = list_versions_fn()
    if from_version is not None:
        versions = [v for v in versions if _cmp_versions(v, from_version) >= 0]
    if to_version is not None:
        versions = [v for v in versions if _cmp_versions(v, to_version) <= 0]

    state = {} if force else _read_label_history_state(state_path)
    label_map: dict[str, str] = state.get("label_map", {})
    transition_rows: list[tuple[str, str, str, str | None]] = state.get("transitions", [])
    last_version: str | None = state.get("last_version")

    pending = versions
    if last_version is not None:
        pending = [v for v in versions if _cmp_versions(v, last_version) > 0]

    iterator: Iterable[str] = pending
    if show_progress:
        from tqdm import tqdm

        iterator = tqdm(pending, desc=f"Walking {datasource.upper()} releases for label history")

    for v in iterator:
        tmpdir = Path(tempfile.mkdtemp(prefix=f"pysec2pri_labelhistory_{datasource}_"))
        try:
            files, _ = download_datasource_with_release(
                datasource, tmpdir, version=v, species=species, keys=["gene", "xref"]
            )
            parser = EnsemblParser(version=v, show_progress=False, species=species)
            curr_map = parser.current_label_snapshot(files["gene"], files["xref"])
        except Exception:
            logger.warning(
                "Skipping %s version %s during label-history walk", datasource, v, exc_info=True
            )
            continue
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

        release_date = resolve_release_date(datasource, v, species=species)
        # Unlike _consolidate_by_release's plain-dict cache, these rows become
        # real Mapping objects -- mapping_date must be a valid date or None,
        # never the raw (non-date-shaped) version string.
        date_str = release_date.date().isoformat() if release_date else None

        for stable_id, prev_label, curr_label in _label_transitions(label_map, curr_map):
            transition_rows.append((stable_id, prev_label, curr_label, date_str))

        label_map = curr_map
        last_version = v
        _write_label_history_state(
            state_path,
            {"label_map": label_map, "last_version": last_version, "transitions": transition_rows},
        )

        parser = EnsemblParser(version=last_version, show_progress=False, species=species)
        mapping_set = parser.parse_label_history(
            (rid, prev, curr, date) for rid, prev, curr, date in transition_rows
        )
        mapping_set.save("sssom", output_path)

    parser = EnsemblParser(version=last_version, show_progress=False, species=species)
    return parser.parse_label_history(
        (rid, prev, curr, date) for rid, prev, curr, date in transition_rows
    )
