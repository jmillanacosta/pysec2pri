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

The cache I/O and release/date walk loop shapes are datasource-agnostic and
live in :mod:`mapkgsutils.consolidate`; this module binds them to
pysec2pri's datasource registries, parsers, and path layout.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any

import mapkgsutils.consolidate as _consolidate

from pysec2pri.constants import ALL_DATASOURCES
from pysec2pri.constants import CONSOLIDATE_DATASOURCES as SUPPORTED_DATASOURCES
from pysec2pri.logging import logger
from pysec2pri.parsers.base import BaseMappingSet, _cmp_versions, product_slug_values

__all__ = [
    "SUPPORTED_DATASOURCES",
    "build_label_history",
    "consolidate_mapping_dates",
    "default_cache_dir",
    "load_mapping_dates",
]


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
    from pysec2pri.downloads import ChEBIDownloader

    config = ALL_DATASOURCES["chebi"]
    subset = kwargs.get("subset") or config.default_subset() or "3star"
    return ChEBIDownloader(subset=subset).list_versions()


def _ensembl_versions(**kwargs: Any) -> list[str]:
    from pysec2pri.downloads import EnsemblDownloader

    return EnsemblDownloader(show_progress=False).list_versions()


def _hgnc_versions(**kwargs: Any) -> list[str]:
    from pysec2pri.downloads import HGNCDownloader

    return HGNCDownloader().list_versions()


def _uniprot_versions(**kwargs: Any) -> list[str]:
    from pysec2pri.downloads import UniProtDownloader

    return UniProtDownloader().list_versions()


# Datasources with a versioned archive that can be walked release by release.
# Datasources absent here (NCBI, VGNC) take the single-parse path instead.
_LIST_VERSIONS_FNS: dict[str, Callable[..., list[str]]] = {
    "chebi": _chebi_versions,
    "ensembl": _ensembl_versions,
    "hgnc": _hgnc_versions,
    "uniprot": _uniprot_versions,
}


def _user_cache_root() -> Path:
    """Return where this OS keeps per-user caches."""
    if sys.platform == "win32":
        local = os.environ.get("LOCALAPPDATA")
        return Path(local) if local else Path.home() / "AppData" / "Local"
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Caches"
    xdg = os.environ.get("XDG_CACHE_HOME")
    return Path(xdg) if xdg else Path.home() / ".cache"


def default_cache_dir() -> Path:
    """Return the default cache directory for consolidated mapping-date indexes.

    Returns:
        ``$PYSEC2PRI_CACHE_DIR`` when set, otherwise a ``pysec2pri`` directory
        under this OS's per-user cache: ``%LOCALAPPDATA%`` on Windows,
        ``~/Library/Caches`` on macOS, ``$XDG_CACHE_HOME`` or ``~/.cache``
        elsewhere.
    """
    env = os.environ.get("PYSEC2PRI_CACHE_DIR")
    return Path(env) if env else _user_cache_root() / "pysec2pri"


def _product_slugs(datasource: str, **kwargs: Any) -> tuple[str, ...]:
    """Return the slug(s) naming which of *datasource*'s datasets this is."""
    config = ALL_DATASOURCES.get(datasource)
    if config is None:
        return ()
    return product_slug_values(config, **kwargs)


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


_sssom_output_path = _consolidate.sssom_output_path
_read_cache = _consolidate.read_cache
_write_cache = _consolidate.write_cache


def _save_optional_output(mapping_set: BaseMappingSet, output: Path | None) -> None:
    """Write *mapping_set* to *output* as SSSOM when an output path is given."""
    if output is not None:
        mapping_set.save("sssom", output)


def load_mapping_dates(
    datasource: str,
    cache_dir: Path | None = None,
    *,
    mapping_sets: str = "ids",
    **kwargs: Any,
) -> dict[str, str]:
    """Load the consolidated ``record_id -> first_seen_date`` index for *datasource*.

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
    return _consolidate.load_mapping_dates(cache_path)


def _parse_mapping_set(
    datasource: str,
    files: dict[str, Path],
    version: str | None,
    mapping_sets: str,
    **kwargs: Any,
) -> Any:
    """Parse one downloaded release into a mapping set, dispatched by datasource."""
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
    inputs: dict[str, Path | str] | None = None,
    **kwargs: Any,
) -> Any:
    """Download one release into a scratch tmpdir, parse it, then clean up.

    The tmpdir is always removed afterwards.
    """
    from pysec2pri.download import download_datasource_with_release

    supplied = {k: Path(v) for k, v in (inputs or {}).items()}
    # Era-aware
    input_map = ALL_DATASOURCES[datasource].inputs_for(mapping_sets, version)
    missing = [key for key in input_map if key not in supplied]

    tmpdir = Path(tempfile.mkdtemp(prefix=f"pysec2pri_consolidate_{datasource}_"))
    try:
        downloaded, version, _ = download_datasource_with_release(
            datasource, tmpdir, version=version, keys=missing, **kwargs
        )
        files = {**supplied, **downloaded}
        return _parse_mapping_set(datasource, files, version, mapping_sets, **kwargs)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _write_consolidated_sssom(
    datasource: str, mapping_sets: str, cache_path: Path, meta_path: Path
) -> tuple[Path, Any]:
    """Build and save the companion SSSOM mapping set next to the cache file.

    Returns ``(output_path, mapping_set)`` -- see
    :func:`mapkgsutils.consolidate.write_consolidated_sssom`.
    """
    from pysec2pri.parsers.base import IdMappingSet, LabelMappingSet

    config = ALL_DATASOURCES[datasource]
    cls = LabelMappingSet if mapping_sets == "labels" else IdMappingSet
    return _consolidate.write_consolidated_sssom(
        cache_path,
        meta_path,
        mapping_set_class=cls,
        record_namespace=str(config.mapping_metadata.get("record_id") or ""),
        mapping_set_metadata=config.mappingset_metadata,
        cardinality_on="label" if mapping_sets == "labels" else "id",
    )


def consolidate_mapping_dates(
    datasource: str,
    *,
    cache_dir: Path | None = None,
    mapping_sets: str = "ids",
    show_progress: bool = True,
    force: bool = False,
    output: Path | None = None,
    inputs: dict[str, Path | str] | None = None,
    **kwargs: Any,
) -> tuple[Path, BaseMappingSet]:
    """Build/update the first-seen-date index for *datasource*.

    Collects every mapping using all available provenance. The walk shape is
    chosen automatically by data availability:

    - Datasources with a versioned archive (e.g., ChEBI, Ensembl, HGNC,
      UniProt): walk every historical release once (oldest first).
      For every mapping seen, records the version/date it first appeared
      and keeps bumping the version/date it was last seen.
    - Datasources without: a single current parse.

    Args:
        datasource: One of :data:`SUPPORTED_DATASOURCES`.
        cache_dir: Directory to write the cache file. Defaults to
            :func:`default_cache_dir`.
        mapping_sets: ``"ids"`` or ``"labels"``.
        show_progress: Whether to show a progress bar over releases
            (versioned-archive datasources only).
        force: Re-scan every release from scratch, ignoring any existing
            cache/resume state (versioned-archive datasources only).
        output: Optional path to also write the consolidated SSSOM mapping
            set to. The internal cache-adjacent copy (see
            :func:`_sssom_output_path`) is always written; when *output* is
            given, the same full mapping set is additionally saved there.
        inputs: Local files (keyed as in the datasource's ``download_urls``)
            to use for the current release instead of downloading. Only
            applies to datasources with no versioned archive (NCBI, VGNC),
            whose consolidate walk is a single current-release parse;
            ignored for the multi-release walk, where *version* is never
            the current release.
        **kwargs: Datasource-specific knobs (``subset`` for ChEBI, ``species``
            for NCBI/Ensembl); ignored for datasources with no such config
            block.

    Returns:
        ``(cache_path, mapping_set)``: the path to the written cache TSV
        (see :func:`_sssom_output_path` for the companion SSSOM mapping set
        written alongside it) and the in-memory consolidated mapping set.

    Raises:
        ValueError: For an unsupported *datasource* or an unsupported
            *mapping_sets* for *datasource*.
    """
    if datasource not in SUPPORTED_DATASOURCES:
        raise ValueError(
            f"Unsupported datasource: {datasource!r}. Supported: {SUPPORTED_DATASOURCES}"
        )
    if mapping_sets not in _SUPPORTED_MAPPING_SETS[datasource]:
        raise ValueError(
            f"{datasource!r} does not support mapping_sets={mapping_sets!r}. "
            f"Supported: {_SUPPORTED_MAPPING_SETS[datasource]}"
        )
    config = ALL_DATASOURCES.get(datasource)
    if config is not None and "species" in config.products:
        species = kwargs.get("species") or config.default_species()
        if species == "all":
            raise ValueError(
                f"{datasource} consolidation requires an explicit single species= taxon ID."
            )

    cache_dir = cache_dir or default_cache_dir()
    cache_path = _cache_path(cache_dir, datasource, mapping_sets, **kwargs)
    meta_path = _meta_path(cache_dir, datasource, mapping_sets, **kwargs)

    list_versions_fn = _LIST_VERSIONS_FNS.get(datasource)
    run_inputs = inputs if list_versions_fn is None else None

    def _run(version: str | None) -> Any:
        return _run_one_version(datasource, version, mapping_sets, inputs=run_inputs, **kwargs)

    if list_versions_fn is None:
        # No versioned archive (NCBI, VGNC): single current parse, full set.
        _consolidate.consolidate(
            cache_path,
            meta_path,
            label=datasource,
            run_one_version=_run,
            show_progress=show_progress,
            force=force,
        )
    else:
        # Versioned archive: walk every historical release for first-seen dates.
        from pysec2pri.download import resolve_release_date

        fn = list_versions_fn
        _consolidate.consolidate(
            cache_path,
            meta_path,
            label=datasource,
            run_one_version=_run,
            list_versions=lambda: fn(**kwargs),
            resolve_release_date=lambda v: resolve_release_date(datasource, v, **kwargs),
            show_progress=show_progress,
            force=force,
        )

    _, mapping_set = _write_consolidated_sssom(datasource, mapping_sets, cache_path, meta_path)
    _save_optional_output(mapping_set, output)
    return cache_path, mapping_set


# Cross-release label history (e.g. Ensembl, whose core schema has no
# previous-gene-symbol table: previous->current symbol transitions
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
) -> BaseMappingSet:
    """Derive previous-label -> current-label transitions across releases.

    Ensembl's core schema has no previous-gene-symbol table, so genuine
    previous -> current symbol transitions are recovered by diffing each
    release's current-label snapshot (``gene``+``xref`` only) against the
    previous release's, oldest to newest. Network-heavy and resumable --
    meant to be run on demand, not as part of normal mapping generation
    (mirrors the release-walk pattern in
    :func:`mapkgsutils.consolidate.consolidate`, but tracks a carried label
    map instead of first/last-seen dates).

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
    if datasource == "ensembl" and species == "all":
        raise ValueError(
            "ensembl label history requires an explicit single species= taxon ID. "
            "Its config default is species='all', which the per-version "
            "download/parse step used here (unlike pysec2pri.api's bulk "
            "generate_ids) has no support for."
        )

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
            files, _, _ = download_datasource_with_release(
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
        # Unlike the consolidate walk's plain-dict cache, these rows become
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


#: Per-(datasource, kind) strategies that recover mappings a source's current
#: release no longer states, from its past releases. Keyed like
#: :data:`_LIST_VERSIONS_FNS`: a source without an entry simply has no extra
#: history to recover for that kind.
_RECOVERY_FNS: dict[tuple[str, str], Callable[..., BaseMappingSet]] = {
    ("ensembl", "labels"): build_label_history,
}


def supports_recovery(datasource: str, kind: str) -> bool:
    """Whether extra mappings can be recovered for *datasource*'s *kind*."""
    return (datasource, kind) in _RECOVERY_FNS


def recover_mapping_set(datasource: str, kind: str, **kwargs: Any) -> BaseMappingSet | None:
    """Recover mappings *datasource*'s current release no longer states.

    Args:
        datasource: Datasource config id.
        kind: Mapping-set kind, e.g. ``"labels"``.
        **kwargs: Forwarded to the strategy; anything it does not declare is
            dropped.

    Returns:
        The recovered mapping set, or ``None`` when this datasource/kind has no
        extra history to recover.
    """
    import inspect

    fn = _RECOVERY_FNS.get((datasource, kind))
    if fn is None:
        return None
    params = inspect.signature(fn).parameters
    accepted = {k: v for k, v in kwargs.items() if k in params and v is not None}
    return fn(datasource, **accepted)
