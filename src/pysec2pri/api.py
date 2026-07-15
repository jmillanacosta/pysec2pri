"""Main functions for pysec2pri.

This module provides functions for parsing biological database
secondary-to-primary mapping files and generating and using the standardized Mapping sets.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

from mapkgsutils.context import ContextSpec, load_xref_mapping
from mapkgsutils.parsers.config import get_datasource_config

from pysec2pri.exports import (
    write_json,
    write_label_sec2pri,
    write_name2synonym,
    write_output,
    write_owl,
    write_pri_ids,
    write_rdf,
    write_sec2pri,
    write_sssom,
)

if TYPE_CHECKING:
    from datetime import datetime

    import pandas as pd
    from mapkgsutils.context import XrefMapping
    from mapkgsutils.diff import MappingDiff

    from pysec2pri.parsers.base import BaseMappingSet, BaseParser

from pysec2pri.parsers.base import (
    ALL_SPECIES,
    AmbiguousMappingSet,
    IdMappingSet,
    LabelMappingSet,
)

_CONFIG_PACKAGE = "pysec2pri.config"

#: Sources that publish one dataset per species instead of one per release, so
#: ``species="all"`` must download and combine each species in turn rather than
#: filter a single file.
_BULK_SPECIES_SOURCES = ("ensembl",)

__all__ = [
    "ContextSpec",
    "combine_mapping_sets",
    "find_ambiguous",
    "generate_ids",
    "generate_labels",
    "generate_primary_ids",
    "generate_primary_labels",
    "list_versions",
    "load_label_mapping",
    "load_mapping",
    "load_xref_mapping",
    "resolve_ids",
    "resolve_labels",
    "save",
    "sources",
    "supports_consolidate",
    "write_all_formats",
    "write_diff_output",
    "write_json",
    "write_label_sec2pri",
    "write_name2synonym",
    "write_output",
    "write_owl",
    "write_rdf",
    "write_sec2pri",
    "write_sssom",
]


def _auto_download(
    datasource: str,
    version: str | None = None,
    keys: list[str] | None = None,
    show_progress: bool = True,
    **options: Any,
) -> tuple[dict[str, Path], datetime | None]:
    """Download files for *datasource* into a temp dir.

    Args:
        datasource: Datasource name (e.g. ``"hgnc"``).
        version: Optional specific version to download.
        keys: Optional list of file-key names to download. When given, only
            those keys are fetched (e.g. ``["complete"]``).
        show_progress: Whether to show download/decompression progress bars.
        **options: Forwarded to the datasource's downloader. Some sources pick
            their files by more than version (e.g. a species selector).

    Returns:
        Tuple of (file-key -> downloaded path mapping, source release date or
        None). The release date is set on the parser so the generated mapping
        set's ``mapping_date`` reflects the upstream release.
    """
    import tempfile

    from pysec2pri.download import download_datasource_with_release

    tmpdir = Path(tempfile.mkdtemp(prefix=f"pysec2pri_{datasource}_"))
    return download_datasource_with_release(
        datasource, tmpdir, version=version, keys=keys, show_progress=show_progress, **options
    )


def _accepted(fn: Any, **candidates: Any) -> dict[str, Any]:
    """Return the subset of *candidates* that *fn* declares and that are set.

    A callable's signature is the source of truth for which options apply to
    it, so callers can offer every option without knowing which datasource
    takes ``subset``, ``species`` or ``entity_type``.
    """
    import inspect

    params = inspect.signature(fn).parameters
    return {k: v for k, v in candidates.items() if k in params and v is not None}


def _resolve_parser_class(config_id: str) -> type[BaseParser]:
    """Return the parser class named by a datasource config's ``parser_class``.

    Args:
        config_id: Datasource config id, e.g. ``"hgnc"``.

    Returns:
        The parser class declared in ``config/<config_id>.yaml``.

    Raises:
        ValueError: If the declared class cannot be found.
    """
    import importlib

    cfg = get_datasource_config(config_id, config_package=_CONFIG_PACKAGE)
    name = cfg.parser_class
    parsers = importlib.import_module("pysec2pri.parsers")
    cls = getattr(parsers, name, None)
    if cls is None:  # parsers that are not re-exported at package level
        hmdb = importlib.import_module("pysec2pri.parsers.hmdb")
        cls = getattr(hmdb, name, None)
    if cls is None:
        raise ValueError(f"{config_id}: unknown parser_class {name!r}")
    return cast("type[BaseParser]", cls)


def sources(kind: str | None = None) -> list[str]:
    """Return the datasources the config files declare.

    Args:
        kind: Restrict to sources declaring this mapping-set kind, e.g.
            ``"labels"``. ``None`` returns every source.

    Returns:
        Sorted datasource names accepted by :func:`generate_ids` and friends.
    """
    from pysec2pri.parsers.base import CONFIG_DIR

    names = []
    for path in sorted(CONFIG_DIR.glob("*.yaml")):
        cfg = get_datasource_config(path.stem, config_package=_CONFIG_PACKAGE)
        if kind is None or kind in cfg.mapping_sets:
            names.append(path.stem)
    return names


def supports_consolidate(source: str, kind: str = "ids") -> bool:
    """Whether *source* can recover extra history for *kind* (see ``--consolidate``)."""
    cfg = get_datasource_config(source, config_package=_CONFIG_PACKAGE)
    return bool((cfg.mapping_sets.get(kind) or {}).get("consolidate"))


def _refresh_consolidated(
    datasource: str,
    kind: str,
    *,
    cache_dir: Path | None = None,
    force: bool = False,
    **options: Any,
) -> None:
    """Build/refresh the cross-release index the parser reads mapping dates from."""
    from pysec2pri.consolidate import consolidate_mapping_dates

    consolidate_mapping_dates(
        datasource,
        cache_dir=cache_dir,
        mapping_sets=kind,
        force=force,
        **_accepted(consolidate_mapping_dates, **options),
    )


def _process_one_ensembl_species(
    kind: str,
    version: str,
    token: str,
    assembly: str,
    release_date: datetime | None,
) -> BaseMappingSet:
    """Download and parse one Ensembl species, for the ``species="all"`` bulk path.

    Builds download URLs directly from *token*/*assembly* (as discovered by
    :func:`~pysec2pri.downloads.ensembl.discover_ensembl_species`), bypassing
    the taxon-ID-based resolution :class:`EnsemblDownloader` normally does
    -- the bulk path already knows the exact token to use for every species
    in one shot, so re-resolving each one individually would be redundant.
    """
    import shutil
    import tempfile

    from mapkgsutils.download import download_urls

    from pysec2pri.constants import ENSEMBL
    from pysec2pri.parsers.ensembl import EnsemblParser

    keys = (
        ["stable_id_event", "mapping_session", "gene"]
        if kind == "ids"
        else ["gene", "xref", "external_synonym"]
    )
    urls = {
        key: ENSEMBL.download_urls[key].format(version=version, species=token, assembly=assembly)
        for key in keys
    }
    tmpdir = Path(tempfile.mkdtemp(prefix=f"pysec2pri_ensembl_all_{token}_"))
    try:
        files = download_urls(urls, tmpdir, decompress=True)
        parser = EnsemblParser(version=version, show_progress=False, species=token)
        parser.release_date = release_date
        if kind == "ids":
            return parser.parse(
                files["stable_id_event"],
                mapping_session_path=files.get("mapping_session"),
                gene_path=files.get("gene"),
            )
        return parser.parse_labels(
            files.get("gene"), files.get("xref"), files.get("external_synonym")
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _generate_ensembl_all_species(
    kind: str,
    version: str | None,
    show_progress: bool,
) -> BaseMappingSet:
    """Process every Ensembl species at *version* and combine into one mapping set.

    Args:
        kind: ``"ids"`` or ``"labels"``.
        version: Ensembl release number. Latest release is used when
            ``None``.
        show_progress: Whether to show a progress bar over species.

    Returns:
        Combined :class:`~pysec2pri.parsers.base.BaseMappingSet`.

    Raises:
        ValueError: If no species could be processed at all.
    """
    from pysec2pri.download import check_ensembl_release, resolve_release_date
    from pysec2pri.downloads.ensembl import discover_ensembl_species
    from pysec2pri.logging import logger
    from pysec2pri.parsers.ensembl import ALL_SPECIES, EnsemblParser

    if version is None:
        version = check_ensembl_release().version
        if version is None:
            raise ValueError("Could not determine the latest Ensembl release.")
    release_date = resolve_release_date("ensembl", version, species=9606)

    species_list = discover_ensembl_species(version)
    iterator: Any = species_list
    if show_progress:
        from tqdm import tqdm

        iterator = tqdm(species_list, desc=f"Processing {len(species_list)} Ensembl species")

    all_mappings: list[Any] = []
    primary_ids: set[str] = set()
    primary_labels: dict[str, set[str]] = {}
    processed = 0
    for token, assembly in iterator:
        try:
            ms = _process_one_ensembl_species(kind, version, token, assembly, release_date)
        except Exception:
            logger.warning(
                "Skipping Ensembl species %s during all-species run", token, exc_info=True
            )
            continue
        all_mappings.extend(ms.mappings or [])
        primary_ids |= getattr(ms, "_primary_ids", None) or set()
        for label, ids in (getattr(ms, "_primary_labels", None) or {}).items():
            primary_labels.setdefault(label, set()).update(ids)
        processed += 1

    if processed == 0:
        raise ValueError(f"No Ensembl species could be processed for release {version!r}.")

    combined_parser = EnsemblParser(version=version, show_progress=False, species=ALL_SPECIES)
    combined_parser.release_date = release_date
    combined = combined_parser.create_mapping_set(
        all_mappings, mapping_type="id" if kind == "ids" else "label"
    )
    if primary_ids:
        object.__setattr__(combined, "_primary_ids", primary_ids)
    if primary_labels:
        object.__setattr__(combined, "_primary_labels", primary_labels)
    return combined


def _generate(
    datasource: str,
    kind: str = "ids",
    *,
    version: str | None = None,
    show_progress: bool = True,
    inputs: dict[str, Path | str] | None = None,
    consolidate: bool = False,
    cache_dir: Path | None = None,
    force: bool = False,
    **options: Any,
) -> BaseMappingSet:
    """Build a mapping set from a datasource's config, with no per-source code.

    The datasource's ``config/<datasource>.yaml`` declares everything needed:
    ``parser_class``, and per ``mapping_sets`` kind a ``method`` plus an
    ``inputs`` map binding each ``download_urls`` key to a parameter of that
    method. Any file not supplied via *inputs* is downloaded. Options are
    offered to both the parser constructor and *method*, each taking only what
    its signature declares.

    Args:
        datasource: Datasource config id, e.g. ``"hgnc"``.
        kind: Mapping-set kind declared in the config, e.g. ``"ids"``.
        version: Version string for metadata and downloads.
        show_progress: Whether to show progress bars.
        inputs: Local paths keyed by ``download_urls`` key, to use instead of
            downloading.
        **options: Extra options such as ``species``/``subset``/
            ``entity_type``; ignored where a datasource does not accept them.

    Returns:
        The generated mapping set.

    Raises:
        ValueError: If *datasource* or *kind* is not declared.
    """
    if datasource not in sources():
        raise ValueError(f"Unknown datasource {datasource!r}. Available: {sources()}")
    cfg = get_datasource_config(datasource, config_package=_CONFIG_PACKAGE)
    spec = cfg.mapping_sets.get(kind)
    if spec is None:
        raise ValueError(
            f"{datasource!r} declares no {kind!r} mapping set. "
            f"Available: {sorted(cfg.mapping_sets)}"
        )

    if consolidate:
        if not spec.get("consolidate"):
            raise ValueError(
                f"{datasource!r} {kind!r} cannot be consolidated: its releases carry no "
                "extra history to recover. See `supports_consolidate`."
            )
        _refresh_consolidated(datasource, kind, cache_dir=cache_dir, force=force, **options)

    # A source may publish one dataset per species rather than one per release;
    # "every species" then means downloading and combining each in turn.
    if (
        not inputs
        and options.get("species") == ALL_SPECIES
        and datasource in _BULK_SPECIES_SOURCES
        and kind in ("ids", "labels")
    ):
        return _generate_ensembl_all_species(kind, version, show_progress)

    input_map: dict[str, str] = spec.get("inputs") or {}
    supplied = {k: Path(v) for k, v in (inputs or {}).items()}

    release_date = None
    if missing := [key for key in input_map if key not in supplied]:
        files, release_date = _auto_download(
            datasource, version, keys=missing, show_progress=show_progress, **options
        )
        supplied.update({k: Path(v) for k, v in files.items()})

    parser_cls = _resolve_parser_class(datasource)
    parser = parser_cls(
        **_accepted(parser_cls.__init__, version=version, show_progress=show_progress, **options)
    )
    parser.release_date = release_date

    method = getattr(parser, spec["method"])
    call_kwargs = {param: supplied[key] for key, param in input_map.items() if key in supplied}
    call_kwargs.update(_accepted(method, **options))
    result: BaseMappingSet = method(**call_kwargs)

    if consolidate:
        from pysec2pri.consolidate import recover_mapping_set

        recovered = recover_mapping_set(
            datasource,
            kind,
            cache_dir=cache_dir,
            force=force,
            show_progress=show_progress,
            **options,
        )
        if recovered is not None:
            result = combine_mapping_sets(result, recovered)
    return result


def generate_ids(
    source: str,
    *,
    version: str | None = None,
    show_progress: bool = True,
    inputs: dict[str, Path | str] | None = None,
    consolidate: bool = False,
    cache_dir: Path | None = None,
    force: bool = False,
    **options: Any,
) -> BaseMappingSet:
    """Return *source*'s secondary-to-primary **ID** mappings.

    Args:
        source: Datasource name; see :func:`sources` for what is available.
        version: Release to build. The latest is used when ``None``.
        show_progress: Whether to show progress bars.
        inputs: Local input files keyed as in the config's ``download_urls``
            (e.g. ``{"withdrawn": "withdrawn.txt"}``). Anything omitted is
            downloaded.
        consolidate: Recover mappings the current release's files no longer
            state, by walking the source's historical releases, and stamp every
            mapping with the release it first appeared in. Slow and
            network-heavy; resumable via *cache_dir*. Only for sources whose
            releases carry such history (see :func:`supports_consolidate`).
        cache_dir: Where to keep the resumable cross-release index. Defaults to
            ``$PYSEC2PRI_CACHE_DIR`` or ``~/.cache/pysec2pri``.
        force: Re-walk every release, ignoring any resume state.
        **options: Source-specific options, e.g. ``species`` for NCBI/VGNC/
            Ensembl, ``subset`` for ChEBI, ``entity_type`` for Wikidata. An
            option a source does not accept is ignored.

    Returns:
        An ``IdMappingSet`` of secondary -> primary identifier mappings.
    """
    return _generate(
        source,
        "ids",
        version=version,
        show_progress=show_progress,
        inputs=inputs,
        consolidate=consolidate,
        cache_dir=cache_dir,
        force=force,
        **options,
    )


def generate_labels(
    source: str,
    *,
    version: str | None = None,
    show_progress: bool = True,
    inputs: dict[str, Path | str] | None = None,
    consolidate: bool = False,
    cache_dir: Path | None = None,
    force: bool = False,
    **options: Any,
) -> BaseMappingSet:
    """Return *source*'s previous/alias-to-current **label** mappings.

    Args:
        source: Datasource name; see ``sources("labels")`` for what is
            available.
        version: Release to build. The latest is used when ``None``.
        show_progress: Whether to show progress bars.
        inputs: Local input files keyed as in the config's ``download_urls``.
            Anything omitted is downloaded.
        consolidate: Recover label changes the current release's files no
            longer state, by walking historical releases, and stamp each with
            the release it first appeared in. See :func:`generate_ids`.
        cache_dir: Where to keep the resumable cross-release index.
        force: Re-walk every release, ignoring any resume state.
        **options: Source-specific options; see :func:`generate_ids`.

    Returns:
        A ``LabelMappingSet`` of secondary -> primary label mappings.
    """
    return _generate(
        source,
        "labels",
        version=version,
        show_progress=show_progress,
        inputs=inputs,
        consolidate=consolidate,
        cache_dir=cache_dir,
        force=force,
        **options,
    )


def generate_primary_ids(
    source: str,
    *,
    version: str | None = None,
    show_progress: bool = True,
    inputs: dict[str, Path | str] | None = None,
    **options: Any,
) -> BaseMappingSet:
    """Return a mapping set carrying only *source*'s full current-ID list.

    The mappings list is empty; the set exists to drive ``to_pri_ids()``. Use
    this to get the authoritative ID list without parsing the withdrawn file.

    Args:
        source: Datasource name; see ``sources("primary_ids")``.
        version: Release to build. The latest is used when ``None``.
        show_progress: Whether to show progress bars.
        inputs: Local input files keyed as in the config's ``download_urls``.
        **options: Source-specific options; see :func:`generate_ids`.

    Returns:
        A mapping set with no mappings and ``_primary_ids`` populated.
    """
    return _generate(
        source,
        "primary_ids",
        version=version,
        show_progress=show_progress,
        inputs=inputs,
        **options,
    )


def generate_primary_labels(
    source: str,
    *,
    version: str | None = None,
    show_progress: bool = True,
    inputs: dict[str, Path | str] | None = None,
    **options: Any,
) -> BaseMappingSet:
    """Return a mapping set carrying only *source*'s full current-label list.

    Args:
        source: Datasource name; see ``sources("primary_labels")``.
        version: Release to build. The latest is used when ``None``.
        show_progress: Whether to show progress bars.
        inputs: Local input files keyed as in the config's ``download_urls``.
        **options: Source-specific options; see :func:`generate_ids`.

    Returns:
        A mapping set with no mappings and ``_primary_labels`` populated.
    """
    return _generate(
        source,
        "primary_labels",
        version=version,
        show_progress=show_progress,
        inputs=inputs,
        **options,
    )


def combine_mapping_sets(
    id_mappings: BaseMappingSet | None,
    synonym_mappings: BaseMappingSet | None,
) -> BaseMappingSet:
    """Combine two mapping sets into one.

    Args:
        id_mappings: First mapping set (e.g. ID mappings).
        synonym_mappings: Second mapping set (e.g. synonym mappings).

    Returns:
        Combined mapping set.

    Raises:
        ValueError: If both mapping sets are ``None``.
    """
    if id_mappings is None and synonym_mappings is None:
        msg = "At least one mapping set must be provided"
        raise ValueError(msg)
    if id_mappings is None:
        return synonym_mappings  # type: ignore[return-value]
    if synonym_mappings is None:
        return id_mappings
    combined_mappings = list(id_mappings.mappings or [])
    combined_mappings.extend(synonym_mappings.mappings or [])
    result = copy.copy(id_mappings)
    result.mappings = combined_mappings
    return result


# Output helpers

_FORMAT_EXTENSIONS: dict[str, str] = {
    "rdf": ".ttl",
    "owl": "_owl.ttl",
    "json": ".json",
}


def _output_filename(base_name: str, fmt: str) -> Path:
    """Return the default filename for *base_name* in format *fmt*."""
    ext = _FORMAT_EXTENSIONS.get(fmt, ".tsv")
    if fmt == "owl":
        return Path(f"{base_name}{ext}")
    return Path(f"{base_name}_{fmt}{ext}")


def save(
    mapping_set: BaseMappingSet,
    output_format: str,
    output: Path | str | None = None,
    *,
    base_name: str,
) -> Path:
    """Write *mapping_set* and return the path that was written.

    Delegates to :meth:`~pysec2pri.parsers.base.BaseMappingSet.save` for
    single formats and :func:`write_all_formats` for ``"all"``.

    Args:
        mapping_set: The mapping set to write.
        output_format: One of ``sssom``, ``sec2pri``, ``pri_ids``,
            ``name2synonym``, ``label_sec2pri``, ``pri_labels``,
            ``rdf``, ``json``, ``owl``, or ``all``.
        output: Explicit output path or directory.  When ``None``, a
            default name derived from *base_name* is used.
        base_name: Stem used to derive file names, e.g. ``"hgnc_2026-04-07"``.

    Returns:
        The directory (for ``"all"``) or file path that was written.
    """
    out = Path(output) if output else None

    if output_format == "all":
        if out is None:
            out_dir = Path(base_name)
        elif out.suffix:
            out_dir = out.parent / base_name
        else:
            out_dir = out
        write_all_formats(mapping_set, out_dir, base_name)
        return out_dir

    # Resolve output path, then delegate to the mapping-set method
    if out is None:
        out_path = _output_filename(base_name, output_format)
    elif out.is_dir():
        out_path = out / _output_filename(base_name, output_format).name
    else:
        out_path = out

    return mapping_set.save(output_format, out_path)


def write_all_formats(
    mapping_set: BaseMappingSet,
    output_dir: Path,
    base_name: str,
    include_name2synonym: bool = True,
) -> None:
    """Write mapping set in all output formats to a directory.

    Args:
        mapping_set: The mapping set to write.
        output_dir: Directory to write files to.
        base_name: Base name for output files (e.g., "chebi_3star_245").
        include_name2synonym: Whether to include name2synonym format.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    write_sssom(mapping_set, output_dir / f"{base_name}_sssom.tsv")

    if isinstance(mapping_set, IdMappingSet):
        write_sec2pri(mapping_set, output_dir / f"{base_name}_sec2pri.tsv")
        write_pri_ids(mapping_set, output_dir / f"{base_name}_pri_ids.txt")

    if isinstance(mapping_set, LabelMappingSet):
        write_label_sec2pri(mapping_set, output_dir / f"{base_name}_label_sec2pri.tsv")
        if include_name2synonym:
            write_name2synonym(mapping_set, output_dir / f"{base_name}_name2synonym.tsv")

    write_rdf(mapping_set, output_dir / f"{base_name}.ttl")
    write_json(mapping_set, output_dir / f"{base_name}.json")
    write_owl(mapping_set, output_dir / f"{base_name}_owl.ttl")


def write_diff_output(
    result: MappingDiff,
    output_path: Path,
) -> None:
    """Write diff results to a TSV file.

    Args:
        result: MappingDiff object with added/removed/changed mappings.
        output_path: Path to write the TSV file.
    """
    import polars as pl

    dfs = []

    if result.added_count > 0:
        added_df = result.added.with_columns(
            pl.lit("added").alias("change_type"),
            pl.lit(None).alias("old_subject_id"),
        ).select(
            [
                "change_type",
                "object_id",
                pl.col("subject_id").alias("new_subject_id"),
                "old_subject_id",
            ]
        )
        dfs.append(added_df)

    if result.removed_count > 0:
        removed_df = result.removed.with_columns(
            pl.lit("removed").alias("change_type"),
            pl.lit(None).alias("new_subject_id"),
        ).select(
            [
                "change_type",
                "object_id",
                "new_subject_id",
                pl.col("subject_id").alias("old_subject_id"),
            ]
        )
        dfs.append(removed_df)

    if result.changed_count > 0:
        changed_df = result.changed.with_columns(
            pl.lit("changed").alias("change_type"),
        ).select(
            [
                "change_type",
                "object_id",
                "new_subject_id",
                "old_subject_id",
            ]
        )
        dfs.append(changed_df)

    if dfs:
        combined = pl.concat(dfs)
        combined.write_csv(output_path, separator="\t")


def load_mapping(path: Path | str) -> IdMappingSet:
    """Load an ID mapping set from an SSSOM TSV file.

    Produces the same :class:`~pysec2pri.parsers.base.IdMappingSet` a fresh
    parse would, ready to pass to :func:`resolve_ids`.

    Args:
        path: Path to the SSSOM TSV file to load.
    """
    from mapkgsutils.exports import read_sssom

    return read_sssom(path, mapping_set_class=IdMappingSet, on="id")


def load_label_mapping(path: Path | str) -> LabelMappingSet:
    """Load a label mapping set from an SSSOM TSV file.

    Produces the same :class:`~pysec2pri.parsers.base.LabelMappingSet` a fresh
    parse would, ready to pass to :func:`resolve_labels`.

    Args:
        path: Path to the SSSOM TSV file to load.
    """
    from mapkgsutils.exports import read_sssom

    return read_sssom(path, mapping_set_class=LabelMappingSet, on="label")


def resolve_ids(
    input_path: Path | str | list[str],
    mapping_set: BaseMappingSet,
    at: str | list[str] | None = None,
    *,
    output_path: Path | str | None = None,
    suffix: str = "_primary",
    sep: str | None = None,
    synonyms: str | None = None,
    label_mapping_set: BaseMappingSet | None = None,
    xref: str | None = None,
    xref_mapping: XrefMapping | None = None,
    report_path: Path | str | None = None,
) -> pd.DataFrame | str | list[str]:
    r"""Resolve secondary IDs to primary IDs.

    Direct lookup: when *input_path* is a plain identifier string or a list
    of identifier strings (i.e. not a path to an existing file), the function
    returns the resolved primary ID(s).  *at*, *output_path*, *suffix*, and
    *sep* are ignored in this mode::

        resolve_ids("HMDB00001", hmdb_ms)  # -> "HMDB:HMDB0000001"
        resolve_ids(["HMDB00001", "HMDB00002"], hmdb_ms)  # -> ["...", "..."]

    DataFrame mode: when *input_path* points to an existing TSV/CSV
    file, *at* is required.  The file is read with
    ``pandas.read_csv`` and for each column named in *at* a new column
    ``<col><suffix>`` is appended containing the resolved primary IDs.
    Identifiers not present in *mapping_set* are kept unchanged.

    Args:
        input_path: An identifier string, a list of identifier strings, or
            the path to a TSV/CSV file.
        mapping_set: A :class:`~pysec2pri.parsers.base.BaseMappingSet`
            (e.g. the result of ``generate_ids("hgnc")``).
        at: Column name(s) to resolve.  Required in DataFrame mode;
            ignored in direct-lookup mode.
        output_path: If given, the resulting DataFrame is written to this
            path (DataFrame mode only).
        suffix: Suffix appended to each resolved column name
            (default ``"_primary"``).
        sep: Delimiter for reading the file.  Inferred from the extension
            when ``None`` (``"\\t"`` for ``.tsv``, ``","`` otherwise).
        xref: *DataFrame mode only.* Column with a per-row cross-reference
            token, passed through to :func:`~pysec2pri.update.update_ids`.
        xref_mapping: The :class:`~mapkgsutils.context.XrefMapping` crosswalk
            table to resolve *xref* tokens against. Required when *xref* is
            given.
        report_path: When given, every disambiguation attempt (from
            *synonyms* and/or *xref*) is logged to this TSV.

    Returns:
        A resolved identifier string, a list of resolved strings (direct-lookup
        mode), or a :class:`pandas.DataFrame` with one additional column per
        entry in *at* (DataFrame mode).
    """
    import pandas as pd

    from pysec2pri.update import build_lookup, update_ids

    # list direct-lookup mode
    if isinstance(input_path, list):
        lkp = build_lookup(mapping_set)
        return [lkp.get(v, v) for v in input_path]

    # single-value lookup mode
    input_path_obj = Path(input_path)
    if not input_path_obj.exists():
        lkp = build_lookup(mapping_set)
        return lkp.get(str(input_path), str(input_path))

    # DataFrame mode
    if at is None:
        raise TypeError("resolve_ids() requires 'at' when input_path is a file")

    if sep is None:
        sep = "\t" if input_path_obj.suffix.lower() == ".tsv" else ","

    df = pd.read_csv(input_path_obj, sep=sep, dtype=str)
    lkp = build_lookup(mapping_set)
    result = update_ids(
        df,
        mapping_set,
        at=at,
        suffix=suffix,
        lookup=lkp,
        synonyms=synonyms,
        label_mapping_set=label_mapping_set,
        xref=xref,
        xref_mapping=xref_mapping,
        report_path=report_path,
    )

    if output_path is not None:
        output_path = Path(output_path)
        out_sep = "\t" if output_path.suffix.lower() == ".tsv" else ","
        result.to_csv(output_path, sep=out_sep, index=False)

    return result


def resolve_labels(
    input_path: Path | str | list[str],
    mapping_set: BaseMappingSet,
    at: str | list[str] | None = None,
    *,
    output_path: Path | str | None = None,
    suffix: str = "_current",
    sep: str | None = None,
    synonyms: str | None = None,
    xref: str | None = None,
    xref_mapping: XrefMapping | None = None,
    report_path: Path | str | None = None,
) -> pd.DataFrame | str | list[str]:
    r"""Resolve previous/alias labels to current labels.

    Direct lookup: when *input_path* is a plain label string or a list
    of label strings (i.e. not a path to an existing file), the function
    returns the resolved current label(s).  *at*, *output_path*, *suffix*,
    and *sep* are ignored in this mode::

        resolve_labels("Ibuprofen", chebi_ms)  # -> "ibuprofen"
        resolve_labels(["Ibuprofen", "Glucose"], chebi_ms)  # -> ["...", "..."]

    DataFrame mode: when *input_path* points to an existing TSV/CSV
    file, *at* is required.  For each column named in *at* a new
    column ``<col><suffix>`` is appended containing the resolved current
    labels.  Symbols not present in *mapping_set* are kept unchanged.

    Args:
        input_path: A label string, a list of label strings, or the path
            to a TSV/CSV file.
        mapping_set: A :class:`~pysec2pri.parsers.base.LabelMappingSet`
            (e.g. the result of ``generate_labels("hgnc")``).
        at: Column name(s) to resolve.  Required in DataFrame mode;
            ignored in direct-lookup mode.
        output_path: If given, the resulting DataFrame is written to this
            path (DataFrame mode only).
        suffix: Suffix appended to each resolved column name
            (default ``"_current"``).
        sep: Delimiter for reading the file.  Inferred from the extension
            when ``None`` (``"\\t"`` for ``.tsv``, ``","`` otherwise).
        xref: *DataFrame mode only.* Column with a per-row cross-reference
            token, passed through to
            :func:`~pysec2pri.update.update_labels`.
        xref_mapping: The :class:`~mapkgsutils.context.XrefMapping` crosswalk
            table to resolve *xref* tokens against. Required when *xref* is
            given.
        report_path: When given, every disambiguation attempt (from
            *synonyms* and/or *xref*) is logged to this TSV.

    Returns:
        A resolved label string, a list of resolved strings (direct-lookup
        mode), or a :class:`pandas.DataFrame` with one additional column per
        entry in *at* (DataFrame mode).
    """
    import pandas as pd

    from pysec2pri.update import build_label_lookup, update_labels

    # list direct-lookup mode
    if isinstance(input_path, list):
        lkp = build_label_lookup(mapping_set)
        return [lkp.get(v, v) for v in input_path]

    # single-value direct-lookup mode
    input_path_obj = Path(input_path)
    if not input_path_obj.exists():
        lkp = build_label_lookup(mapping_set)
        return lkp.get(str(input_path), str(input_path))

    # DataFrame mode
    if at is None:
        raise TypeError("resolve_labels() requires 'at' when input_path is a file")

    if sep is None:
        sep = "\t" if input_path_obj.suffix.lower() == ".tsv" else ","

    df = pd.read_csv(input_path_obj, sep=sep, dtype=str)
    lkp = build_label_lookup(mapping_set)
    result = pd.DataFrame(
        update_labels(
            df,
            mapping_set,
            at=at,
            suffix=suffix,
            lookup=lkp,
            synonyms=synonyms,
            xref=xref,
            xref_mapping=xref_mapping,
            report_path=report_path,
        )
    )

    if output_path is not None:
        output_path = Path(output_path)
        out_sep = "\t" if output_path.suffix.lower() == ".tsv" else ","
        result.to_csv(output_path, sep=out_sep, index=False)

    return result


def list_versions(datasource: str) -> Any:
    """List all available archive versions for a datasource.

    For datasources that publish versioned archives (ChEBI, HGNC, UniProt),
    this queries the remote archive index and returns all available version
    strings sorted in ascending order.

    NCBI and HMDB do not maintain versioned archives; calling this function
    for those datasources raises :class:`ValueError`.

    Args:
        datasource: Datasource name, one of ``"chebi"``, ``"hgnc"``, or
            ``"uniprot"``.

    Returns:
        Sorted list of version strings.  Format depends on the datasource:

        - **chebi**: integer release numbers, e.g. ``["200", ..., "245"]``
        - **hgnc**: ISO dates, e.g. ``["2023-01-01", ..., "2026-04-07"]``
        - **uniprot**: release IDs, e.g. ``["2024_01", "2024_02", ...]``

    Raises:
        ValueError: If *datasource* is unknown or has no versioned archive.
    """
    from pysec2pri.download import list_versions as _list_versions

    return _list_versions(datasource)


def find_ambiguous(
    mapping_set: BaseMappingSet,
) -> AmbiguousMappingSet:
    """Find identifiers that are ambiguous in *mapping_set*.

    An identifier is ambiguous when it appears both as a ``subject_id`` (i.e. a
    secondary/previous term) and as a current primary identifier.  Such
    entries cannot be automatically resolved without risk of corrupting
    references that are already current.

    This is a convenience wrapper around
    :meth:`~pysec2pri.parsers.base.BaseMappingSet.find_ambiguous`.

    Args:
        mapping_set: A :class:`~pysec2pri.parsers.base.BaseMappingSet`
            (e.g. the result of ``generate_ids("hgnc")``).

    Returns:
        An :class:`~pysec2pri.parsers.base.AmbiguousMappingSet` whose
        ``mappings`` list contains one entry for each conflicting subject, with a
        ``comment`` explaining the conflict.
    """
    return mapping_set.find_ambiguous()
