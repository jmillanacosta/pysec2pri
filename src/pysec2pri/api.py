"""Main API functions for pysec2pri.

This module provides high-level functions for parsing biological database
secondary-to-primary mapping files and generating output.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from pysec2pri.exports import (
    write_json,
    write_name2synonym,
    write_output,
    write_owl,
    write_pri_ids,
    write_rdf,
    write_sec2pri,
    write_sssom,
    write_symbol2prev,
)

if TYPE_CHECKING:
    import pandas as pd

    from pysec2pri.diff import MappingDiff
    from pysec2pri.parsers.base import Sec2PriMappingSet

__all__ = [
    "combine_mapping_sets",
    "generate_chebi",
    "generate_chebi_synonyms",
    "generate_hgnc",
    "generate_hgnc_symbols",
    "generate_hmdb",
    "generate_hmdb_proteins",
    "generate_ncbi",
    "generate_ncbi_symbols",
    "generate_uniprot",
    "generate_wikidata",
    "resolve_ids",
    "save",
    "write_all_formats",
    "write_diff_output",
    "write_json",
    "write_name2synonym",
    "write_output",
    "write_owl",
    "write_rdf",
    "write_sec2pri",
    "write_sssom",
    "write_symbol2prev",
]


def _auto_download(datasource: str, version: str | None = None) -> dict[str, Path]:
    """Download all files for *datasource* into a temp dir."""
    import tempfile

    from pysec2pri.download import download_datasource

    tmpdir = Path(tempfile.mkdtemp(prefix=f"pysec2pri_{datasource}_"))
    return download_datasource(datasource, tmpdir, version=version)


def generate_chebi(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
    subset: str = "3star",
    mapping_sets: str = "ids",
) -> Sec2PriMappingSet:
    """Return ChEBI mappings (IDs, synonyms, or both).

    Downloads the latest release automatically when ``input_path`` is omitted.
    Pass an SDF file (releases < 245) or a directory of TSV flat files
    (releases >= 245) to use a local copy.

    Args:
        input_path: Local SDF file or TSV directory. Auto-downloaded if ``None``.
        version: Release number (e.g. ``"245"``).
        show_progress: Whether to show progress bars.
        subset: ``"3star"`` (default) or ``"complete"``.
        mapping_sets: ``"ids"`` (default), ``"synonyms"``, or ``"all"``.
    """
    import tempfile

    from pysec2pri.download import check_chebi_release
    from pysec2pri.parsers import ChEBIParser
    from pysec2pri.parsers.chebi import ChEBIDownloader

    if input_path is None:
        if version is None:
            version = check_chebi_release().version or "245"
        downloader = ChEBIDownloader(version=version, subset=subset)
        tmpdir = Path(tempfile.mkdtemp(prefix=f"pysec2pri_chebi_{version}_"))
        downloader.download(tmpdir, version=version)
        input_path = tmpdir

    parser = ChEBIParser(version=version, show_progress=show_progress, subset=subset)

    if mapping_sets == "synonyms":
        return parser.parse_synonyms(Path(input_path))
    if mapping_sets == "all":
        ids = parser.parse(Path(input_path))
        syns = parser.parse_synonyms(Path(input_path))
        return combine_mapping_sets(ids, syns)
    return parser.parse(Path(input_path))


def generate_chebi_synonyms(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
    subset: str = "3star",
) -> Sec2PriMappingSet:
    """Return ChEBI synonym (name) mappings."""
    return generate_chebi(
        input_path=input_path,
        version=version,
        show_progress=show_progress,
        subset=subset,
        mapping_sets="synonyms",
    )


def generate_hgnc(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return HGNC secondary to primary ID mappings.

    Downloads the withdrawn file automatically when ``input_path`` is omitted.

    Args:
        input_path: Local HGNC withdrawn TSV. Auto-downloaded if ``None``.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import HGNCParser

    if input_path is None:
        input_path = _auto_download("hgnc", version)["withdrawn"]

    parser = HGNCParser(version=version, show_progress=show_progress)
    return parser.parse(Path(input_path))


def generate_hgnc_symbols(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
    statuses: list[str] | None = None,
) -> Sec2PriMappingSet:
    """Return HGNC symbol to previous-symbol mappings.

    Downloads the complete set file automatically when ``input_path`` is omitted.

    Args:
        input_path: Local HGNC complete set TSV. Auto-downloaded if ``None``.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
        statuses: Entry statuses to include (e.g. ``["Approved"]``).
    """
    from pysec2pri.parsers import HGNCParser

    if input_path is None:
        input_path = _auto_download("hgnc", version)["complete"]

    parser = HGNCParser(version=version, show_progress=show_progress)
    return parser.parse_symbols(Path(input_path), statuses=statuses)


def generate_ncbi(
    input_path: Path | str | None = None,
    tax_id: str = "9606",
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return NCBI Gene secondary to primary ID mappings.

    Downloads the gene history file automatically when ``input_path`` is omitted.

    Args:
        input_path: Local gene_history file. Auto-downloaded if ``None``.
        tax_id: NCBI taxonomy ID to filter (default: ``"9606"`` for human).
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import NCBIParser

    if input_path is None:
        input_path = _auto_download("ncbi", version)["gene_history"]

    parser = NCBIParser(version=version, show_progress=show_progress)
    return parser.parse(Path(input_path), tax_id=tax_id)


def generate_ncbi_symbols(
    input_path: Path | str | None = None,
    tax_id: str = "9606",
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return NCBI Gene symbol to previous-symbol mappings.

    Downloads the gene_info file automatically when ``input_path`` is omitted.

    Args:
        input_path: Local gene_info file. Auto-downloaded if ``None``.
        tax_id: NCBI taxonomy ID to filter (default: ``"9606"`` for human).
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import NCBIParser

    if input_path is None:
        input_path = _auto_download("ncbi", version)["gene_info"]

    parser = NCBIParser(version=version, show_progress=show_progress)
    return parser.parse_symbols(Path(input_path), tax_id=tax_id)


def generate_uniprot(
    input_path: Path | str | None = None,
    delac_file: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return UniProt secondary to primary accession mappings.

    Downloads sec_ac.txt and delac_sp.txt automatically when ``input_path``
    is omitted.

    Args:
        input_path: Local sec_ac.txt. Auto-downloaded if ``None``.
        delac_file: Local delac_sp.txt.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import UniProtParser

    if input_path is None:
        files = _auto_download("uniprot", version)
        input_path = files.get("sec_ac") or next(iter(files.values()))
        if delac_file is None:
            delac_file = files.get("delac_sp")

    parser = UniProtParser(version=version, show_progress=show_progress)
    return parser.parse(
        Path(input_path),
        delac_path=Path(delac_file) if delac_file else None,
    )


def generate_hmdb(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return HMDB metabolite secondary to primary accession mappings.

    Downloads hmdb_metabolites.xml automatically when ``input_path`` is omitted.

    Args:
        input_path: Local hmdb_metabolites.xml (or .zip/.gz). Auto-downloaded if ``None``.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import HMDBParser

    if input_path is None:
        input_path = _auto_download("hmdb", version)["metabolites"]

    parser = HMDBParser(version=version, show_progress=show_progress)
    return parser.parse(Path(input_path))


def generate_hmdb_proteins(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return HMDB protein secondary to primary accession mappings.

    Downloads hmdb_proteins.xml automatically when ``input_path`` is omitted.

    Args:
        input_path: Local hmdb_proteins.xml (or .zip/.gz). Auto-downloaded if ``None``.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import HMDBParser

    if input_path is None:
        input_path = _auto_download("hmdb", version)["proteins"]

    parser = HMDBParser(version=version, show_progress=show_progress)
    return parser.parse_proteins(Path(input_path))


def generate_wikidata(
    input_path: Path | str | None = None,
    entity_type: str | None = None,
    version: str | None = None,
    endpoint: str | None = None,
    show_progress: bool = True,
    test_subset: bool = False,
) -> Sec2PriMappingSet:
    """Return Wikidata redirect mappings via SPARQL (or a pre-downloaded TSV).

    Queries the QLever Wikidata endpoint when ``input_path`` is omitted.
    If ``entity_type`` is ``None``, all entity types (metabolites, genes,
    proteins) are queried and combined.

    Args:
        input_path: Pre-downloaded TSV file. Queries SPARQL if ``None``.
        entity_type: ``"metabolites"``, ``"chemicals"``, ``"genes"``, or
        `"proteins"``. Queries all types when ``None``.
        version: Version string for metadata (defaults to today's date).
        endpoint: Custom SPARQL endpoint URL.
        show_progress: Whether to show progress bars.
        test_subset: Use test queries limited to 10 results.
    """
    from pysec2pri.parsers import WikidataParser

    parser = WikidataParser(
        version=version,
        show_progress=show_progress,
        entity_type=entity_type or "metabolites",
        endpoint=endpoint,
        test_subset=test_subset,
    )

    if input_path is not None:
        return parser.parse_from_file(Path(input_path))

    if entity_type is None:
        return parser.parse_all()

    return parser.parse()


def combine_mapping_sets(
    id_mappings: Sec2PriMappingSet | None,
    synonym_mappings: Sec2PriMappingSet | None,
) -> Sec2PriMappingSet:
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
    combined = list(id_mappings.mappings or [])
    combined.extend(synonym_mappings.mappings or [])
    id_mappings.mappings = combined
    return id_mappings


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
    mapping_set: Sec2PriMappingSet,
    output_format: str,
    output: Path | str | None = None,
    *,
    base_name: str,
) -> Path:
    """Write *mapping_set* and return the path that was written.

    Delegates to :meth:`~pysec2pri.parsers.base.Sec2PriMappingSet.save` for
    single formats and :func:`write_all_formats` for ``"all"``.

    Args:
        mapping_set: The mapping set to write.
        output_format: One of ``sssom``, ``sec2pri``, ``pri_ids``,
            ``name2synonym``, ``symbol2prev``, ``rdf``, ``json``, ``owl``,
            or ``all``.
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
    mapping_set: Sec2PriMappingSet,
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
    write_sec2pri(mapping_set, output_dir / f"{base_name}_sec2pri.tsv")
    write_pri_ids(mapping_set, output_dir / f"{base_name}_pri_ids.tsv")
    write_rdf(mapping_set, output_dir / f"{base_name}.ttl")
    write_json(mapping_set, output_dir / f"{base_name}.json")
    write_owl(mapping_set, output_dir / f"{base_name}_owl.ttl")

    if include_name2synonym:
        write_name2synonym(mapping_set, output_dir / f"{base_name}_name2synonym.tsv")


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


def resolve_ids(
    input_path: Path | str,
    mapping_set: Sec2PriMappingSet,
    at: str | list[str],
    *,
    output_path: Path | str | None = None,
    suffix: str = "_primary",
    sep: str | None = None,
) -> pd.DataFrame:
    r"""Read a delimited file, resolve secondary IDs to primary IDs, and return a DataFrame.

    The file is read with ``pandas.read_csv`` (auto-detecting TSV vs CSV by
    extension unless *sep* is supplied).  For each column named in *at*, a new
    column ``<col><suffix>`` is appended containing the resolved primary IDs.
    Identifiers not present in *mapping_set* are kept unchanged.

    Args:
        input_path: Path to the input TSV or CSV file.
        mapping_set: A :class:`~pysec2pri.parsers.base.Sec2PriMappingSet`
            (e.g. the result of ``generate_hgnc()``).
        at: Column name or list of column names that contain identifiers to
            resolve.
        output_path: If given, the resulting DataFrame is written back to this
            path in the same delimiter format (TSV when the extension is
            ``.tsv``, CSV otherwise).
        suffix: Suffix appended to each resolved column name
            (default ``"_primary"``).
        sep: Delimiter to use when reading the file.  Inferred from the file
            extension when ``None`` (``"\\t"`` for ``.tsv``, ``","`` for all
            others).

    Returns:
        :class:`pandas.DataFrame` with one additional column per entry in *at*.
    """
    import pandas as pd

    from pysec2pri.update_ids import build_lookup, update_ids

    input_path = Path(input_path)

    if sep is None:
        sep = "\t" if input_path.suffix.lower() == ".tsv" else ","

    df = pd.read_csv(input_path, sep=sep, dtype=str)

    lkp = build_lookup(mapping_set)
    result: pd.DataFrame = update_ids(df, mapping_set, at=at, suffix=suffix, lookup=lkp)

    if output_path is not None:
        output_path = Path(output_path)
        out_sep = "\t" if output_path.suffix.lower() == ".tsv" else ","
        result.to_csv(output_path, sep=out_sep, index=False)

    return result
