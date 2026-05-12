"""Main API functions for pysec2pri.

This module provides high-level functions for parsing biological database
secondary-to-primary mapping files and generating output.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from pysec2pri.exports import (
    write_name2synonym,
    write_pri_ids,
    write_sec2pri,
    write_sssom,
    write_symbol2prev,
)

if TYPE_CHECKING:
    from pysec2pri.diff import MappingDiff
    from pysec2pri.parsers.base import Sec2PriMappingSet

__all__ = [
    "combine_mapping_sets",
    "parse_chebi",
    "parse_chebi_synonyms",
    "parse_hgnc",
    "parse_hgnc_symbols",
    "parse_hmdb",
    "parse_hmdb_proteins",
    "parse_ncbi",
    "parse_ncbi_symbols",
    "parse_uniprot",
    "parse_wikidata",
    "write_all_formats",
    "write_diff_output",
    "write_name2synonym",
    "write_sec2pri",
    "write_sssom",
    "write_symbol2prev",
]


def _download_chebi(
    version: str | None,
    subset: str,
) -> tuple[Path, str]:
    """Download ChEBI files to a temp directory and return (dir, version)."""
    import tempfile

    from pysec2pri.download import check_chebi_release
    from pysec2pri.parsers.chebi import ChEBIDownloader

    if version is None:
        release_info = check_chebi_release()
        version = release_info.version or "245"
    downloader = ChEBIDownloader(version=version, subset=subset)
    tmpdir = Path(tempfile.mkdtemp(prefix=f"pysec2pri_chebi_{version}_"))
    downloader.download(tmpdir, version=version)
    return tmpdir, version


def parse_chebi(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
    subset: str = "3star",
) -> Sec2PriMappingSet:
    """Parse ChEBI data and return ID mappings only.

    Downloads the latest release automatically when ``input_path`` is omitted.
    The ``input_path`` can be an SDF file (releases < 245) or a directory of TSV
    flat files (releases >= 245, default).

    Args:
        input_path: Path to an SDF file or a directory of TSV flat files.
            If ``None``, the latest release is downloaded automatically.
        version: Release number (e.g. ``"245"``).
        show_progress: Whether to show progress bars during parsing.
        subset: ``"3star"`` (default) or ``"complete"``.

    Returns:
        Sec2PriMappingSet with ID mappings only.
    """
    from pysec2pri.parsers import ChEBIParser

    if input_path is None:
        input_path, version = _download_chebi(version, subset)

    parser = ChEBIParser(version=version, show_progress=show_progress, subset=subset)
    return parser.parse(Path(input_path))


def parse_chebi_synonyms(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
    subset: str = "3star",
) -> Sec2PriMappingSet:
    """Parse ChEBI data and return synonym (name) mappings only.

    Downloads the latest release automatically when ``input_path`` is omitted.

    Args:
        input_path: Path to an SDF file or a directory of TSV flat files.
            If ``None``, the latest release is downloaded automatically.
        version: Release number (e.g. ``"245"``).
        show_progress: Whether to show progress bars during parsing.
        subset: ``"3star"`` (default) or ``"complete"``.

    Returns:
        Sec2PriMappingSet with synonym mappings.
    """
    from pysec2pri.parsers import ChEBIParser

    if input_path is None:
        input_path, version = _download_chebi(version, subset)

    parser = ChEBIParser(version=version, show_progress=show_progress, subset=subset)
    return parser.parse_synonyms(Path(input_path))


def parse_hmdb(
    input_file: Path | str,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Parse HMDB metabolites XML and extract secondary-to-primary ID mappings.

    Args:
        input_file: Path to ``hmdb_metabolites.xml`` (or ``.zip``/``.gz``).
        version: Version string for metadata.
        show_progress: Whether to show progress bars.

    Returns:
        Sec2PriMappingSet for metabolite accessions.
    """
    from pysec2pri.parsers import HMDBParser

    parser = HMDBParser(version=version, show_progress=show_progress)
    return parser.parse(Path(input_file))


def parse_hmdb_proteins(
    input_file: Path | str,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Parse HMDB proteins XML and extract secondary-to-primary ID mappings.

    Args:
        input_file: Path to ``hmdb_proteins.xml`` (or ``.zip``/``.gz``).
        version: Version string for metadata.
        show_progress: Whether to show progress bars.

    Returns:
        Sec2PriMappingSet for protein accessions (``HMDBP`` prefix).
    """
    from pysec2pri.parsers import HMDBParser

    parser = HMDBParser(version=version, show_progress=show_progress)
    return parser.parse_proteins(Path(input_file))


def parse_hgnc(
    withdrawn_file: Path | str,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Parse HGNC withdrawn file and extract ID mappings."""
    from pysec2pri.parsers import HGNCParser

    parser = HGNCParser(version=version, show_progress=show_progress)
    return parser.parse(Path(withdrawn_file))


def parse_hgnc_symbols(
    complete_set_file: Path | str,
    version: str | None = None,
    show_progress: bool = True,
    statuses: list[str] | None = None,
) -> Sec2PriMappingSet:
    """Parse HGNC complete set file and extract symbol mappings.

    Args:
        complete_set_file: Path to the HGNC complete set TSV file.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
        statuses: Entry statuses to include (e.g. ``["Approved"]``).
        If ``None`` (default), all entries are included.
    """
    from pysec2pri.parsers import HGNCParser

    parser = HGNCParser(version=version, show_progress=show_progress)
    return parser.parse_symbols(Path(complete_set_file), statuses=statuses)


def parse_ncbi(
    history_file: Path | str,
    tax_id: str = "9606",
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Parse NCBI Gene history file and extract ID mappings."""
    from pysec2pri.parsers import NCBIParser

    parser = NCBIParser(version=version, show_progress=show_progress)
    return parser.parse(Path(history_file), tax_id=tax_id)


def parse_ncbi_symbols(
    gene_info_file: Path | str,
    tax_id: str = "9606",
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Parse NCBI Gene info file and extract symbol mappings."""
    from pysec2pri.parsers import NCBIParser

    parser = NCBIParser(version=version, show_progress=show_progress)
    return parser.parse_symbols(Path(gene_info_file), tax_id=tax_id)


def parse_uniprot(
    sec_ac_file: Path | str,
    delac_file: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Parse UniProt secondary accession files and extract ID mappings."""
    from pysec2pri.parsers import UniProtParser

    parser = UniProtParser(version=version, show_progress=show_progress)
    return parser.parse(
        Path(sec_ac_file),
        delac_path=Path(delac_file) if delac_file else None,
    )


def parse_wikidata(
    input_file: Path | str | None = None,
    entity_type: str | None = None,
    version: str | None = None,
    endpoint: str | None = None,
    show_progress: bool = True,
    test_subset: bool = False,
) -> Sec2PriMappingSet:
    """Parse Wikidata redirects via SPARQL or from a pre-downloaded file.

    Queries the QLever Wikidata endpoint (faster than official endpoint)
    for redirect mappings.

    If entity_type is None (default), queries ALL entity types defined
    in the config (metabolites, genes, proteins) and combines results.

    Args:
        input_file: Optional path to pre-downloaded TSV file. If None,
                   queries SPARQL endpoint directly.
        entity_type: Type of entities to query. One of:
                    'metabolites', 'chemicals', 'genes', 'proteins'.
                    If None, queries all entity types.
        version: Version string for the mappings (defaults to today's date).
        endpoint: Optional custom SPARQL endpoint URL.
        show_progress: Whether to show progress bars.
        test_subset: Whether to use test queries (LIMIT 10 results).

    Returns:
        Sec2PriMappingSet with computed cardinalities.

    Example:
        >>> from pysec2pri.api import parse_wikidata
        >>> # Query all entity types
        >>> mappings = parse_wikidata()
        >>> # Query only metabolites
        >>> mappings = parse_wikidata(entity_type="metabolites")
    """
    from pysec2pri.parsers import WikidataParser

    parser = WikidataParser(
        version=version,
        show_progress=show_progress,
        entity_type=entity_type or "metabolites",
        endpoint=endpoint,
        test_subset=test_subset,
    )

    if input_file is not None:
        return parser.parse_from_file(Path(input_file))

    # If no entity_type specified, query all types
    if entity_type is None:
        return parser.parse_all()

    return parser.parse()


def combine_mapping_sets(
    id_mappings: Sec2PriMappingSet | None,
    synonym_mappings: Sec2PriMappingSet | None,
) -> Sec2PriMappingSet:
    """Combine ID and synonym mapping sets into one.

    Args:
        id_mappings: Mapping set with ID mappings.
        synonym_mappings: Mapping set with synonym mappings.

    Returns:
        Combined mapping set with both ID and synonym mappings.

    Raises:
        ValueError: If both mapping sets are None.
    """
    if id_mappings is None and synonym_mappings is None:
        msg = "At least one mapping set must be provided"
        raise ValueError(msg)

    if id_mappings is None:
        return synonym_mappings  # type: ignore[return-value]
    if synonym_mappings is None:
        return id_mappings

    # Merge mappings from both sets
    combined_mappings = list(id_mappings.mappings or [])
    combined_mappings.extend(synonym_mappings.mappings or [])

    # Use id_mappings as base and update with combined
    result = id_mappings
    result.mappings = combined_mappings
    return result


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
