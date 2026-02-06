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


def parse_chebi(
    input_file: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
    subset: str = "3star",
    *,
    secondary_ids_path: Path | str | None = None,
    compounds_path: Path | str | None = None,
) -> Sec2PriMappingSet:
    """Parse ChEBI data files and extract ID mappings.

    For releases >= 245: use secondary_ids_path and compounds_path (TSV).
    For releases < 245: use input_file (SDF format).

    Args:
        input_file: Path to the ChEBI SDF file (legacy format < 245).
        version: Version/release identifier for the datasource.
        show_progress: Whether to show progress bars during parsing.
        subset: "3star" or "complete" - which compounds to include.
        secondary_ids_path: Path to secondary_ids.tsv (new format >= 245).
        compounds_path: Path to compounds.tsv for 3-star filtering.

    Returns:
        Sec2PriMappingSet with computed cardinalities.
    """
    from pysec2pri.parsers import ChEBIParser

    parser = ChEBIParser(version=version, show_progress=show_progress, subset=subset)
    return parser.parse(
        input_path=Path(input_file) if input_file else None,
        secondary_ids_path=secondary_ids_path,
        compounds_path=compounds_path,
    )


def parse_chebi_synonyms(
    input_file: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
    subset: str = "3star",
    *,
    names_path: Path | str | None = None,
    compounds_path: Path | str | None = None,
) -> Sec2PriMappingSet:
    """Parse ChEBI data files and extract synonym mappings.

    For releases >= 245, use names_path and compounds_path (TSV format).
    For releases < 245, use input_file (SDF format).

    Args:
        input_file: Path to the ChEBI SDF file (legacy format < 245).
        version: Version/release identifier for the datasource.
        show_progress: Whether to show progress bars during parsing.
        subset: "3star" or "complete" - which compounds to include.
        names_path: Path to names.tsv (new format >= 245).
        compounds_path: Path to compounds.tsv for 3-star filtering.

    Returns:
        Sec2PriMappingSet with computed cardinalities based on labels.
    """
    from pysec2pri.parsers import ChEBIParser

    parser = ChEBIParser(version=version, show_progress=show_progress, subset=subset)
    return parser.parse_synonyms(
        input_path=Path(input_file) if input_file else None,
        names_path=names_path,
        compounds_path=compounds_path,
    )


def parse_hmdb(
    input_file: Path | str,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Parse HMDB XML files and extract ID mappings."""
    from pysec2pri.parsers import HMDBParser

    parser = HMDBParser(version=version, show_progress=show_progress)
    return parser.parse(Path(input_file))


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
) -> Sec2PriMappingSet:
    """Parse HGNC complete set file and extract symbol mappings."""
    from pysec2pri.parsers import HGNCParser

    parser = HGNCParser(version=version, show_progress=show_progress)
    return parser.parse_symbols(Path(complete_set_file))


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
