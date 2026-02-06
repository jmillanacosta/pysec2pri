"""Main API functions for pysec2pri.

This module provides high-level functions for parsing biological database
secondary-to-primary mapping files and generating output.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from pysec2pri.exports import (
    write_name2synonym,
    write_sec2pri,
    write_sssom,
    write_symbol2prev,
)

if TYPE_CHECKING:
    from pysec2pri.parsers.base import Sec2PriMappingSet

__all__ = [
    "parse_chebi",
    "parse_chebi_synonyms",
    "parse_hgnc",
    "parse_hgnc_symbols",
    "parse_hmdb",
    "parse_ncbi",
    "parse_ncbi_symbols",
    "parse_uniprot",
    "parse_wikidata",
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

    parser = ChEBIParser(
        version=version, show_progress=show_progress, subset=subset
    )
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

    parser = ChEBIParser(
        version=version, show_progress=show_progress, subset=subset
    )
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
