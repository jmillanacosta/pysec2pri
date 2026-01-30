"""Main API functions for pysec2pri.

This module provides high-level functions for parsing biological database
secondary-to-primary mapping files files and generating SSSOM output.
These functions can be imported directly from pysec2pri.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

# Re-export legacy output functions from exports module
from pysec2pri.exports import (
    write_name2synonym,
    write_outputs,
    write_sec2pri,
    write_sssom,
    write_subject_ids,
    write_symbol2prev,
)
from pysec2pri.models import MappingSet

if TYPE_CHECKING:
    from sssom import MappingSetDocument
    from sssom.util import MappingSetDataFrame
    from sssom_schema import MappingSet as SSSOMMapingSet

__all__ = [
    # ChEBI
    "parse_chebi",
    # HGNC
    "parse_hgnc",
    # HMDB
    "parse_hmdb",
    # NCBI
    "parse_ncbi",
    # UniProt
    "parse_uniprot",
    "to_sssom_dataframe",
    "to_sssom_document",
    # SSSOM conversion
    "to_sssom_mapping_set",
    "write_name2synonym",
    "write_outputs",
    "write_sec2pri",
    # Legacy output formats (re-exported from exports module)
    "write_sssom",
    "write_subject_ids",
    "write_symbol2prev",
]


# =============================================================================
# ChEBI Functions
# =============================================================================


def parse_chebi(
    input_file: Path | str,
    version: str | None = None,
    show_progress: bool = True,
) -> MappingSet:
    """Parse a ChEBI SDF file and extract mappings.

    Args:
        input_file: Path to ChEBI SDF file (e.g., ChEBI_complete_3star.sdf).
        version: Version/release of the source database.
        show_progress: Whether to show progress bars during parsing.

    Returns:
        MappingSet containing ChEBI identifier and synonym mappings.

    Example:
        >>> from pysec2pri import parse_chebi  # doctest: +SKIP
        >>> mapping_set = parse_chebi("ChEBI_complete_3star.sdf")  # doctest: +SKIP
        >>> print(f"Found {len(mapping_set)} mappings")  # doctest: +SKIP
    """
    from pysec2pri.parsers import ChEBIParser

    parser = ChEBIParser(version=version, show_progress=show_progress)
    return parser.parse(Path(input_file))


# =============================================================================
# HMDB Functions
# =============================================================================


def parse_hmdb(
    input_file: Path | str,
    version: str | None = None,
    show_progress: bool = True,
) -> MappingSet:
    """Parse HMDB XML files and extract mappings.

    Args:
        input_file: Path to HMDB ZIP file or directory with XML files.
        version: Version/release of the source database.
        show_progress: Whether to show progress bars during parsing.

    Returns:
        MappingSet containing HMDB identifier and synonym mappings.

    Example:
        >>> from pysec2pri import parse_hmdb  # doctest: +SKIP
        >>> mapping_set = parse_hmdb("hmdb_metabolites.zip")  # doctest: +SKIP
    """
    from pysec2pri.parsers import HMDBParser

    parser = HMDBParser(version=version, show_progress=show_progress)
    return parser.parse(Path(input_file))


# =============================================================================
# HGNC Functions
# =============================================================================


def parse_hgnc(
    withdrawn_file: Path | str,
    complete_set_file: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
    include_unmapped_genes: bool = False,
) -> MappingSet:
    """Parse HGNC files and extract mappings.

    Args:
        withdrawn_file: Path to HGNC withdrawn IDs file.
        complete_set_file: Path to HGNC complete set file (optional).
        version: Version/release of the source database.
        show_progress: Whether to show progress bars during parsing.
        include_unmapped_genes: If True, include entries for priID only genes.

    Returns:
        MappingSet containing HGNC identifier and symbol mappings.
    """
    from pysec2pri.parsers import HGNCParser

    parser = HGNCParser(
        version=version,
        show_progress=show_progress,
        include_unmapped_genes=include_unmapped_genes,
    )
    return parser.parse(
        Path(withdrawn_file),
        complete_set_path=Path(complete_set_file) if complete_set_file else None,
    )


# =============================================================================
# NCBI Functions
# =============================================================================


def parse_ncbi(
    history_file: Path | str,
    gene_info_file: Path | str | None = None,
    tax_id: str = "9606",
    version: str | None = None,
    show_progress: bool = True,
) -> MappingSet:
    """Parse NCBI Gene files and extract mappings.

    Args:
        history_file: Path to gene_history file (can be .gz compressed).
        gene_info_file: Path to gene_info file for symbol information.
        tax_id: Taxonomy ID to filter (default: "9606" for human).
        version: Version/release of the source database.
        show_progress: Whether to show progress bars during parsing.

    Returns:
        MappingSet containing NCBI Gene identifier and symbol mappings.

    Example:
        >>> from pysec2pri import parse_ncbi  # doctest: +SKIP
        >>> mapping_set = parse_ncbi(
        ...     "gene_history.gz",
        ...     gene_info_file="gene_info.gz",
        ...     tax_id="9606",  # doctest: +SKIP
        ... )
    """
    from pysec2pri.parsers import NCBIParser

    parser = NCBIParser(version=version, show_progress=show_progress)
    return parser.parse(
        Path(history_file),
        gene_info_path=Path(gene_info_file) if gene_info_file else None,
        tax_id=tax_id,
    )


# =============================================================================
# UniProt Functions
# =============================================================================


def parse_uniprot(
    sec_ac_file: Path | str,
    delac_file: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> MappingSet:
    """Parse UniProt files and extract mappings.

    Args:
        sec_ac_file: Path to sec_ac.txt (secondary accessions file).
        delac_file: Path to delac_sp.txt (deleted accessions file).
        version: Version/release of the source database.
        show_progress: Whether to show progress bars during parsing.

    Returns:
        MappingSet containing UniProt identifier mappings.

    Example:
        >>> from pysec2pri import parse_uniprot  # doctest: +SKIP
        >>> mapping_set = parse_uniprot("sec_ac.txt", delac_file="delac_sp.txt")  # doctest: +SKIP
    """
    from pysec2pri.parsers import UniProtParser

    parser = UniProtParser(version=version, show_progress=show_progress)
    return parser.parse(
        Path(sec_ac_file),
        delac_path=Path(delac_file) if delac_file else None,
    )


# =============================================================================
# Wikidata Functions
# =============================================================================


def parse_wikidata(
    entity_type: str = "metabolites",
    version: str | None = None,
    endpoint: str | None = None,
    show_progress: bool = True,
    test_subset: bool = False,
) -> MappingSet:
    """Parse Wikidata redirects for a specific entity type, with test subset support."""
    from pysec2pri.parsers import WikidataParser

    parser = WikidataParser(
        version=version,
        show_progress=show_progress,
        entity_type=entity_type,
        endpoint=endpoint,
        test_subset=test_subset,
    )
    return parser.parse(None)


# =============================================================================
# SSSOM Conversion Functions
# =============================================================================


def to_sssom_mapping_set(
    mapping_set: MappingSet,
    mapping_date: str | None = None,
) -> SSSOMMapingSet:
    """Convert a pysec2pri MappingSet to an sssom_schema MappingSet.

    Args:
        mapping_set: The pysec2pri MappingSet to convert.
        mapping_date: Date string (YYYY-MM-DD). Defaults to today.

    Returns:
        An sssom_schema MappingSet object.
    """
    return mapping_set.to_sssom_mapping_set(mapping_date)


def to_sssom_document(
    mapping_set: MappingSet,
    mapping_date: str | None = None,
) -> MappingSetDocument:
    """Convert a pysec2pri MappingSet to an SSSOM MappingSetDocument.

    Args:
        mapping_set: The pysec2pri MappingSet to convert.
        mapping_date: Date string (YYYY-MM-DD). Defaults to today.

    Returns:
        An SSSOM MappingSetDocument with converter.
    """
    return mapping_set.to_sssom_document(mapping_date)


def to_sssom_dataframe(
    mapping_set: MappingSet,
    mapping_date: str | None = None,
) -> MappingSetDataFrame:
    """Convert a pysec2pri MappingSet to an SSSOM MappingSetDataFrame.

    Args:
        mapping_set: The pysec2pri MappingSet to convert.
        mapping_date: Date string (YYYY-MM-DD). Defaults to today.

    Returns:
        An SSSOM MappingSetDataFrame.

    Raises:
        ImportError: If pandas is not installed.
    """
    return mapping_set.to_sssom_dataframe(mapping_date)
