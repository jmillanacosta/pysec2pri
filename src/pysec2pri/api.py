"""Main API functions for pysec2pri.

This module provides high-level functions for parsing biological database 
secondary-to-primary mapping files files and generating SSSOM output.
These functions can be imported directly from pysec2pri.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path
from typing import TYPE_CHECKING

from curies import Converter
from sssom import MappingSetDocument
from sssom_schema import Mapping as SSSOMMapping
from sssom_schema import MappingSet as SSSOMMapingSet

from pysec2pri.constants import MAPPING_JUSTIFICATION, STANDARD_PREFIX_MAP
from pysec2pri.models import MappingSet

# Re-export legacy output functions from exports module
from pysec2pri.exports import (
    write_sssom,
    write_primary_ids,
    write_sec2pri,
    write_name2synonym,
    write_symbol2prev,
    write_outputs,
)

if TYPE_CHECKING:
    from sssom.util import MappingSetDataFrame

__all__ = [
    # ChEBI
    "parse_chebi",
    # HMDB
    "parse_hmdb",
    # HGNC
    "parse_hgnc",
    # NCBI
    "parse_ncbi",
    # UniProt
    "parse_uniprot",
    # SSSOM conversion
    "to_sssom_mapping_set",
    "to_sssom_document",
    "to_sssom_dataframe",
    # Legacy output formats (re-exported from exports module)
    "write_sssom",
    "write_primary_ids",
    "write_sec2pri",
    "write_name2synonym",
    "write_symbol2prev",
    "write_outputs",
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
        >>> from pysec2pri import parse_chebi
        >>> mapping_set = parse_chebi("ChEBI_complete_3star.sdf")
        >>> print(f"Found {len(mapping_set)} mappings")
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
        >>> from pysec2pri import parse_hmdb
        >>> mapping_set = parse_hmdb("hmdb_metabolites.zip")
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
        include_unmapped_genes: If True, include entries for genes that
            have no alias or previous symbols (primary ID only).
            Default is False (only include actual mappings).

    Returns:
        MappingSet containing HGNC identifier and symbol mappings.

    Example:
        >>> from pysec2pri import parse_hgnc
        >>> mapping_set = parse_hgnc(
        ...     "withdrawn_2024-01-01.txt",
        ...     complete_set_file="hgnc_complete_set_2024-01-01.txt"
        ... )
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
        >>> from pysec2pri import parse_ncbi
        >>> mapping_set = parse_ncbi(
        ...     "gene_history.gz",
        ...     gene_info_file="gene_info.gz",
        ...     tax_id="9606"
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
        >>> from pysec2pri import parse_uniprot
        >>> mapping_set = parse_uniprot(
        ...     "sec_ac.txt",
        ...     delac_file="delac_sp.txt"
        ... )
    """
    from pysec2pri.parsers import UniProtParser

    parser = UniProtParser(version=version, show_progress=show_progress)
    return parser.parse(
        Path(sec_ac_file),
        delac_path=Path(delac_file) if delac_file else None,
    )


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
    if mapping_date is None:
        mapping_date = date.today().isoformat()

    # Convert individual mappings
    sssom_mappings = []
    for m in mapping_set.mappings:
        sssom_dict = m.to_sssom_dict()

        # Build SSSOM Mapping with required fields
        predicate_id = sssom_dict.get("predicate_id")
        if predicate_id is None:
            predicate_id = "skos:relatedMatch"  # fallback

        sssom_mapping = SSSOMMapping(
            subject_id=sssom_dict.get("subject_id"),
            subject_label=sssom_dict.get("subject_label"),
            predicate_id=predicate_id,
            object_id=sssom_dict.get("object_id"),
            object_label=sssom_dict.get("object_label"),
            mapping_justification=MAPPING_JUSTIFICATION,
            mapping_cardinality=sssom_dict.get("mapping_cardinality"),
            comment=sssom_dict.get("comment"),
        )
        sssom_mappings.append(sssom_mapping)

    # Build the MappingSet
    sssom_mapping_set = SSSOMMapingSet(
        mapping_set_id=mapping_set.mapping_set_id or "https://w3id.org/sssom/mappings",
        license=mapping_set.license_url,
        mappings=sssom_mappings,
        mapping_set_version=mapping_set.version,
        mapping_set_description=mapping_set.mapping_set_description,
        mapping_date=mapping_date,
        comment=mapping_set.comment,
    )

    return sssom_mapping_set


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
    sssom_ms = to_sssom_mapping_set(mapping_set, mapping_date)

    # Build prefix map for converter (start with standard prefixes)
    prefix_map = dict(STANDARD_PREFIX_MAP)

    # Add datasource-specific prefixes
    for prefix, url in mapping_set.curie_map.items():
        # Remove trailing colon if present
        clean_prefix = prefix.rstrip(":")
        prefix_map[clean_prefix] = url

    converter = Converter.from_prefix_map(prefix_map)

    return MappingSetDocument(
        mapping_set=sssom_ms,
        converter=converter,
    )


def to_sssom_dataframe(
    mapping_set: MappingSet,
    mapping_date: str | None = None,
) -> "MappingSetDataFrame":
    """Convert a pysec2pri MappingSet to an SSSOM MappingSetDataFrame.

    Args:
        mapping_set: The pysec2pri MappingSet to convert.
        mapping_date: Date string (YYYY-MM-DD). Defaults to today.

    Returns:
        An SSSOM MappingSetDataFrame.

    Raises:
        ImportError: If pandas is not installed.
    """
    from sssom.parsers import to_mapping_set_dataframe

    doc = to_sssom_document(mapping_set, mapping_date)
    return to_mapping_set_dataframe(doc)
