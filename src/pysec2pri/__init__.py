"""Secondary to primary identifier mapping.

This package provides tools for converting secondary (retired/withdrawn)
identifiers to primary (current) identifiers across various biological
databases including ChEBI, HMDB, HGNC, NCBI Gene, and UniProt.
"""

from pysec2pri.api import (
    parse_chebi,
    parse_chebi_synonyms,
    parse_hgnc,
    parse_hgnc_symbols,
    parse_hmdb,
    parse_ncbi,
    parse_ncbi_symbols,
    parse_uniprot,
    write_name2synonym,
    write_sec2pri,
    write_sssom,
    write_symbol2prev,
)
from pysec2pri.parsers.base import (
    WITHDRAWN_ENTRY,
    WITHDRAWN_ENTRY_LABEL,
    IdMappingSet,
    LabelMappingSet,
    Sec2PriMappingSet,
)

__all__ = [
    # Sentinel values
    "WITHDRAWN_ENTRY",
    "WITHDRAWN_ENTRY_LABEL",
    # Mapping Set classes
    "IdMappingSet",
    "LabelMappingSet",
    "Sec2PriMappingSet",
    # Parsing functions
    "parse_chebi",
    "parse_chebi_synonyms",
    "parse_hgnc",
    "parse_hgnc_symbols",
    "parse_hmdb",
    "parse_ncbi",
    "parse_ncbi_symbols",
    "parse_uniprot",
    # Export functions
    "write_name2synonym",
    "write_sec2pri",
    "write_sssom",
    "write_symbol2prev",
]
