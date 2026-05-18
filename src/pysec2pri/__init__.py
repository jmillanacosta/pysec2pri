"""Secondary to primary identifier mapping.

This package provides tools for converting secondary (retired/withdrawn)
identifiers to primary (current) identifiers across various biological
databases including ChEBI, HMDB, HGNC, NCBI Gene, and UniProt.
"""

from pysec2pri.api import (
    generate_chebi,
    generate_chebi_synonyms,
    generate_hgnc,
    generate_hgnc_symbols,
    generate_hmdb,
    generate_hmdb_proteins,
    generate_ncbi,
    generate_ncbi_symbols,
    generate_uniprot,
    generate_wikidata,
    resolve_ids,
    write_json,
    write_name2synonym,
    write_output,
    write_owl,
    write_rdf,
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
from pysec2pri.update_ids import build_lookup, update_ids

__all__ = [
    # Sentinel values
    "WITHDRAWN_ENTRY",
    "WITHDRAWN_ENTRY_LABEL",
    # Mapping Set classes
    "IdMappingSet",
    "LabelMappingSet",
    "Sec2PriMappingSet",
    # ID resolution
    "build_lookup",
    # generate_* (download + parse in one call)
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
    "update_ids",
    # Export functions
    "write_json",
    "write_name2synonym",
    "write_output",
    "write_owl",
    "write_rdf",
    "write_sec2pri",
    "write_sssom",
    "write_symbol2prev",
]
