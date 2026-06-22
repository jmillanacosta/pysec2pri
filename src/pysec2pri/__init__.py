"""Secondary to primary identifier mapping.

This package provides tools for converting secondary (retired/withdrawn)
identifiers to primary (current) identifiers across various biological
databases including ChEBI, HMDB, HGNC, NCBI Gene, and UniProt.
"""

from pysec2pri.api import (
    crosswalk,
    find_ambiguous,
    generate_chebi,
    generate_chebi_synonyms,
    generate_hgnc,
    generate_hgnc_labels,
    generate_hgnc_primary_ids,
    generate_hmdb,
    generate_hmdb_proteins,
    generate_ncbi,
    generate_ncbi_labels,
    generate_uniprot,
    generate_wikidata,
    generate_wikidata_labels,
    list_versions,
    resolve_ids,
    resolve_labels,
    write_json,
    write_label_sec2pri,
    write_name2synonym,
    write_output,
    write_owl,
    write_rdf,
    write_sec2pri,
    write_sssom,
)
from pysec2pri.context import (
    ContextSpec,
    DecisionRecord,
    XrefMapping,
    XrefRecord,
    load_xref_mapping,
)
from pysec2pri.parsers.base import (
    WITHDRAWN_ENTRY,
    WITHDRAWN_ENTRY_LABEL,
    AmbiguousMappingSet,
    IdMappingSet,
    LabelMappingSet,
    Sec2PriMappingSet,
)
from pysec2pri.update_ids import (
    build_alias_index,
    build_label_lookup,
    build_lookup,
    build_primary_token_to_id,
    resolve_ambiguous_with_hints,
    update_ids,
    update_labels,
)

__all__ = [
    # Sentinel values
    "WITHDRAWN_ENTRY",
    "WITHDRAWN_ENTRY_LABEL",
    # Mapping Set classes
    "AmbiguousMappingSet",
    # Disambiguation context (symbol/id/xref)
    "ContextSpec",
    "DecisionRecord",
    "IdMappingSet",
    "LabelMappingSet",
    "Sec2PriMappingSet",
    "XrefMapping",
    "XrefRecord",
    # ID resolution
    "build_alias_index",
    "build_label_lookup",
    "build_lookup",
    "build_primary_token_to_id",
    "crosswalk",
    "find_ambiguous",
    # generate_* (download + parse in one call)
    "generate_chebi",
    "generate_chebi_synonyms",
    "generate_hgnc",
    "generate_hgnc_labels",
    "generate_hgnc_primary_ids",
    "generate_hmdb",
    "generate_hmdb_proteins",
    "generate_ncbi",
    "generate_ncbi_labels",
    "generate_uniprot",
    "generate_wikidata",
    "generate_wikidata_labels",
    "list_versions",
    "load_xref_mapping",
    "resolve_ambiguous_with_hints",
    "resolve_ids",
    "resolve_labels",
    "update_ids",
    "update_labels",
    # Export functions
    "write_json",
    "write_label_sec2pri",
    "write_name2synonym",
    "write_output",
    "write_owl",
    "write_rdf",
    "write_sec2pri",
    "write_sssom",
]
