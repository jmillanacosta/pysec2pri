"""Secondary to primary identifier mapping.

This package provides tools for converting secondary (retired/withdrawn)
identifiers to primary (current) identifiers across various biological
databases including ChEBI, HMDB, HGNC, NCBI Gene, UniProt, and Wikidata.

Example usage:

    >>> from pysec2pri import parse_chebi, write_sssom
    >>> mapping_set = parse_chebi("ChEBI_complete_3star.sdf")
    >>> write_sssom(mapping_set, "chebi_mappings.sssom.tsv")
"""

# API functions - the main interface for users
from pysec2pri.api import (
    parse_chebi,
    parse_hgnc,
    parse_hmdb,
    parse_ncbi,
    parse_uniprot,
    to_sssom_dataframe,
    to_sssom_document,
    to_sssom_mapping_set,
    write_sssom,
)

# Constants
from pysec2pri.constants import (
    ALL_DATASOURCES,
    CHEBI,
    HGNC,
    HMDB,
    MAPPING_JUSTIFICATION,
    NCBI,
    STANDARD_PREFIX_MAP,
    UNIPROT,
    WIKIDATA,
    DatasourceConfig,
    get_datasource_config,
    get_prefix_map_for_datasource,
    validate_identifier,
)

# Diff utilities
from pysec2pri.diff import (
    MappingDiff,
    diff_mapping_sets,
    diff_sssom_files,
    summarize_diff,
)

# Logging utilities
from pysec2pri.logging import logger, set_log_level

# Download utilities
from pysec2pri.download import (
    ReleaseInfo,
    check_release,
    download_datasource,
    download_file,
    get_latest_release_info,
)

# Models
from pysec2pri.models import (
    BaseMapping,
    ChEBIMapping,
    HGNCMapping,
    HMDBMapping,
    IdMapping,
    MappingCardinality,
    MappingSet,
    NCBIGeneMapping,
    PredicateID,
    SymbolMapping,
    UniProtMapping,
    compute_cardinality,
)

# Wikidata
from pysec2pri.parsers.wikidata import WikidataParser, parse_wikidata

__all__ = [
    # =========================================================================
    # High-level API functions (recommended for most users)
    # =========================================================================
    # Parsing functions
    "parse_chebi",
    "parse_hmdb",
    "parse_hgnc",
    "parse_ncbi",
    "parse_uniprot",
    "parse_wikidata",
    # SSSOM output functions
    "write_sssom",
    "to_sssom_mapping_set",
    "to_sssom_document",
    "to_sssom_dataframe",
    # =========================================================================
    # Download utilities
    # =========================================================================
    "ReleaseInfo",
    "check_release",
    "download_datasource",
    "download_file",
    "get_latest_release_info",
    # =========================================================================
    # Diff utilities
    # =========================================================================
    "MappingDiff",
    "diff_mapping_sets",
    "diff_sssom_files",
    "summarize_diff",
    # =========================================================================
    # Constants
    # =========================================================================
    "DatasourceConfig",
    "CHEBI",
    "HMDB",
    "HGNC",
    "NCBI",
    "UNIPROT",
    "WIKIDATA",
    "ALL_DATASOURCES",
    "STANDARD_PREFIX_MAP",
    "MAPPING_JUSTIFICATION",
    "get_datasource_config",
    "get_prefix_map_for_datasource",
    "validate_identifier",
    # =========================================================================
    # Models (for advanced usage)
    # =========================================================================
    # Enums
    "MappingCardinality",
    "PredicateID",
    # Base classes
    "BaseMapping",
    "IdMapping",
    "SymbolMapping",
    # Datasource-specific mapping classes
    "ChEBIMapping",
    "HMDBMapping",
    "UniProtMapping",
    "HGNCMapping",
    "NCBIGeneMapping",
    # Container
    "MappingSet",
    # Utilities
    "compute_cardinality",
    # Wikidata parser
    "WikidataParser",
    # =========================================================================
    # Logging utilities
    # =========================================================================
    "logger",
    "set_log_level",
]
