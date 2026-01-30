"""Secondary to primary identifier mapping.

This package provides tools for converting secondary (retired/withdrawn)
identifiers to primary (current) identifiers across various biological
databases including ChEBI, HMDB, HGNC, NCBI Gene, UniProt, and Wikidata.
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

# Download utilities
from pysec2pri.download import (
    ReleaseInfo,
    check_release,
    download_datasource,
    download_file,
    get_latest_release_info,
)

# Logging utilities
from pysec2pri.logging import logger, set_log_level

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
    "ALL_DATASOURCES",
    "CHEBI",
    "HGNC",
    "HMDB",
    "MAPPING_JUSTIFICATION",
    "NCBI",
    "STANDARD_PREFIX_MAP",
    "UNIPROT",
    "WIKIDATA",
    # Base classes
    "BaseMapping",
    # Datasource-specific mapping classes
    "ChEBIMapping",
    # =========================================================================
    # Constants
    # =========================================================================
    "DatasourceConfig",
    "HGNCMapping",
    "HMDBMapping",
    "IdMapping",
    # =========================================================================
    # Models (for advanced usage)
    # =========================================================================
    # Enums
    "MappingCardinality",
    # =========================================================================
    # Diff utilities
    # =========================================================================
    "MappingDiff",
    # Container
    "MappingSet",
    "NCBIGeneMapping",
    "PredicateID",
    # =========================================================================
    # Download utilities
    # =========================================================================
    "ReleaseInfo",
    "SymbolMapping",
    "UniProtMapping",
    # Wikidata parser
    "WikidataParser",
    "check_release",
    # Utilities
    "compute_cardinality",
    "diff_mapping_sets",
    "diff_sssom_files",
    "download_datasource",
    "download_file",
    "get_datasource_config",
    "get_latest_release_info",
    "get_prefix_map_for_datasource",
    # =========================================================================
    # Logging utilities
    # =========================================================================
    "logger",
    # =========================================================================
    # High-level API functions (recommended for most users)
    # =========================================================================
    # Parsing functions
    "parse_chebi",
    "parse_hgnc",
    "parse_hmdb",
    "parse_ncbi",
    "parse_uniprot",
    "parse_wikidata",
    "set_log_level",
    "summarize_diff",
    "to_sssom_dataframe",
    "to_sssom_document",
    "to_sssom_mapping_set",
    "validate_identifier",
    # SSSOM output functions
    "write_sssom",
]
