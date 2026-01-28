"""Centralized constants and configuration for pysec2pri.

This module contains all configuration values including CURIE prefixes,
URL mappings, output file naming conventions, download URLs, and other
constants used throughout the package.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum


class OutputType(str, Enum):
    """Available output types for pysec2pri.

    These correspond to the outputs from the legacy R scripts:
    - SSSOM: Standard SSSOM TSV format (default)
    - PRIMARY_IDS: List of primary IDs
    - SEC2PRI: Secondary to primary ID mappings (legacy format)
    - NAME2SYNONYM: Name to synonym mappings (for ChEBI, HMDB)
    - SYMBOL2PREV: Symbol to alias/previous symbol mappings (for HGNC, NCBI)
    """

    SSSOM = "sssom"
    PRIMARY_IDS = "priIDs"
    SEC2PRI = "secID2priID"
    NAME2SYNONYM = "name2synonym"
    SYMBOL2PREV = "symbol2prev"


@dataclass(frozen=True)
class DatasourceConfig:
    """Configuration for a biological database datasource."""

    name: str
    prefix: str
    curie_base_url: str
    default_output_filename: str
    # Available output types for this datasource
    available_outputs: tuple[OutputType, ...] = field(default_factory=tuple)
    # Download URLs for source files
    download_urls: dict[str, str] = field(default_factory=dict)
    # Regex pattern for validating identifiers
    id_pattern: str = ""
    # Archive/index URL for checking new releases
    archive_url: str = ""


# =============================================================================
# Datasource configurations
# =============================================================================

CHEBI = DatasourceConfig(
    name="ChEBI",
    prefix="CHEBI",
    curie_base_url="http://purl.obolibrary.org/obo/CHEBI_",
    default_output_filename="chebi_sec2pri.sssom.tsv",
    available_outputs=(
        OutputType.SSSOM,
        OutputType.PRIMARY_IDS,
        OutputType.SEC2PRI,
        OutputType.NAME2SYNONYM,
    ),
    download_urls={
        "sdf": "http://ftp.ebi.ac.uk/pub/databases/chebi/SDF/chebi_3_stars.sdf.gz",  # noqa: E501
    },
    archive_url="http://ftp.ebi.ac.uk/pub/databases/chebi/archive/",
    id_pattern=r"^\d+$",
)

HMDB = DatasourceConfig(
    name="HMDB",
    prefix="HMDB",
    curie_base_url="http://www.hmdb.ca/metabolites/",
    default_output_filename="hmdb_sec2pri.sssom.tsv",
    available_outputs=(
        OutputType.SSSOM,
        OutputType.PRIMARY_IDS,
        OutputType.SEC2PRI,
        OutputType.NAME2SYNONYM,
    ),
    download_urls={
        "metabolites": "http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip",  # noqa: E501
    },
    id_pattern=r"^HMDB\d+$",
)

HGNC = DatasourceConfig(
    name="HGNC",
    prefix="HGNC",
    curie_base_url="http://identifiers.org/hgnc/",
    default_output_filename="hgnc_sec2pri.sssom.tsv",
    available_outputs=(
        OutputType.SSSOM,
        OutputType.PRIMARY_IDS,
        OutputType.SEC2PRI,
        OutputType.SYMBOL2PREV,
    ),
    download_urls={  # DELETE
        "complete": "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt",  # noqa: E501
        "withdrawn": "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/withdrawn.tsv",  # DELETE # noqa: E501
    },
    archive_url="https://www.genenames.org/download/archive/quarterly/tsv/",
    id_pattern=r"^HGNC:\d+$",
)

NCBI = DatasourceConfig(
    name="NCBI",
    prefix="NCBIGene",
    curie_base_url="http://www.ncbi.nlm.nih.gov/gene/",
    default_output_filename="ncbi_sec2pri.sssom.tsv",
    available_outputs=(
        OutputType.SSSOM,
        OutputType.PRIMARY_IDS,
        OutputType.SEC2PRI,
        OutputType.SYMBOL2PREV,
    ),
    download_urls={
        "gene_history": "https://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz",
        "gene_info": "https://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz",
    },
    id_pattern=r"^\d+$",
)

UNIPROT = DatasourceConfig(
    name="UniProt",
    prefix="UniProt",
    curie_base_url="http://identifiers.org/uniprot/",
    default_output_filename="uniprot_sec2pri.sssom.tsv",
    available_outputs=(
        OutputType.SSSOM,
        OutputType.PRIMARY_IDS,
        OutputType.SEC2PRI,
    ),
    download_urls={
        "secondary": "https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/docs/sec_ac.txt",  # noqa: E501
        "deleted": "https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/docs/delac_sp.txt",  # noqa: E501
    },
    archive_url="https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/",  # noqa: E501
    id_pattern=r"^[A-Z][A-Z0-9]{5}$|^[A-Z][A-Z0-9]{9}$",
)

# Wikidata configuration (special - uses SPARQL queries)
WIKIDATA = DatasourceConfig(
    name="Wikidata",
    prefix="wikidata",
    curie_base_url="http://www.wikidata.org/entity/",
    default_output_filename="wikidata_sec2pri.sssom.tsv",
    available_outputs=(
        OutputType.SSSOM,
        OutputType.SEC2PRI,
    ),
    download_urls={
        "qlever_endpoint": "https://qlever.dev/api/wikidata",
    },
    id_pattern=r"^Q\d+$",
)

# =============================================================================
# sec2pri project namespace
# =============================================================================

SEC2PRI_NAMESPACE = "https://github.com/sec2pri/sec2pri/vocab/"
SEC2PRI_PREFIX = "sec2pri"

# CURIE for withdrawn/deleted entries with no replacement
WITHDRAWN_CURIE = f"{SEC2PRI_PREFIX}:WithdrawnEntry"

# =============================================================================
# Standard prefix map for SSSOM output
# =============================================================================

STANDARD_PREFIX_MAP: dict[str, str] = {
    "IAO": "http://purl.obolibrary.org/obo/IAO_",
    "oboInOwl": "http://www.geneontology.org/formats/oboInOwl#",
    "semapv": "https://w3id.org/semapv/vocab/",
    "skos": "http://www.w3.org/2004/02/skos/core#",
    "sssom": "https://w3id.org/sssom/",
    SEC2PRI_PREFIX: SEC2PRI_NAMESPACE,
}

# =============================================================================
# Mapping justification URIs
# =============================================================================

MAPPING_JUSTIFICATION = "semapv:BackgroundKnowledgeBasedMatching"

# =============================================================================
# Default input file patterns (for reference)
# =============================================================================

DEFAULT_INPUT_PATTERNS: dict[str, list[str]] = {
    "chebi": ["ChEBI_complete.sdf", "ChEBI_complete_3star.sdf"],
    "hmdb": ["hmdb_metabolites.zip", "hmdb_metabolites.xml"],
    "hgnc_withdrawn": ["withdrawn.txt"],
    "hgnc_complete": ["hgnc_complete_set.txt"],
    "ncbi_history": ["gene_history.gz", "gene_history"],
    "ncbi_info": ["gene_info.gz", "gene_info"],
    "uniprot_secondary": ["sec_ac.txt"],
    "uniprot_deleted": ["delac_sp.txt"],
}

# =============================================================================
# Lookup helpers
# =============================================================================

ALL_DATASOURCES: dict[str, DatasourceConfig] = {
    "chebi": CHEBI,
    "hmdb": HMDB,
    "hgnc": HGNC,
    "ncbi": NCBI,
    "uniprot": UNIPROT,
    "wikidata": WIKIDATA,
}


def get_datasource_config(name: str) -> DatasourceConfig | None:
    """Get datasource configuration by name (case-insensitive).

    Args:
        name: Datasource name (e.g., "chebi", "HMDB", "hgnc").

    Returns:
        DatasourceConfig if found, None otherwise.
    """
    return ALL_DATASOURCES.get(name.lower())


def get_prefix_map_for_datasource(config: DatasourceConfig) -> dict[str, str]:
    """Get a complete prefix map including standard and datasource prefixes.

    Args:
        config: The datasource configuration.

    Returns:
        Dictionary mapping prefixes to URLs.
    """
    return {**STANDARD_PREFIX_MAP, config.prefix: config.curie_base_url}


def validate_identifier(identifier: str, config: DatasourceConfig) -> bool:
    """Validate an identifier against the datasource's ID pattern.

    Args:
        identifier: The identifier to validate.
        config: The datasource configuration with the pattern.

    Returns:
        True if the identifier matches the pattern.
    """
    if not config.id_pattern:
        return True
    return bool(re.match(config.id_pattern, identifier))


def get_output_filename(
    config: DatasourceConfig,
    output_type: OutputType,
) -> str:
    """Get the output filename for a specific output type.

    Follows the legacy naming convention: {database}_{output_type}.tsv

    Args:
        config: The datasource configuration.
        output_type: The type of output file.

    Returns:
        The filename string.
    """
    db_name = config.name
    if output_type == OutputType.SSSOM:
        return f"{db_name}_sec2pri.sssom.tsv"
    return f"{db_name}_{output_type.value}.tsv"


def get_output_directory_name(datasource_name: str) -> str:
    """Get the default output directory name with date suffix.

    Format: {datasource}{ddmmyy}

    Args:
        datasource_name: Name of the datasource.

    Returns:
        Directory name string (e.g., "chebi280126").
    """
    from datetime import date as date_module

    today = date_module.today()
    date_suffix = today.strftime("%d%m%y")
    return f"{datasource_name.lower()}{date_suffix}"
