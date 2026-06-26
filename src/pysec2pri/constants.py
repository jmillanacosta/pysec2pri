"""Constants for supported datasources."""

from pysec2pri.parsers.base import DatasourceConfig, get_datasource_config

__all__ = [
    "ALL_DATASOURCES",
    "CHEBI",
    "ENSEMBL",
    "HGNC",
    "NCBI",
    "UNIPROT",
    "VGNC",
    "WIKIDATA",
]

# Pre-loaded datasource configurations
_CONFIG_PACKAGE = "pysec2pri.config"
CHEBI: DatasourceConfig = get_datasource_config("chebi", config_package=_CONFIG_PACKAGE)
ENSEMBL: DatasourceConfig = get_datasource_config("ensembl", config_package=_CONFIG_PACKAGE)
HGNC: DatasourceConfig = get_datasource_config("hgnc", config_package=_CONFIG_PACKAGE)
HMDB_PROT: DatasourceConfig = get_datasource_config("hmdb_proteins", config_package=_CONFIG_PACKAGE)
HMDB_MET: DatasourceConfig = get_datasource_config(
    "hmdb_metabolites", config_package=_CONFIG_PACKAGE
)
NCBI: DatasourceConfig = get_datasource_config("ncbi", config_package=_CONFIG_PACKAGE)
UNIPROT: DatasourceConfig = get_datasource_config("uniprot", config_package=_CONFIG_PACKAGE)
VGNC: DatasourceConfig = get_datasource_config("vgnc", config_package=_CONFIG_PACKAGE)
WIKIDATA: DatasourceConfig = get_datasource_config("wikidata", config_package=_CONFIG_PACKAGE)

# Mapping of datasource names to their configurations
ALL_DATASOURCES: dict[str, DatasourceConfig] = {
    "chebi": CHEBI,
    "ensembl": ENSEMBL,
    "hgnc": HGNC,
    "hmdb_metabolites": HMDB_MET,
    "hmdb_proteins": HMDB_PROT,
    "ncbi": NCBI,
    "uniprot": UNIPROT,
    "vgnc": VGNC,
    "wikidata": WIKIDATA,
}
