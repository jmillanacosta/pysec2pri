"""Constants for supported datasources."""

from pysec2pri.parsers.base import DatasourceConfig, get_datasource_config

__all__ = [
    "ALL_DATASOURCES",
    "CHEBI",
    "HGNC",
    "HMDB",
    "NCBI",
    "UNIPROT",
    "WIKIDATA",
]

# Pre-loaded datasource configurations
CHEBI: DatasourceConfig = get_datasource_config("chebi")
HGNC: DatasourceConfig = get_datasource_config("hgnc")
HMDB: DatasourceConfig = get_datasource_config("hmdb")
NCBI: DatasourceConfig = get_datasource_config("ncbi")
UNIPROT: DatasourceConfig = get_datasource_config("uniprot")
WIKIDATA: DatasourceConfig = get_datasource_config("wikidata")

# Mapping of datasource names to their configurations
ALL_DATASOURCES: dict[str, DatasourceConfig] = {
    "chebi": CHEBI,
    "hgnc": HGNC,
    "hmdb": HMDB,
    "ncbi": NCBI,
    "uniprot": UNIPROT,
    "wikidata": WIKIDATA,
}
