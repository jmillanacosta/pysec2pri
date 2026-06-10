"""Constants for supported datasources."""

from pysec2pri.parsers.base import DatasourceConfig, get_datasource_config

__all__ = [
    "ALL_DATASOURCES",
    "CHEBI",
    "HGNC",
    "NCBI",
    "UNIPROT",
    "WIKIDATA",
]

# Pre-loaded datasource configurations
CHEBI: DatasourceConfig = get_datasource_config("chebi")
HGNC: DatasourceConfig = get_datasource_config("hgnc")
HMDB_PROT: DatasourceConfig = get_datasource_config("hmdb_proteins")
HMDB_MET: DatasourceConfig = get_datasource_config("hmdb_metabolites")
NCBI: DatasourceConfig = get_datasource_config("ncbi")
UNIPROT: DatasourceConfig = get_datasource_config("uniprot")
WIKIDATA: DatasourceConfig = get_datasource_config("wikidata")

# Mapping of datasource names to their configurations
ALL_DATASOURCES: dict[str, DatasourceConfig] = {
    "chebi": CHEBI,
    "hgnc": HGNC,
    "hmdb_metabolites": HMDB_MET,
    "hmdb_proteins": HMDB_PROT,
    "ncbi": NCBI,
    "uniprot": UNIPROT,
    "wikidata": WIKIDATA,
}
