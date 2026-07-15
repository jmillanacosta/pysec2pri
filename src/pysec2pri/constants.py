"""Constants for supported datasources.

Datasource configs are built on access.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from mapkgsutils.parsers.config import get_datasource_config

if TYPE_CHECKING:
    from mapkgsutils.parsers.config import DatasourceConfig

__all__ = [
    "ALL_DATASOURCES",
    "CHEBI",
    "CONSOLIDATE_DATASOURCES",
    "ENSEMBL",
    "HGNC",
    "NCBI",
    "UNIPROT",
    "VGNC",
    "WIKIDATA",
]

_CONFIG_PACKAGE = "pysec2pri.config"

#: Datasource config ids that support ``consolidate`` (first-seen-date).
CONSOLIDATE_DATASOURCES: tuple[str, ...] = ("chebi", "ensembl", "hgnc", "ncbi", "uniprot", "vgnc")

# Public constant name is (``<key>.yaml``).
_DATASOURCE_KEYS: dict[str, str] = {
    "CHEBI": "chebi",
    "ENSEMBL": "ensembl",
    "HGNC": "hgnc",
    "HMDB_PROT": "hmdb_proteins",
    "HMDB_MET": "hmdb_metabolites",
    "NCBI": "ncbi",
    "UNIPROT": "uniprot",
    "VGNC": "vgnc",
    "WIKIDATA": "wikidata",
}


def __getattr__(name: str) -> object:
    """Build a datasource config (or the ``ALL_DATASOURCES`` map) on first access."""
    if name == "ALL_DATASOURCES":
        value: object = {
            key: get_datasource_config(key, config_package=_CONFIG_PACKAGE)
            for key in _DATASOURCE_KEYS.values()
        }
    elif name in _DATASOURCE_KEYS:
        value = get_datasource_config(_DATASOURCE_KEYS[name], config_package=_CONFIG_PACKAGE)
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    globals()[name] = value  # cache
    return value


def __dir__() -> list[str]:
    """Expose lazily-available names to ``dir()`` and tab completion."""
    return sorted({*__all__, *_DATASOURCE_KEYS})


if TYPE_CHECKING:
    # For type checking.
    CHEBI: DatasourceConfig
    ENSEMBL: DatasourceConfig
    HGNC: DatasourceConfig
    HMDB_PROT: DatasourceConfig
    HMDB_MET: DatasourceConfig
    NCBI: DatasourceConfig
    UNIPROT: DatasourceConfig
    VGNC: DatasourceConfig
    WIKIDATA: DatasourceConfig
    ALL_DATASOURCES: dict[str, DatasourceConfig]
