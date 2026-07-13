"""Parsers for biological database formats."""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

# name -> submodule providing it (resolved on first attribute access)
_LAZY_EXPORTS: dict[str, str] = {
    # Framework / base classes and config helpers
    "WITHDRAWN_ENTRY": "pysec2pri.parsers.base",
    "WITHDRAWN_ENTRY_LABEL": "pysec2pri.parsers.base",
    "BaseMappingSet": "pysec2pri.parsers.base",
    "BaseParser": "pysec2pri.parsers.base",
    "DatasourceConfig": "pysec2pri.parsers.base",
    "IdMappingSet": "pysec2pri.parsers.base",
    "LabelMappingSet": "pysec2pri.parsers.base",
    "get_datasource_config": "pysec2pri.parsers.base",
    "load_config": "pysec2pri.parsers.base",
    # Concrete datasource parsers
    "ChEBIParser": "pysec2pri.parsers.chebi",
    "EnsemblParser": "pysec2pri.parsers.ensembl",
    "HGNCParser": "pysec2pri.parsers.hgnc",
    "HMDBMetaboliteParser": "pysec2pri.parsers.hmdb",
    "HMDBProteinParser": "pysec2pri.parsers.hmdb",
    "NCBIParser": "pysec2pri.parsers.ncbi",
    "UniProtParser": "pysec2pri.parsers.uniprot",
    "VGNCParser": "pysec2pri.parsers.vgnc",
    "WikidataParser": "pysec2pri.parsers.wikidata",
}


def __getattr__(name: str) -> object:
    """Import and return a parser/base symbol on first access (PEP 562)."""
    module_name = _LAZY_EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(importlib.import_module(module_name), name)
    globals()[name] = value  # cache so subsequent access skips __getattr__
    return value


def __dir__() -> list[str]:
    """Expose lazily-available names to ``dir()`` and tab completion."""
    return sorted(__all__)


if TYPE_CHECKING:
    # Eager imports for type checkers / IDEs only; never executed at runtime.
    from pysec2pri.parsers.base import (
        WITHDRAWN_ENTRY,
        WITHDRAWN_ENTRY_LABEL,
        BaseMappingSet,
        BaseParser,
        DatasourceConfig,
        IdMappingSet,
        LabelMappingSet,
        get_datasource_config,
        load_config,
    )
    from pysec2pri.parsers.chebi import ChEBIParser
    from pysec2pri.parsers.ensembl import EnsemblParser
    from pysec2pri.parsers.hgnc import HGNCParser
    from pysec2pri.parsers.hmdb import HMDBMetaboliteParser, HMDBProteinParser
    from pysec2pri.parsers.ncbi import NCBIParser
    from pysec2pri.parsers.uniprot import UniProtParser
    from pysec2pri.parsers.vgnc import VGNCParser
    from pysec2pri.parsers.wikidata import WikidataParser

__all__ = [
    "WITHDRAWN_ENTRY",
    "WITHDRAWN_ENTRY_LABEL",
    "BaseMappingSet",
    "BaseParser",
    "ChEBIParser",
    "DatasourceConfig",
    "EnsemblParser",
    "HGNCParser",
    "HMDBMetaboliteParser",
    "HMDBProteinParser",
    "IdMappingSet",
    "LabelMappingSet",
    "NCBIParser",
    "UniProtParser",
    "VGNCParser",
    "WikidataParser",
    "get_datasource_config",
    "load_config",
]
