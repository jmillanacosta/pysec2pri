"""Parsers for various biological database formats.

This package provides parsers for extracting secondary-to-primary identifier
mappings from various biological database formats.
"""

from pysec2pri.parsers.base import BaseParser
from pysec2pri.parsers.chebi import ChEBIParser
from pysec2pri.parsers.hgnc import HGNCParser
from pysec2pri.parsers.hmdb import HMDBParser
from pysec2pri.parsers.ncbi import NCBIParser
from pysec2pri.parsers.uniprot import UniProtParser
from pysec2pri.parsers.wikidata import WikidataParser, parse_wikidata

__all__ = [
    "BaseParser",
    "ChEBIParser",
    "HGNCParser",
    "HMDBParser",
    "NCBIParser",
    "UniProtParser",
    "WikidataParser",
    "parse_wikidata",
]
