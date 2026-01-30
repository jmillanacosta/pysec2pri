"""SPARQL queries for Wikidata redirect mappings.

This module loads and provides access to SPARQL queries stored as .rq files.
The queries are sourced from the sec2pri/mapping_preprocessing repository.
"""

from __future__ import annotations

from pathlib import Path

__all__ = [
    "CHEMICAL_REDIRECTS_QUERY",
    "GENE_REDIRECTS_QUERY",
    "PROTEIN_REDIRECTS_QUERY",
    "WIKIDATA_QUERIES",
    "WIKIDATA_TEST_QUERIES",
    "get_query",
]

# Directory containing .rq query files
_QUERIES_DIR = Path(__file__).parent


def _load_query(filename: str) -> str:
    """Load a SPARQL query from a .rq file."""
    query_path = _QUERIES_DIR / filename
    if not query_path.exists():
        raise FileNotFoundError(f"Query file not found: {query_path}")
    return query_path.read_text(encoding="utf-8")


# Load queries from .rq files
CHEMICAL_REDIRECTS_QUERY = _load_query("chemical_redirects.rq")
GENE_REDIRECTS_QUERY = _load_query("gene_redirects.rq")
PROTEIN_REDIRECTS_QUERY = _load_query("protein_redirects.rq")
CHEMICAL_REDIRECTS_TEST_QUERY = _load_query("chemical_redirects_test.rq")
GENE_REDIRECTS_TEST_QUERY = _load_query("gene_redirects_test.rq")
PROTEIN_REDIRECTS_TEST_QUERY = _load_query("protein_redirects_test.rq")
# Query mapping by entity type
WIKIDATA_QUERIES: dict[str, str] = {
    "metabolites": CHEMICAL_REDIRECTS_QUERY,
    "chemicals": CHEMICAL_REDIRECTS_QUERY,
    "genes": GENE_REDIRECTS_QUERY,
    "proteins": PROTEIN_REDIRECTS_QUERY,
}

WIKIDATA_TEST_QUERIES: dict[str, str] = {
    "metabolites": CHEMICAL_REDIRECTS_TEST_QUERY,
    "chemicals": CHEMICAL_REDIRECTS_TEST_QUERY,
    "genes": GENE_REDIRECTS_TEST_QUERY,
    "proteins": PROTEIN_REDIRECTS_TEST_QUERY,
}
# Column name mapping for each query type
# Maps entity type to the expected column names in query results
QUERY_COLUMNS: dict[str, dict[str, str]] = {
    "metabolites": {
        "subject_id": "primaryID",
        "object_id": "secondaryID",
        "primary_label": "name",
        "secondary_label": "synonym",
    },
    "chemicals": {
        "subject_id": "primaryID",
        "object_id": "secondaryID",
        "primary_label": "name",
        "secondary_label": "synonym",
    },
    "genes": {
        "subject_id": "primaryID",
        "object_id": "secondaryID",
        "primary_label": "primarySymbol",
        "secondary_label": "secondarySymbol",
    },
    "proteins": {
        "subject_id": "primaryID",
        "object_id": "secondaryID",
        "primary_label": "name",
        "secondary_label": "synonym",
    },
}


def get_query(entity_type: str) -> str:
    """Get the SPARQL query for a specific entity type.

    Args:
        entity_type: Type of entities (metabolites, chemicals, genes, proteins).

    Returns:
        The SPARQL query string.

    Raises:
        ValueError: If entity type is unknown.
    """
    query = WIKIDATA_QUERIES.get(entity_type.lower())
    if query is None:
        available = list(WIKIDATA_QUERIES.keys())
        raise ValueError(f"Unknown entity type: {entity_type}. Available: {available}")
    return query


def get_column_mapping(entity_type: str) -> dict[str, str]:
    """Get the column name mapping for a specific entity type.

    Args:
        entity_type: Type of entities (metabolites, chemicals, genes, proteins).

    Returns:
        Dictionary mapping semantic names to actual column names.

    Raises:
        ValueError: If entity type is unknown.
    """
    mapping = QUERY_COLUMNS.get(entity_type.lower())
    if mapping is None:
        available = list(QUERY_COLUMNS.keys())
        raise ValueError(f"Unknown entity type: {entity_type}. Available: {available}")
    return mapping
