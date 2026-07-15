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
# Query mapping by entity type, keyed as in wikidata.yaml's `queries` block.
WIKIDATA_QUERIES: dict[str, str] = {
    "chemicals": CHEMICAL_REDIRECTS_QUERY,
    "genes": GENE_REDIRECTS_QUERY,
    "proteins": PROTEIN_REDIRECTS_QUERY,
}

WIKIDATA_TEST_QUERIES: dict[str, str] = {
    "chemicals": CHEMICAL_REDIRECTS_TEST_QUERY,
    "genes": GENE_REDIRECTS_TEST_QUERY,
    "proteins": PROTEIN_REDIRECTS_TEST_QUERY,
}

# Result column names per entity type. Only genes name their label columns
# after symbols; the rest share the name/synonym pair.
_NAME_SYNONYM_COLUMNS = {
    "subject_id": "secondaryID",
    "object_id": "primaryID",
    "primary_label": "name",
    "secondary_label": "synonym",
}
QUERY_COLUMNS: dict[str, dict[str, str]] = {
    "chemicals": _NAME_SYNONYM_COLUMNS,
    "genes": {
        "subject_id": "secondaryID",
        "object_id": "primaryID",
        "primary_label": "primarySymbol",
        "secondary_label": "secondarySymbol",
    },
    "proteins": _NAME_SYNONYM_COLUMNS,
}


def get_column_mapping(entity_type: str) -> dict[str, str]:
    """Get the column name mapping for a specific entity type.

    Args:
        entity_type: Entity type, e.g. ``"chemicals"``/``"genes"``/``"proteins"``.

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
