"""Wikidata parser for redirect mappings via SPARQL queries."""

from __future__ import annotations

import io
from datetime import date
from pathlib import Path
from typing import TYPE_CHECKING
import regex as re

import httpx
import polars as pl

from pysec2pri.constants import WIKIDATA
from pysec2pri.models import IdMapping, MappingCardinality, MappingSet
from pysec2pri.parsers.base import BaseParser
from pysec2pri.queries import get_column_mapping, get_query

if TYPE_CHECKING:
    pass

__all__ = ["WikidataParser", "parse_wikidata", "query_wikidata"]


# QLever endpoint for Wikidata (much faster than official endpoint)
QLEVER_ENDPOINT = "https://qlever.dev/api/wikidata"


def query_wikidata(
    query: str,
    endpoint: str | None = None,
    timeout: float = 100000.0,
) -> pl.DataFrame:
    """Execute a SPARQL query against Wikidata/QLever endpoint.

    The QLever endpoint is used by default as it's much faster than the
    official Wikidata SPARQL endpoint for large queries.

    Args:
        query: The SPARQL query to execute.
        endpoint: SPARQL endpoint URL. Defaults to QLever Wikidata.
        timeout: Request timeout in seconds.

    Returns:
        Polars DataFrame with query results.
    """
    if endpoint is None:
        endpoint = QLEVER_ENDPOINT

    headers = {"Accept": "text/tab-separated-values"}

    with httpx.Client(timeout=timeout) as client:
        response = client.get(
            endpoint,
            params={"query": query},
            headers=headers,
        )
        response.raise_for_status()

    # Parse TSV response using Polars
    tsv_content = response.text

    if not tsv_content.strip():
        return pl.DataFrame()

    # Clean Literals
    cleaned = re.sub(
        r'"([^"\n]*)"@[a-zA-Z-]+',
        r'"\1"',
        tsv_content,
    )

    # Read TSV with Polars
    df = pl.read_csv(
        io.StringIO(cleaned),
        separator="\t",
        has_header=True,
        truncate_ragged_lines=True,
    )

    # Clean up Wikidata URIs in all string columns
    # Format: <http://www.wikidata.org/entity/Q12345> -> Q12345
    for col in df.columns:
        if df[col].dtype == pl.Utf8:
            df = df.with_columns(
                pl.col(col)
                .str.replace(
                    r"<http://www\.wikidata\.org/entity/([^>]+)>",
                    "$1",
                )
                .str.replace(r"@en$", "")  # Remove language tags
                .str.strip_chars()
                .alias(col)
            )

    # Clean column names (remove ? prefix if present)
    new_names = {col: col.lstrip("?") for col in df.columns}
    df = df.rename(new_names)

    return df


class WikidataParser(BaseParser):
    """Parser for Wikidata redirect mappings via SPARQL."""

    datasource_name = "Wikidata"
    default_source_url = "https://www.wikidata.org/"

    def __init__(
        self,
        version: str | None = None,
        show_progress: bool = True,
        entity_type: str = "metabolites",
        endpoint: str | None = None,
    ):
        """Initialize the Wikidata parser.

        Args:
            version: Version/date string for the mappings.
            show_progress: Whether to show progress.
            entity_type: Type of entities to query (metabolites/genes/proteins).
            endpoint: Optional custom SPARQL endpoint.
        """
        super().__init__(version=version, show_progress=show_progress)
        self.entity_type = entity_type
        self.endpoint = endpoint

    def parse(self, input_path: Path | str | None = None) -> MappingSet:
        """Query Wikidata and return a MappingSet.

        Note: input_path is ignored for Wikidata (uses SPARQL queries).

        Args:
            input_path: Ignored for Wikidata.

        Returns:
            A MappingSet containing Wikidata redirect mappings.
        """
        # Get the appropriate query for the entity type
        query = get_query(self.entity_type)
        col_map = get_column_mapping(self.entity_type)

        df = query_wikidata(query, endpoint=self.endpoint)

        mappings: list[IdMapping] = []

        if df.is_empty():
            version = self.version or date.today().isoformat()
            return MappingSet(
                datasource_name=self.datasource_name,
                version=version,
                mappings=mappings,
                curie_map={WIKIDATA.prefix: WIKIDATA.curie_base_url},
            )

        # Get actual column names from the mapping
        primary_id_col = col_map["primary_id"]
        secondary_id_col = col_map["secondary_id"]
        primary_label_col = col_map.get("primary_label")
        secondary_label_col = col_map.get("secondary_label")

        # Process results using Polars DataFrame
        rows = df.iter_rows(named=True)
        if self.show_progress:
            from tqdm import tqdm

            rows = tqdm(
                rows, total=len(df), desc="Processing Wikidata results"
            )

        for row in rows:
            primary_id = row.get(primary_id_col, "")
            secondary_id = row.get(secondary_id_col, "")
            primary_label = row.get(primary_label_col) if primary_label_col else None
            secondary_label = (
                row.get(secondary_label_col) if secondary_label_col else None
            )

            # Create mapping for redirect (secondary -> primary ID mapping)
            if secondary_id and primary_id and secondary_id != primary_id:
                mapping = IdMapping(
                    primary_id=primary_id,
                    secondary_id=secondary_id,
                    primary_label=primary_label,
                    secondary_label=secondary_label,
                    mapping_cardinality=MappingCardinality.ONE_TO_ONE,
                    source_url=self.default_source_url,
                )
                mappings.append(mapping)

            # Create mapping for synonyms (altLabel -> label mapping)
            elif secondary_label and primary_label and secondary_label != primary_label:
                mapping = IdMapping(
                    primary_id=primary_id,
                    secondary_id=None,
                    primary_label=primary_label,
                    secondary_label=secondary_label,
                    source_url=self.default_source_url,
                )
                mappings.append(mapping)

        version = self.version or date.today().isoformat()

        return MappingSet(
            datasource_name=self.datasource_name,
            version=version,
            mappings=mappings,
            curie_map={WIKIDATA.prefix: WIKIDATA.curie_base_url},
        )


def parse_wikidata(
    entity_type: str = "metabolites",
    version: str | None = None,
    endpoint: str | None = None,
    show_progress: bool = True,
) -> MappingSet:
    """Parse Wikidata redirects for a specific entity type.

    Args:
        entity_type: Type of entities (metabolites, genes, proteins).
        version: Version string for the mappings.
        endpoint: Optional custom SPARQL endpoint.
        show_progress: Whether to show progress bars.

    Returns:
        MappingSet containing Wikidata redirect mappings.

    Example:
        >>> from pysec2pri.parsers.wikidata import parse_wikidata
        >>> mappings = parse_wikidata("metabolites")
        >>> print(f"Found {len(mappings)} redirects")
    """
    parser = WikidataParser(
        version=version,
        show_progress=show_progress,
        entity_type=entity_type,
        endpoint=endpoint,
    )
    return parser.parse(None)
