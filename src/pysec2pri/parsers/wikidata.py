"""Wikidata parser for redirect mappings via SPARQL queries."""

from __future__ import annotations

import io
from datetime import date
from pathlib import Path

import httpx
import polars as pl
import regex as re

from pysec2pri.constants import WIKIDATA
from pysec2pri.models import (
    IdMapping,
    MappingCardinality,
    MappingSet,
    SymbolMapping,
)
from pysec2pri.parsers.base import BaseParser
from pysec2pri.queries import WIKIDATA_QUERIES, WIKIDATA_TEST_QUERIES, get_column_mapping

# QLever endpoint for Wikidata (much faster than official endpoint)
QLEVER_ENDPOINT = "https://qlever.dev/api/wikidata"


def query_wikidata(
    query: str,
    endpoint: str | None = None,
    timeout: float = 100000.0,
) -> pl.DataFrame:
    """Execute a SPARQL query against Wikidata/QLever endpoint.

    The QLever endpoint is used by default.

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
        test_subset: bool = False,
    ):
        """Initialize the Wikidata parser.

        Args:
            version: Version/date string for the mappings.
            show_progress: Whether to show progress.
            entity_type: Type of entities to query (metabolites/genes/proteins).
            endpoint: Optional custom SPARQL endpoint.
            test_subset: Whether to use test queries (LIMIT 10)
        """
        super().__init__(version=version, show_progress=show_progress)
        self.entity_type = entity_type
        self.endpoint = endpoint
        self.test_subset = test_subset

    def parse(
        self,
        input_path: Path | str | None = None,
        entity_type: str = "metabolites",
        version: str | None = None,
        endpoint: str | None = None,
        show_progress: bool = True,
    ) -> MappingSet:
        """Query Wikidata and return a MappingSet.

        If test_subset is True, use the test query for the entity type.

        Args:
            entity_type: Type of entities (metabolites, genes, proteins).
            input_path: None, from endpoint.
            entity_type: Query to be selected.
            version: Version string for the mappings.
            endpoint: Optional custom SPARQL endpoint.
            show_progress: Whether to show progress bars.
            test_subset: Whether to use the test queries (LIMIT 10).

        Returns:
            MappingSet containing Wikidata redirect mappings.
        """
        if self.test_subset:
            query_str = WIKIDATA_TEST_QUERIES.get(self.entity_type)
            if query_str is None:
                available = list(WIKIDATA_TEST_QUERIES.keys())
                raise ValueError(
                    f"Unknown entity type for test subset: {self.entity_type}. "
                    f"Available: {available}"
                )
        else:
            query_str = WIKIDATA_QUERIES.get(self.entity_type)
            if query_str is None:
                available = list(WIKIDATA_QUERIES.keys())
                raise ValueError(f"Unknown entity type: {self.entity_type}. Available: {available}")
        df = query_wikidata(query_str, endpoint=self.endpoint)

        if df.is_empty():
            return self._empty_mappingset()

        df = self._normalize_ids(df)
        mappings = self._build_mappings(df)
        version = self.version or date.today().isoformat()

        return MappingSet(
            datasource_name=self.datasource_name,
            version=version,
            mappings=mappings,
            curie_map={WIKIDATA.prefix: WIKIDATA.curie_base_url},
            comment=self._build_comment(f"Wikidata mappings for entity type: {self.entity_type}"),
        )

    def _empty_mappingset(self) -> MappingSet:
        version = self.version or date.today().isoformat()
        return MappingSet(
            datasource_name=self.datasource_name,
            version=version,
            mappings=[],
            curie_map={WIKIDATA.prefix: WIKIDATA.curie_base_url},
            comment=self._build_comment("No mappings found for this entity type."),
        )

    def _normalize_ids(self, df: pl.DataFrame) -> pl.DataFrame:
        col_map = get_column_mapping(self.entity_type)
        return df.with_columns(
            [
                pl.col(col_map["subject_id"])
                .map_elements(self.normalize_subject_id)
                .alias("subject_id_norm"),
                pl.col(col_map["object_id"])
                .map_elements(self.normalize_subject_id)
                .alias("object_id_norm"),
            ]
        )

    def _build_mappings(self, df: pl.DataFrame) -> list[IdMapping | SymbolMapping]:
        col_map = get_column_mapping(self.entity_type)
        primary_label_col = col_map.get("primary_label")
        secondary_label_col = col_map.get("secondary_label")

        # Redirect mappings - filter for valid id pairs
        redirect_mask = (
            (pl.col("object_id_norm").is_not_null())
            & (pl.col("subject_id_norm").is_not_null())
            & (pl.col("object_id_norm") != pl.col("subject_id_norm"))
        )
        redirects_df = df.filter(redirect_mask)
        id_pairs = list(
            zip(
                redirects_df["object_id_norm"].to_list(),
                redirects_df["subject_id_norm"].to_list(),
                strict=False,
            )
        )
        cardinality_map = self.compute_cardinality(id_pairs)

        mappings: list[IdMapping | SymbolMapping] = []
        for row in redirects_df.iter_rows(named=True):
            # Get labels if columns exist
            subj_label = row.get(primary_label_col) if primary_label_col else None
            obj_label = row.get(secondary_label_col) if secondary_label_col else None
            mappings.append(
                IdMapping(
                    subject_id=row["subject_id_norm"],
                    object_id=row["object_id_norm"],
                    subject_label=subj_label,
                    object_label=obj_label,
                    mapping_cardinality=cardinality_map.get(
                        (row["object_id_norm"], row["subject_id_norm"]),
                        MappingCardinality.ONE_TO_ONE,
                    ),
                    predicate_id="oboInOwl:consider",
                    source_url=self.default_source_url,
                    comment=self._build_comment("Wikidata redirect mapping."),
                )
            )

        # Synonym mappings (symbols without ID redirects)
        if primary_label_col and secondary_label_col:
            synonym_mask = (
                (pl.col("subject_id_norm").is_not_null())
                & (pl.col(primary_label_col).is_not_null())
                & (pl.col(secondary_label_col).is_not_null())
                & (pl.col(primary_label_col) != pl.col(secondary_label_col))
            )
            synonyms_df = df.filter(synonym_mask)
            for row in synonyms_df.iter_rows(named=True):
                mappings.append(
                    SymbolMapping(
                        subject_id=row["subject_id_norm"],
                        object_id=row.get("object_id_norm"),
                        subject_label=row[primary_label_col],
                        object_label=row[secondary_label_col],
                        predicate_id="oboInOwl:consider",
                        source_url=self.default_source_url,
                        comment=self._build_comment("Wikidata symbol/synonym mapping."),
                    )
                )

        return mappings


def parse_wikidata(
    entity_type: str = "metabolites",
    version: str | None = None,
    endpoint: str | None = None,
    show_progress: bool = True,
    test_subset: bool = False,
) -> MappingSet:
    """Parse Wikidata redirects for a specific entity type.

    Args:
        entity_type: Type of entities (metabolites, genes, proteins).
        version: Version string for the mappings.
        endpoint: Optional custom SPARQL endpoint.
        show_progress: Whether to show progress bars.
        test_subset: whether to use test queries (LIMIT 10)

    Returns:
        MappingSet containing Wikidata redirect mappings.

    Example:
        >>> from pysec2pri.parsers.wikidata import parse_wikidata
        >>> mappings = parse_wikidata("metabolites")
        >>> print(f"Found {len(mappings)} redirects")
    """
    return WikidataParser(
        version=version,
        show_progress=show_progress,
        entity_type=entity_type,
        endpoint=endpoint,
        test_subset=test_subset,
    ).parse()
