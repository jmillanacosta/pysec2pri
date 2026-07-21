"""Wikidata parser for redirect mappings via SPARQL queries.

This parser extracts ID-to-ID mappings for chemicals, genes, and proteins
from Wikidata via SPARQL queries.

Uses SSSOM-compliant IdMappingSet with cardinality computation.
"""

from __future__ import annotations

import io
from pathlib import Path

import httpx
import polars as pl
from sssom_schema import Mapping

from pysec2pri.logging import logger
from pysec2pri.parsers.base import (
    BaseMappingSet,
    BaseParser,
    LabelMappingSet,
    get_datasource_config,
)
from pysec2pri.queries import (
    WIKIDATA_QUERIES,
    WIKIDATA_TEST_QUERIES,
    get_column_mapping,
)

# Default QLever endpoint (fallback if not in config)
DEFAULT_QLEVER_ENDPOINT = "https://qlever.dev/api/wikidata"

__all__ = ["WikidataParser", "query_wikidata"]


def query_wikidata(
    query: str,
    endpoint: str | None = None,
    timeout: float = 6000.0,
) -> pl.DataFrame:
    """Execute a SPARQL query against Wikidata/QLever endpoint.

    Args:
        query: The SPARQL query to execute.
        endpoint: SPARQL endpoint URL. Defaults to QLever Wikidata.
        timeout: Request timeout in seconds.

    Returns:
        Polars DataFrame with query results.
    """
    if endpoint is None:
        endpoint = DEFAULT_QLEVER_ENDPOINT

    headers = {"Accept": "text/tab-separated-values"}

    with httpx.Client(timeout=timeout) as client:
        response = client.get(
            endpoint,
            params={"query": query},
            headers=headers,
        )
        try:
            response.raise_for_status()
        except httpx.HTTPStatusError as e:
            logger.warning("Skipping query - %s", e)
            return pl.DataFrame()

    # Parse TSV response using Polars
    tsv_content = response.text

    if not tsv_content.strip():
        return pl.DataFrame()

    uri_pat = r"<http://www\.wikidata\.org/entity/([^>]+)>"
    df = (
        pl.scan_csv(
            io.StringIO(tsv_content),
            separator="\t",
            has_header=True,
            truncate_ragged_lines=True,
            infer_schema_length=0,
            quote_char=None,
        )
        .with_columns(
            pl.col(pl.Utf8)
            .str.replace(uri_pat, "$1")
            .str.replace(r'^"(.*)"@[a-zA-Z]+(?:-[a-zA-Z0-9]+)*$', r"$1")
            .str.strip_chars()
        )
        .rename(lambda c: c.lstrip("?"))
        .collect()
    )

    return df


class WikidataParser(BaseParser):
    """Parser for Wikidata redirect mappings via SPARQL.

    Queries the QLever Wikidata endpoint to find redirect mappings
    for chemicals, genes, and proteins.

    Returns IdMappingSet for all mappings.
    """

    datasource_name = "wikidata"
    default_source_url = "https://www.wikidata.org/"

    @classmethod
    def entity_types(cls) -> list[str]:
        """Return the entity types declared by ``wikidata.yaml``'s ``queries`` block.

        Each key names both a ``--entity-type`` choice and its redirect query,
        so the config is the only place they are listed.
        """
        cfg = get_datasource_config(cls.datasource_name, config_package=cls.config_package)
        return list(cfg.queries)

    def __init__(
        self,
        version: str | None = None,
        show_progress: bool = True,
        entity_type: str | None = None,
        endpoint: str | None = None,
        test_subset: bool = False,
    ):
        """Initialize the Wikidata parser.

        Args:
            version: Version/date string for the mappings.
            show_progress: Whether to show progress.
            entity_type: Entity type to query; one of :meth:`entity_types`.
                Defaults to the first one the config declares.
            endpoint: Optional custom SPARQL endpoint.
            test_subset: Whether to use test queries (LIMIT 10).
        """
        super().__init__(version=version, show_progress=show_progress)
        self.entity_type = entity_type if entity_type is not None else self.entity_types()[0]
        self.test_subset = test_subset

        # Use provided endpoint, or fall back to config, or default
        if endpoint:
            self.endpoint = endpoint
        elif self._config and self._config.sparql_endpoint:
            self.endpoint = self._config.sparql_endpoint
        else:
            self.endpoint = DEFAULT_QLEVER_ENDPOINT

    def parse(self, input_path: Path | str | None = None) -> BaseMappingSet:
        """Query Wikidata and return a MappingSet.

        Args:
            input_path: Ignored for Wikidata (queries endpoint directly).

        Returns:
            IdMappingSet containing Wikidata redirect mappings.
        """
        # Select query based on test_subset flag
        if self.test_subset:
            query_str = WIKIDATA_TEST_QUERIES.get(self.entity_type)
            if query_str is None:
                available = list(WIKIDATA_TEST_QUERIES.keys())
                raise ValueError(f"Unknown entity type: {self.entity_type}. Available: {available}")
        else:
            query_str = WIKIDATA_QUERIES.get(self.entity_type)
            if query_str is None:
                available = list(WIKIDATA_QUERIES.keys())
                raise ValueError(f"Unknown entity type: {self.entity_type}. Available: {available}")

        logger.info(
            "Querying Wikidata for %s redirects%s...",
            self.entity_type,
            " (test subset)" if self.test_subset else "",
        )

        df = query_wikidata(query_str, endpoint=self.endpoint)

        if df.is_empty():
            self._resolve_version()
            return self.create_mapping_set([], mapping_type="id")

        df = self._normalize_ids(df)
        mappings = self._build_redirect_mappings(df)
        self._resolve_version()

        return self.create_mapping_set(mappings, mapping_type="id")

    def parse_all(self) -> BaseMappingSet:
        """Query all entity types from config and return combined MappingSet.

        Runs all SPARQL queries defined in the config file's 'queries' section
        (e.g., chemical_redirects, gene_redirects, protein_redirects) and
        combines the results into a single MappingSet.

        Returns:
            IdMappingSet containing all Wikidata redirect mappings.
        """
        all_mappings: list[Mapping] = []

        for entity_type in self.entity_types():
            # Select appropriate query
            if self.test_subset:
                query_str = WIKIDATA_TEST_QUERIES.get(entity_type)
            else:
                query_str = WIKIDATA_QUERIES.get(entity_type)

            if not query_str:
                logger.warning("No query found for %s, skipping", entity_type)
                continue

            logger.info(
                "Querying Wikidata for %s%s...",
                entity_type,
                " (test subset)" if self.test_subset else "",
            )

            original_entity = self.entity_type
            self.entity_type = entity_type
            try:
                df = query_wikidata(query_str, endpoint=self.endpoint)

                if df.is_empty():
                    logger.info("No results for %s", entity_type)
                    continue

                df = self._normalize_ids(df)
                mappings = self._build_redirect_mappings(df)
                all_mappings.extend(mappings)
            except Exception as e:
                logger.warning("Failed to query %s: %s", entity_type, e)
                continue
            finally:
                self.entity_type = original_entity

        self._resolve_version()
        return self.create_mapping_set(all_mappings, mapping_type="id")

    def parse_from_file(self, input_path: Path | str) -> BaseMappingSet:
        """Parse Wikidata redirects from a pre-downloaded TSV file.

        Args:
            input_path: Path to TSV file with SPARQL results.

        Returns:
            IdMappingSet with computed cardinalities.
        """
        input_path = Path(input_path)

        df = pl.read_csv(input_path, separator="\t", has_header=True)
        df = self._normalize_ids(df)
        mappings = self._build_redirect_mappings(df)
        self._resolve_version()

        return self.create_mapping_set(mappings, mapping_type="id")

    def parse_labels(self, input_path: Path | str | None = None) -> LabelMappingSet:
        """Return a LabelMappingSet of previous-label  to current-label mappings.

        Queries the SPARQL endpoint (or reads *input_path*) exactly like
        :meth:`parse`, but wraps the result in a :class:`LabelMappingSet`
        so label-specific exports (``label_sec2pri``, ``pri_labels``) work.

        Args:
            input_path: Pre-downloaded TSV file. Queries SPARQL if ``None``.

        Returns:
            :class:`LabelMappingSet` with label-based mappings.
        """
        if input_path is not None:
            input_path = Path(input_path)
            df = pl.read_csv(input_path, separator="\t", has_header=True)
        else:
            if self.test_subset:
                query_str = WIKIDATA_TEST_QUERIES.get(self.entity_type)
            else:
                query_str = WIKIDATA_QUERIES.get(self.entity_type)
            if query_str is None:
                available = list(WIKIDATA_QUERIES.keys())
                raise ValueError(f"Unknown entity type: {self.entity_type}. Available: {available}")
            df = query_wikidata(query_str, endpoint=self.endpoint)

        if df.is_empty():
            self._resolve_version()
            return self.create_mapping_set([], mapping_type="label")

        df = self._normalize_ids(df)
        mappings = self._build_redirect_mappings(df)
        self._resolve_version()
        return self.create_mapping_set(mappings, mapping_type="label")

    def _normalize_ids(self, df: pl.DataFrame) -> pl.DataFrame:
        """Normalize IDs to WD:Qxxx format."""
        col_map = get_column_mapping(self.entity_type)
        subject_col = col_map["subject_id"]
        object_col = col_map["object_id"]

        if subject_col not in df.columns:
            logger.warning("Column %s not found in results", subject_col)
            return df

        result = df.with_columns(
            self._normalize_qid_expr(subject_col).alias("subject_id_norm"),
        )

        if object_col in df.columns:
            result = result.with_columns(
                self._normalize_qid_expr(object_col).alias("object_id_norm"),
            )
        else:
            result = result.with_columns(
                pl.lit(None, dtype=pl.Utf8).alias("object_id_norm"),
            )

        return result

    @staticmethod
    def _normalize_qid_expr(col_name: str) -> pl.Expr:
        """Normalize QID column."""
        c = pl.col(col_name).str.strip_chars()
        return (
            pl.when(pl.col(col_name).is_null() | (c == ""))
            .then(pl.lit(None, dtype=pl.Utf8))
            .when(c.str.starts_with("WD:"))
            .then(c)
            .when(c.str.contains(r"wikidata\.org/entity/"))
            .then(pl.lit("WD:") + c.str.split("/").list.last())
            .otherwise(pl.lit("WD:") + c)
        )

    def _build_redirect_mappings(self, df: pl.DataFrame) -> list[Mapping]:
        col_map = get_column_mapping(self.entity_type)
        primary_label_col = col_map.get("primary_label")
        secondary_label_col = col_map.get("secondary_label")
        m_meta = self.get_mapping_metadata()

        fixed = {
            "predicate_id": m_meta.get("predicate_id", "oboInOwl:consider"),
            "predicate_label": m_meta.get("predicate_label"),
            "mapping_justification": m_meta.get(
                "mapping_justification",
                "semapv:BackgroundKnowledgeBasedMatching",
            ),
            "subject_source": self.default_source_url,
            "object_source": self.default_source_url,
            "mapping_tool": m_meta.get("mapping_tool"),
            "license": m_meta.get("license"),
            "comment": self._build_comment(f"Wikidata {self.entity_type} redirect."),
        }

        # Build the minimal Polars projection
        cols: list[str] = ["subject_id_norm", "object_id_norm"]
        if primary_label_col and primary_label_col in df.columns:
            cols.append(primary_label_col)
        if secondary_label_col and secondary_label_col in df.columns:
            cols.append(secondary_label_col)

        redirects_df = (
            df.lazy()
            .filter(
                pl.col("object_id_norm").is_not_null()
                & pl.col("subject_id_norm").is_not_null()
                & (pl.col("object_id_norm") != pl.col("subject_id_norm"))
            )
            .select(cols)
            .collect()
        )

        col_idx: dict[str, int] = {c: i for i, c in enumerate(redirects_df.columns)}
        subj_i = col_idx["subject_id_norm"]
        obj_i = col_idx["object_id_norm"]
        subj_label_i = col_idx.get(secondary_label_col) if secondary_label_col else None
        obj_label_i = col_idx.get(primary_label_col) if primary_label_col else None

        rows_data = [
            {
                "subject_id": row[subj_i],
                "object_id": row[obj_i],
                "subject_label": row[subj_label_i] if subj_label_i is not None else None,
                "object_label": row[obj_label_i] if obj_label_i is not None else None,
            }
            for row in redirects_df.iter_rows()
        ]

        return self._build_mappings(
            rows_data,
            fixed,
            desc=f"Building {self.entity_type} mappings",
            total=len(rows_data),
        )
