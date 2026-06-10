"""Wikidata parser for redirect mappings via SPARQL queries.

This parser extracts ID-to-ID mappings for chemicals, genes, and proteins
from Wikidata via SPARQL queries.

Uses SSSOM-compliant IdMappingSet with cardinality computation.
"""

from __future__ import annotations

import io
from datetime import date
from pathlib import Path

import httpx
import polars as pl
from sssom_schema import Mapping

from pysec2pri.logging import logger
from pysec2pri.parsers.base import (
    BaseParser,
    IdMappingSet,
    LabelMappingSet,
    Sec2PriMappingSet,
)
from pysec2pri.queries import (
    WIKIDATA_QUERIES,
    WIKIDATA_TEST_QUERIES,
    get_column_mapping,
)
from pysec2pri.version import VERSION

# Default QLever endpoint (fallback if not in config)
DEFAULT_QLEVER_ENDPOINT = "https://qlever.dev/api/wikidata"

__all__ = ["WikidataParser", "query_wikidata"]


def query_wikidata(
    query: str,
    endpoint: str | None = None,
    timeout: float = 100000.0,
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
            entity_type: Type of entities to query.
            endpoint: Optional custom SPARQL endpoint.
            test_subset: Whether to use test queries (LIMIT 10).
        """
        super().__init__(version=version, show_progress=show_progress)
        self.entity_type = entity_type
        self.test_subset = test_subset

        # Use provided endpoint, or fall back to config, or default
        if endpoint:
            self.endpoint = endpoint
        elif self._config and self._config.sparql_endpoint:
            self.endpoint = self._config.sparql_endpoint
        else:
            self.endpoint = DEFAULT_QLEVER_ENDPOINT

    def parse(self, input_path: Path | str | None = None) -> Sec2PriMappingSet:
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
            return self._empty_mappingset()

        df = self._normalize_ids(df)
        mappings = self._build_redirect_mappings(df)
        version = self._resolve_version()

        return self._create_mapping_set(mappings, version)

    def parse_all(self) -> Sec2PriMappingSet:
        """Query all entity types from config and return combined MappingSet.

        Runs all SPARQL queries defined in the config file's 'queries' section
        (e.g., chemical_redirects, gene_redirects, protein_redirects) and
        combines the results into a single MappingSet.

        Returns:
            IdMappingSet containing all Wikidata redirect mappings.
        """
        all_mappings: list[Mapping] = []

        # Get query types from config or use defaults
        if self._config and self._config.queries:
            query_types = list(self._config.queries.keys())
        else:
            query_types = [
                "chemical_redirects",
                "gene_redirects",
                "protein_redirects",
            ]

        # Map config query names to entity types
        query_to_entity = {
            "chemical_redirects": "metabolites",
            "gene_redirects": "genes",
            "protein_redirects": "proteins",
        }

        for query_name in query_types:
            entity_type = query_to_entity.get(query_name)
            if not entity_type:
                logger.warning("Unknown query type: %s, skipping", query_name)
                continue

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

            try:
                df = query_wikidata(query_str, endpoint=self.endpoint)

                if df.is_empty():
                    logger.info("No results for %s", entity_type)
                    continue

                # Temporarily set entity_type for normalization
                original_entity = self.entity_type
                self.entity_type = entity_type

                df = self._normalize_ids(df)
                mappings = self._build_redirect_mappings(df)
                all_mappings.extend(mappings)

                self.entity_type = original_entity

            except Exception as e:
                logger.warning("Failed to query %s: %s", entity_type, e)
                continue

        version = self._resolve_version()
        return self._create_mapping_set(all_mappings, version)

    def parse_from_file(self, input_path: Path | str) -> Sec2PriMappingSet:
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
        version = self._resolve_version()

        return self._create_mapping_set(mappings, version)

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
            return self._create_label_mapping_set([], self._resolve_version())

        df = self._normalize_ids(df)
        mappings = self._build_redirect_mappings(df)
        version = self._resolve_version()
        return self._create_label_mapping_set(mappings, version)

    def _create_label_mapping_set(
        self,
        mappings: list[Mapping],
        version: str,
    ) -> LabelMappingSet:
        """Create a LabelMappingSet with metadata."""
        ms_meta = self.get_mappingset_metadata()

        curie_map = ms_meta.get("curie_map", {})
        curie_map_str: dict[str, str] = {}
        for k, v in curie_map.items():
            if hasattr(v, "prefix_reference"):
                curie_map_str[str(k)] = str(v.prefix_reference)
            else:
                curie_map_str[str(k)] = str(v)

        mapping_set = LabelMappingSet(
            mapping_set_id=ms_meta.get("mapping_set_id", ""),
            mapping_set_version=version,
            mapping_set_title=ms_meta.get("mapping_set_title"),
            mapping_set_description=ms_meta.get("mapping_set_description"),
            curie_map=curie_map_str,
            license=ms_meta.get("license"),
            creator_id=ms_meta.get("creator_id"),
            creator_label=ms_meta.get("creator_label"),
            mapping_provider=ms_meta.get("mapping_provider"),
            mapping_tool=ms_meta.get("mapping_tool"),
            mapping_tool_version=VERSION,
            mapping_date=date.today().isoformat(),
            subject_source=ms_meta.get("subject_source"),
            object_source=ms_meta.get("object_source"),
            subject_source_version=version,
            object_source_version=version,
            comment=self._build_comment(f"Wikidata label mappings for {self.entity_type}."),
            mappings=mappings,
        )
        mapping_set.compute_cardinalities()
        return mapping_set

    def _empty_mappingset(self) -> Sec2PriMappingSet:
        """Create an empty mapping set."""
        version = self._resolve_version()
        ms_meta = self.get_mappingset_metadata()

        curie_map = ms_meta.get("curie_map", {})
        curie_map_str: dict[str, str] = {}
        for k, v in curie_map.items():
            if hasattr(v, "prefix_reference"):
                curie_map_str[str(k)] = str(v.prefix_reference)
            else:
                curie_map_str[str(k)] = str(v)

        mapping_set = IdMappingSet(
            mapping_set_id=ms_meta.get("mapping_set_id", ""),
            mapping_set_version=version,
            curie_map=curie_map_str,
            license=ms_meta.get("license"),
            mapping_tool=ms_meta.get("mapping_tool"),
            mapping_tool_version=VERSION,
            mapping_date=date.today().isoformat(),
            subject_source=ms_meta.get("subject_source"),
            object_source=ms_meta.get("object_source"),
            subject_source_version=version,
            object_source_version=version,
            mappings=[],
        )
        mapping_set.compute_cardinalities()
        return mapping_set

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

    @staticmethod
    def _normalize_qid(qid: str | None) -> str | None:
        """Normalize a QID to WD:Qxxx format."""
        if not qid or qid == "":
            return None
        qid = qid.strip()
        if qid.startswith("WD:"):
            return qid
        if "wikidata.org/entity/" in qid:
            qid = qid.split("/")[-1]
        return f"WD:{qid}"

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

    def _create_mapping_set(
        self,
        mappings: list[Mapping],
        version: str,
    ) -> Sec2PriMappingSet:
        """Create the final MappingSet with metadata."""
        ms_meta = self.get_mappingset_metadata()

        curie_map = ms_meta.get("curie_map", {})
        curie_map_str: dict[str, str] = {}
        for k, v in curie_map.items():
            if hasattr(v, "prefix_reference"):
                curie_map_str[str(k)] = str(v.prefix_reference)
            else:
                curie_map_str[str(k)] = str(v)

        mapping_set = IdMappingSet(
            mapping_set_id=ms_meta.get("mapping_set_id", ""),
            mapping_set_version=version,
            mapping_set_title=ms_meta.get("mapping_set_title"),
            mapping_set_description=ms_meta.get("mapping_set_description"),
            curie_map=curie_map_str,
            license=ms_meta.get("license"),
            creator_id=ms_meta.get("creator_id"),
            creator_label=ms_meta.get("creator_label"),
            mapping_provider=ms_meta.get("mapping_provider"),
            mapping_tool=ms_meta.get("mapping_tool"),
            mapping_tool_version=VERSION,
            mapping_date=date.today().isoformat(),
            subject_source=ms_meta.get("subject_source"),
            object_source=ms_meta.get("object_source"),
            subject_source_version=version,
            object_source_version=version,
            comment=self._build_comment(f"Wikidata mappings for {self.entity_type}."),
            mappings=mappings,
        )
        mapping_set.compute_cardinalities()
        return mapping_set
