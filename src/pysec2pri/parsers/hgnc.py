"""HGNC TSV file parser for secondary-to-primary identifier mappings.

This parser extracts:
1. ID-to-ID mappings: withdrawn/merged HGNC IDs -> current HGNC IDs
2. Label-to-label mappings: previous/alias symbols -> current symbols

Uses SSSOM-compliant MappingSet classes with cardinality computation.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
from sssom_schema import Mapping

from pysec2pri.parsers.base import (
    WITHDRAWN_ENTRY,
    WITHDRAWN_ENTRY_LABEL,
    BaseParser,
    Sec2PriMappingSet,
)

# HGNC column names (case-insensitive matching used)
HGNC_ID = "hgnc_id"
SYMBOL = "symbol"
ALIAS_SYMBOL = "alias_symbol"
PREV_SYMBOL = "prev_symbol"
STATUS = "status"

# Merged info column has different naming variants across HGNC file versions
MERGED_INFO_PATTERNS = [
    "merged_into_report(i.e. hgnc_id/symbol/status)",
    "merged_into_report(i.e hgnc_id/symbol/status)",
    "merged_into_report(s) (i.e hgnc_id|symbol|status)",
]


class HGNCParser(BaseParser):
    """Parser for HGNC TSV files using Polars for memory efficiency.

    Extracts secondary-to-primary HGNC identifier mappings and
    symbol mappings from HGNC withdrawn and complete set files.

    Returns:
    - IdMappingSet for ID-to-ID mappings (withdrawn/merged IDs)
    - LabelMappingSet for symbol mappings (alias/previous symbols)
    """

    datasource_name = "hgnc"

    def __init__(
        self,
        version: str | None = None,
        show_progress: bool = True,
    ):
        """Initialize the HGNC parser.

        Args:
            version: Version/release identifier for the datasource.
            show_progress: Whether to show progress bars during parsing.
        """
        super().__init__(version=version, show_progress=show_progress)

    @property
    def withdrawn_source_url(self) -> str:
        """Get the withdrawn file download URL from config."""
        return self.get_download_url("withdrawn") or ""

    @property
    def complete_set_source_url(self) -> str:
        """Get the complete set download URL from config."""
        return self.get_download_url("complete") or ""

    def parse(self, input_path: Path | str | None) -> Sec2PriMappingSet:
        """Parse HGNC withdrawn TSV file into an IdMappingSet.

        Args:
            input_path: Path to the withdrawn HGNC TSV file.

        Returns:
            IdMappingSet with computed cardinalities based on IDs.
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        input_path = Path(input_path)

        # Parse withdrawn file for ID mappings
        mappings = self._parse_withdrawn(input_path)

        # Create IdMappingSet and compute cardinalities
        mapping_set = self._create_mapping_set(mappings, mapping_type="id")
        return mapping_set

    def parse_symbols(self, complete_set_path: Path | str | None) -> Sec2PriMappingSet:
        """Parse HGNC complete set for symbol (label) mappings.

        Args:
            complete_set_path: Path to the complete HGNC set TSV file.

        Returns:
            LabelMappingSet with computed cardinalities based on labels.
        """
        if complete_set_path is None:
            raise ValueError("complete_set_path must not be None")
        complete_set_path = Path(complete_set_path)

        # Parse complete set for symbol mappings
        mappings = self._parse_complete_set(complete_set_path)

        # Create LabelMappingSet and compute cardinalities
        mapping_set = self._create_mapping_set(mappings, mapping_type="label")
        return mapping_set

    def parse_all(
        self,
        withdrawn_path: Path | str | None,
        complete_set_path: Path | str | None,
    ) -> tuple[Sec2PriMappingSet, Sec2PriMappingSet]:
        """Parse both withdrawn and complete set files.

        Args:
            withdrawn_path: Path to the withdrawn HGNC TSV file.
            complete_set_path: Path to the complete HGNC set TSV file.

        Returns:
            Tuple of (IdMappingSet, LabelMappingSet).
        """
        id_mappings = self.parse(withdrawn_path)
        label_mappings = self.parse_symbols(complete_set_path)
        return id_mappings, label_mappings

    def _parse_withdrawn(self, file_path: Path) -> list[Mapping]:
        """Parse withdrawn HGNC file for ID-to-ID mappings.

        Args:
            file_path: Path to the withdrawn HGNC TSV file.

        Returns:
            List of SSSOM Mapping objects.
        """
        df = pl.read_csv(
            file_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=[""],
        )

        merged_col = self._find_merged_column(df.columns, MERGED_INFO_PATTERNS)
        if merged_col is None:
            raise ValueError(f"Could not find merged_into_report column in {file_path}")

        hgnc_id_col = self._find_column(df.columns, HGNC_ID)
        if hgnc_id_col is None:
            raise ValueError(f"Could not find hgnc_id column in {file_path}")

        status_col = self._find_column(df.columns, STATUS)
        symbol_col = self._find_column(df.columns, SYMBOL)

        m_meta = self.get_mapping_metadata()
        mappings: list[Mapping] = []

        rows = list(df.iter_rows(named=True))
        for row in self._progress(rows, desc="Processing withdrawn"):
            hgnc_id = row.get(hgnc_id_col)
            if not hgnc_id:
                continue

            merged_info = row.get(merged_col)
            status = row.get(status_col) if status_col else None
            symbol = row.get(symbol_col) if symbol_col else None

            # Case 1: Withdrawn with no replacement
            if not merged_info and status and "Entry Withdrawn" in str(status):
                mapping = Mapping(
                    subject_id=WITHDRAWN_ENTRY,
                    object_id=hgnc_id,
                    subject_label=WITHDRAWN_ENTRY_LABEL,
                    object_label=symbol or "",
                    predicate_id="oboInOwl:consider",
                    mapping_justification=m_meta["mapping_justification"],
                    subject_source=m_meta.get("subject_source"),
                    object_source=m_meta.get("object_source"),
                    mapping_tool=m_meta.get("mapping_tool"),
                    license=m_meta.get("license"),
                    comment="Withdrawn entry with no replacement.",
                )
                mappings.append(mapping)
                continue

            # Case 2: Merged into another entry
            if merged_info:
                parsed = self._parse_merged_info(merged_info)
                if parsed:
                    target_id, target_symbol = parsed
                    # SSSOM: subject = primary (target), object = secondary (old)
                    mapping = Mapping(
                        subject_id=target_id,
                        object_id=hgnc_id,
                        subject_label=target_symbol or "",
                        object_label=symbol or "",
                        predicate_id=m_meta["predicate_id"],
                        mapping_justification=m_meta["mapping_justification"],
                        subject_source=m_meta.get("subject_source"),
                        object_source=m_meta.get("object_source"),
                        mapping_tool=m_meta.get("mapping_tool"),
                        license=m_meta.get("license"),
                    )
                    mappings.append(mapping)

        return mappings

    def _parse_complete_set(self, file_path: Path) -> list[Mapping]:
        """Parse complete HGNC set for symbol (label) mappings.

        Args:
            file_path: Path to the complete HGNC set TSV file.

        Returns:
            List of SSSOM Mapping objects for symbol mappings.
        """
        df = pl.read_csv(
            file_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=[""],
        )

        status_col = self._find_column(df.columns, STATUS)
        hgnc_id_col = self._find_column(df.columns, HGNC_ID)
        symbol_col = self._find_column(df.columns, SYMBOL)
        alias_col = self._find_column(df.columns, ALIAS_SYMBOL)
        prev_col = self._find_column(df.columns, PREV_SYMBOL)

        if not all([status_col, hgnc_id_col, symbol_col]):
            raise ValueError(f"Missing required columns in {file_path}")

        # Filter to approved entries only
        df_approved = df.filter(
            pl.col(status_col)  # type: ignore[arg-type]
            == "Approved"
        )

        m_meta = self.get_mapping_metadata()
        mappings: list[Mapping] = []

        rows = list(df_approved.iter_rows(named=True))
        for row in self._progress(rows, desc="Processing symbols"):
            hgnc_id = row.get(hgnc_id_col)
            symbol = row.get(symbol_col)
            if not hgnc_id or not symbol:
                continue

            alias_str = row.get(alias_col) if alias_col else None
            prev_str = row.get(prev_col) if prev_col else None
            aliases = self._split_symbols(alias_str) if alias_str else []
            prev_symbols = self._split_symbols(prev_str) if prev_str else []

            # Create label mappings for aliases
            for alias in aliases:
                mapping = Mapping(
                    subject_id=hgnc_id,
                    subject_label=symbol,
                    object_id=hgnc_id,  # Same entity
                    object_label=alias,
                    predicate_id="oboInOwl:hasRelatedSynonym",
                    mapping_justification=m_meta["mapping_justification"],
                    subject_source=m_meta.get("subject_source"),
                    object_source=m_meta.get("object_source"),
                    mapping_tool=m_meta.get("mapping_tool"),
                    license=m_meta.get("license"),
                    comment="Alias symbol mapping.",
                )
                mappings.append(mapping)

            # Create label mappings for previous symbols
            for prev in prev_symbols:
                mapping = Mapping(
                    subject_id=hgnc_id,
                    subject_label=symbol,
                    object_id=hgnc_id,  # Same entity
                    object_label=prev,
                    predicate_id="oboInOwl:hasRelatedSynonym",
                    mapping_justification=m_meta["mapping_justification"],
                    subject_source=m_meta.get("subject_source"),
                    object_source=m_meta.get("object_source"),
                    mapping_tool=m_meta.get("mapping_tool"),
                    license=m_meta.get("license"),
                    comment="Previous symbol mapping.",
                )
                mappings.append(mapping)

        return mappings

    def _create_mapping_set(
        self, mappings: list[Mapping], mapping_type: str = "id"
    ) -> Sec2PriMappingSet:
        """Create an IdMappingSet or LabelMappingSet with config metadata.

        Delegates to BaseParser.create_mapping_set().
        """
        return self.create_mapping_set(mappings, mapping_type)


__all__ = ["HGNCParser"]
