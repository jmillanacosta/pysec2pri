"""HGNC TSV file parser for secondary-to-primary identifier mappings.

This parser extracts:
1. ID-to-ID mappings: withdrawn/merged HGNC IDs -> current HGNC IDs
2. Label-to-label mappings: previous/alias labels -> current labels

Uses SSSOM-compliant MappingSet classes with cardinality computation.
"""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

import polars as pl
from sssom_schema import Mapping

from pysec2pri.parsers._nomenclature import (
    ALIAS_SYMBOL,
    DATE_SYMBOL_CHANGED,
    PREV_SYMBOL,
    STATUS,
    SYMBOL,
    GeneNomenclatureParser,
)
from pysec2pri.parsers.base import BaseMappingSet, LabelMappingSet

# HGNC column names (case-insensitive matching used)
HGNC_ID = "hgnc_id"

# Merged info column has different naming variants across HGNC file versions
MERGED_INFO_PATTERNS = [
    "merged_into_report(i.e. hgnc_id/symbol/status)",
    "merged_into_report(i.e hgnc_id/symbol/status)",
    "merged_into_report(s) (i.e hgnc_id|symbol|status)",
]


class HGNCParser(GeneNomenclatureParser):
    """Parser for HGNC TSV files using Polars for memory efficiency.

    Extracts secondary-to-primary HGNC identifier mappings and
    symbol mappings from HGNC withdrawn and complete set files.

    Returns:
    - IdMappingSet for ID-to-ID mappings (withdrawn/merged IDs)
    - LabelMappingSet for symbol mappings (alias/previous symbols)
    """

    datasource_name = "hgnc"
    id_column: ClassVar[str] = HGNC_ID
    withdrawn_label_column: ClassVar[str] = SYMBOL
    merged_patterns: ClassVar[list[str]] = MERGED_INFO_PATTERNS

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

    def parse(
        self,
        input_path: Path | str | None,
        complete_set_path: Path | str | None = None,
    ) -> BaseMappingSet:
        """Parse HGNC withdrawn TSV file into an IdMappingSet.

        Args:
            input_path: Path to the withdrawn HGNC TSV file.
            complete_set_path: Optional path to the HGNC complete set TSV.
                When supplied, ``all_primary_ids|labels`` on the returned mapping set
                is populated with every current HGNC ID, not just those that
                appear as ``object_id`` in a withdrawn to primary mapping.

        Returns:
            IdMappingSet with computed cardinalities based on IDs.
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        input_path = Path(input_path)
        self._resolve_version(input_path)

        mappings = self._parse_withdrawn(input_path)
        # Populate the full primary ID set when the complete set is available
        primary_ids = (
            self._extract_primary_ids(Path(complete_set_path))
            if complete_set_path is not None
            else None
        )
        return self.create_mapping_set(mappings, mapping_type="id", primary_ids=primary_ids)

    def parse_primary_labels(
        self,
        complete_set_path: Path | str | None,
    ) -> BaseMappingSet:
        """Return a mapping set whose only content is the full primary Symbol list.

        Reads the HGNC complete set to extract every current HGNC Symbol and
        stores it in ``_primary_labels``.  The ``mappings`` list is intentionally
        left empty, this mapping set exists only to drive ``to_pri_labels()``.

        Args:
            complete_set_path: Path to the HGNC complete set TSV file.

        Returns:
            :class:`~pysec2pri.parsers.base.LabelMappingSet` with no mappings and
            ``_primary_labels`` populated with all current HGNC labels.
        """
        if complete_set_path is None:
            raise ValueError("complete_set_path must not be None")
        complete_set_path = Path(complete_set_path)
        self._resolve_version(complete_set_path)

        return self.create_mapping_set(
            [],
            mapping_type="label",
            primary_labels=self._extract_primary_labels(complete_set_path),
        )

    def parse_primary_ids(
        self,
        complete_set_path: Path | str | None,
    ) -> BaseMappingSet:
        """Return a mapping set whose only content is the full primary ID list.

        Reads the HGNC complete set to extract every current HGNC ID and
        stores it in ``_primary_ids``.  The ``mappings`` list is intentionally
        left empty, this mapping set exists only to drive ``to_pri_ids()``.

        Args:
            complete_set_path: Path to the HGNC complete set TSV file.

        Returns:
            :class:`~pysec2pri.parsers.base.IdMappingSet` with no mappings and
            ``_primary_ids`` populated with all current HGNC IDs.
        """
        if complete_set_path is None:
            raise ValueError("complete_set_path must not be None")
        complete_set_path = Path(complete_set_path)
        self._resolve_version(complete_set_path)

        return self.create_mapping_set(
            [],
            mapping_type="id",
            primary_ids=self._extract_primary_ids(complete_set_path),
        )

    def parse_labels(
        self,
        complete_set_path: Path | str | None,
        statuses: list[str] | None = None,
    ) -> LabelMappingSet:
        """Parse HGNC complete set for symbol (label) mappings.

        Args:
            complete_set_path: Path to the complete HGNC set TSV file.
            statuses: Entry statuses to include (e.g. ``["Approved"]``).
                If ``None`` (default), all entries are included.

        Returns:
            LabelMappingSet with computed cardinalities based on labels.
        """
        if complete_set_path is None:
            raise ValueError("complete_set_path must not be None")
        complete_set_path = Path(complete_set_path)
        self._resolve_version(complete_set_path)

        mappings = self._parse_complete_set(complete_set_path, statuses=statuses)
        return self.create_mapping_set(
            mappings,
            mapping_type="label",
            primary_labels=self._extract_primary_labels(complete_set_path),
        )

    def parse_all(
        self,
        withdrawn_path: Path | str | None,
        complete_set_path: Path | str | None,
    ) -> tuple[BaseMappingSet, BaseMappingSet]:
        """Parse both withdrawn and complete set files.

        Args:
            withdrawn_path: Path to the withdrawn HGNC TSV file.
            complete_set_path: Path to the complete HGNC set TSV file.

        Returns:
            Tuple of (IdMappingSet, LabelMappingSet).
        """
        id_mappings = self.parse(withdrawn_path)
        label_mappings = self.parse_labels(complete_set_path)
        return id_mappings, label_mappings

    def _extract_primary_labels(self, file_path: Path) -> dict[str, set[str]]:
        """Extract all current HGNC Symbols from the complete set file.

        Returns a ``dict`` mapping each symbol text to the set of primary HGNC
        IDs that carry that symbol.

        Args:
            file_path: Path to the HGNC complete set TSV file.

        Returns:
            ``dict[label, set[hgnc_id]]``
        """
        df = self._read_tsv(file_path)
        sym_col = self._find_column(df.columns, SYMBOL)
        id_col = self._find_column(df.columns, self.id_column)
        status_col = self._find_column(df.columns, STATUS)
        if sym_col is None:
            raise ValueError(f"Could not find hgnc_symbol column in {file_path}")
        if id_col is None:
            raise ValueError(f"Could not find hgnc_id column in {file_path}")
        approved = df.filter(pl.col(status_col) == "Approved") if status_col else df
        return self._labels_by_id(approved, id_col, sym_col)

    def _parse_complete_set(
        self, file_path: Path, statuses: list[str] | None = None
    ) -> list[Mapping]:
        """Parse complete HGNC set for symbol (label) mappings.

        Args:
            file_path: Path to the complete HGNC set TSV file.
            statuses: Entry statuses to include (e.g. ``["Approved"]``).
                If ``None`` (default), all entries are included.

        Returns:
            List of SSSOM Mapping objects for label mappings.
        """
        df = self._read_tsv(file_path)

        status_col = self._find_column(df.columns, STATUS)
        id_col = self._find_column(df.columns, self.id_column)
        label_col = self._find_column(df.columns, SYMBOL)
        alias_col = self._find_column(df.columns, ALIAS_SYMBOL)
        prev_col = self._find_column(df.columns, PREV_SYMBOL)
        date_changed_col = self._find_column(df.columns, DATE_SYMBOL_CHANGED)

        if not all([status_col, id_col, label_col]):
            raise ValueError(f"Missing required columns in {file_path}")
        assert id_col is not None
        assert label_col is not None

        if statuses is not None and status_col:
            df = df.filter(pl.col(status_col).is_in(statuses))

        return self._build_label_mappings(
            df,
            id_col=id_col,
            label_col=label_col,
            alias_col=alias_col,
            prev_col=prev_col,
            date_changed_col=date_changed_col,
        )


__all__ = ["HGNCParser"]
