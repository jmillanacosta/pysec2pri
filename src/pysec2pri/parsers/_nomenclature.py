"""Shared base for gene-nomenclature parsers (HGNC, VGNC, PGNC?).

Both committees publish the same two-file shape, a withdrawn file of
retired/merged identifiers and a complete gene-set file with the current
symbols with their alias and previous symbols. This base holds the parsing
they have in common; each subclass will specify the identifier column,
withdrawn-symbol column, and merged-info column patterns, and adds whatever
is specific to it e.g. VGNC taxons.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar

import polars as pl
from sssom_schema import Mapping

from pysec2pri.parsers.base import (
    WITHDRAWN_ENTRY,
    WITHDRAWN_ENTRY_LABEL,
    BaseParser,
)

if TYPE_CHECKING:
    from pathlib import Path

# Column names shared across both gene-nomenclature files (case-insensitive).
SYMBOL = "symbol"
ALIAS_SYMBOL = "alias_symbol"
PREV_SYMBOL = "prev_symbol"
DATE_SYMBOL_CHANGED = "date_symbol_changed"
STATUS = "status"


class GeneNomenclatureParser(BaseParser):
    """Base parser for gene-nomenclature datasources that publish HGNC-shaped files.

    Subclasses set :attr:`id_column`, :attr:`withdrawn_label_column`, and
    :attr:`merged_patterns` to bind the shared parsing to their own columns.
    """

    #: Column holding the primary identifier (e.g. ``"hgnc_id"``).
    id_column: ClassVar[str]
    #: Column holding a withdrawn entry symbol.
    withdrawn_label_column: ClassVar[str]
    #: Naming variants of the merged-info column across file versions.
    merged_patterns: ClassVar[list[str]] = []

    def _extract_primary_ids(self, file_path: Path) -> set[str]:
        """Extract every current identifier from the complete gene-set file.

        Args:
            file_path: Path to the complete gene-set TSV file.

        Returns:
            Set of all identifiers present in the file.
        """
        df = self._read_tsv(file_path)
        id_col = self._find_column(df.columns, self.id_column)
        if id_col is None:
            raise ValueError(f"Could not find {self.id_column} column in {file_path}")
        return {str(val) for val in df[id_col].drop_nulls().to_list()}

    def _labels_by_id(self, df: pl.DataFrame, id_col: str, symbol_col: str) -> dict[str, set[str]]:
        """Return ``{symbol: {identifier, ...}}`` for every row of *df*."""
        result: dict[str, set[str]] = {}
        for id_, label in df.select([id_col, symbol_col]).drop_nulls().rows():
            result.setdefault(str(label), set()).add(str(id_))
        return result

    def _parse_withdrawn(
        self, file_path: Path, *, taxon_by_id: dict[str, str] | None = None
    ) -> list[Mapping]:
        """Parse a withdrawn file into ID-to-ID mappings.

        Args:
            file_path: Path to the withdrawn TSV file.
            taxon_by_id: Optional ``{identifier: taxon_id}`` map.

        Returns:
            List of SSSOM Mapping objects.
        """
        df = self._read_tsv(file_path)

        merged_col = self._find_merged_column(df.columns, self.merged_patterns)
        if merged_col is None:
            raise ValueError(f"Could not find merged_into_report column in {file_path}")

        id_col = self._find_column(df.columns, self.id_column)
        if id_col is None:
            raise ValueError(f"Could not find {self.id_column} column in {file_path}")

        status_col = self._find_column(df.columns, STATUS)
        label_col = self._find_column(df.columns, self.withdrawn_label_column)

        fixed = self._fixed_mapping_fields()

        rows_data: list[dict[str, str | None]] = []
        for row in df.iter_rows(named=True):
            identifier = row.get(id_col)
            if not identifier:
                continue

            merged_info = row.get(merged_col)
            status = row.get(status_col) if status_col else None
            label = row.get(label_col) if label_col else None

            # Withdrawn with no replacement
            if not merged_info and status and "Entry Withdrawn" in str(status):
                rows_data.append(
                    {
                        "subject_id": identifier,
                        "object_id": WITHDRAWN_ENTRY,
                        "subject_label": label or "",
                        "object_label": WITHDRAWN_ENTRY_LABEL,
                        "predicate_id": "oboInOwl:consider",
                        "comment": "Withdrawn entry with no replacement.",
                        "record_id": self._record_id(
                            self._record_namespace(),
                            WITHDRAWN_ENTRY,
                            identifier,
                        ),
                    }
                )
                continue

            # Merged into another entry
            if merged_info:
                parsed = self._parse_merged_info(merged_info)
                if parsed:
                    target_id, target_label = parsed
                    m_meta = self.get_mapping_metadata()
                    row_taxon = taxon_by_id.get(target_id) if taxon_by_id else None
                    record_ns = (
                        self._record_namespace(species=row_taxon)
                        if row_taxon is not None
                        else self._record_namespace()
                    )
                    rows_data.append(
                        {
                            "subject_id": identifier,
                            "object_id": target_id,
                            "subject_label": label or "",
                            "object_label": target_label or "",
                            "predicate_id": m_meta["predicate_id"],
                            "predicate_label": m_meta.get("predicate_label"),
                            "record_id": self._record_id(record_ns, target_id, identifier),
                        }
                    )

        return self._build_mappings(
            rows_data, fixed, desc="Processing withdrawn", total=len(rows_data)
        )

    def _build_label_mappings(
        self,
        df: pl.DataFrame,
        *,
        id_col: str,
        label_col: str,
        alias_col: str | None,
        prev_col: str | None,
        date_changed_col: str | None,
        taxon_col: str | None = None,
    ) -> list[Mapping]:
        """Build alias/previous symbol label mappings from a filtered gene-set frame.

        Args:
            df: The gene-set frame, already filtered to the wanted rows.
            id_col: Resolved primary-identifier column.
            label_col: Resolved current-symbol column.
            alias_col: Resolved alias-symbol column, or ``None``.
            prev_col: Resolved previous-symbol column, or ``None``.
            date_changed_col: Resolved symbol-changed-date column, or ``None``.
            taxon_col: Optional resolved per-row taxon column name.

        Returns:
            List of SSSOM Mapping objects for label mappings.
        """
        fixed = self._fixed_mapping_fields()

        rows_data: list[dict[str, str | None]] = []
        for row in df.iter_rows(named=True):
            identifier = row.get(id_col)
            label = row.get(label_col)
            if not identifier or not label:
                continue

            alias_str = row.get(alias_col) if alias_col else None
            prev_str = row.get(prev_col) if prev_col else None
            aliases = self._split_labels(labels_str=alias_str) if alias_str else []
            prev_labels = self._split_labels(labels_str=prev_str) if prev_str else []
            # The date the current symbol was set, i.e. when its previous
            # symbol(s) became secondary. Only the most recent change is
            # recorded, so with multiple prev_symbol entries this date applies
            # exactly to the latest rename and approximately to earlier ones.
            symbol_changed_date = row.get(date_changed_col) if date_changed_col else None
            row_taxon = row.get(taxon_col) if taxon_col else None
            record_ns = (
                self._record_namespace(species=str(row_taxon))
                if row_taxon is not None
                else self._record_namespace()
            )

            for alias in aliases:
                rows_data.append(
                    {
                        "object_id": identifier,
                        "subject_label": alias,
                        "subject_type": "rdfs literal",
                        "object_label": label,
                        "_label_type": "alias",
                        "comment": "Alias symbol mapping.",
                        "record_id": self._record_id(record_ns, identifier, alias),
                    }
                )

            for prev in prev_labels:
                rows_data.append(
                    {
                        "object_id": identifier,
                        "subject_label": prev,
                        "subject_type": "rdfs literal",
                        "object_label": label,
                        "_label_type": "previous",
                        "comment": "Previous symbol mapping.",
                        "mapping_date": symbol_changed_date,
                        "record_id": self._record_id(record_ns, identifier, prev),
                    }
                )

        return self._build_mappings(
            rows_data, fixed, desc="Processing symbols", total=len(rows_data)
        )


__all__ = [
    "ALIAS_SYMBOL",
    "DATE_SYMBOL_CHANGED",
    "PREV_SYMBOL",
    "STATUS",
    "SYMBOL",
    "GeneNomenclatureParser",
]
