"""A simple copy-paste parser example. Copy to ``src/pysec2pri/parsers/<source>.py``.

This example uses two files per release:

- ``withdrawn.tsv``: ``example_id, replacement_id, status``.
- ``complete.tsv``: ``example_id, symbol, prev_symbols, alias_symbols``.

The class name must match ``parser_class`` in the config, and every method the
config's ``mapping_sets`` names must exist here. The specific published files
for a datasource will differ, adapt as needed.
"""

from __future__ import annotations

from pathlib import Path

from pysec2pri.parsers.base import (
    WITHDRAWN_ENTRY,
    WITHDRAWN_ENTRY_LABEL,
    BaseMappingSet,
    BaseParser,
)

# Column names in Example's files. _find_column matches case-insensitively.
EXAMPLE_ID = "example_id"
REPLACEMENT_ID = "replacement_id"
STATUS = "status"
SYMBOL = "symbol"
PREV_SYMBOLS = "prev_symbols"
ALIAS_SYMBOLS = "alias_symbols"


class ExampleParser(BaseParser):
    """Parser for Example's withdrawn and complete files."""

    datasource_name = "example"

    def parse(
        self,
        input_path: Path | str | None,
        complete_set_path: Path | str | None = None,
    ) -> BaseMappingSet:
        """Parse the withdrawn file into retired -> current ID mappings.

        Args:
            input_path: Path to ``withdrawn.tsv``.
            complete_set_path: Path to ``complete.tsv``. Optional: without it
                the mappings are the same, but the full current-ID list and
                ambiguity checks are incomplete.

        Returns:
            An ``IdMappingSet``.

        Raises:
            ValueError: If *input_path* is missing, or lacks its ID column.
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        input_path = Path(input_path)
        self._resolve_version(input_path)

        df = self._read_tsv(input_path)
        id_col = self._find_column(df.columns, EXAMPLE_ID)
        if id_col is None:
            raise ValueError(f"No {EXAMPLE_ID} column in {input_path}")
        new_col = self._find_column(df.columns, REPLACEMENT_ID)
        status_col = self._find_column(df.columns, STATUS)

        # Fields shared by every row (license, sources, justification ...).
        fixed = self._fixed_mapping_fields()
        m_meta = self.get_mapping_metadata()

        rows: list[dict[str, str | None]] = []
        for row in df.iter_rows(named=True):
            retired = row.get(id_col)
            if not retired:
                continue
            replacement = row.get(new_col) if new_col else None
            status = str(row.get(status_col) or "") if status_col else ""

            if replacement:
                rows.append(
                    {
                        "subject_id": retired,
                        "object_id": replacement,
                        "predicate_id": m_meta["predicate_id"],
                        "predicate_label": m_meta.get("predicate_label"),
                    }
                )
            elif "withdrawn" in status.lower():
                # Retired with nothing to point at.
                rows.append(
                    {
                        "subject_id": retired,
                        "object_id": WITHDRAWN_ENTRY,
                        "object_label": WITHDRAWN_ENTRY_LABEL,
                        "predicate_id": "oboInOwl:consider",
                        "comment": "Withdrawn entry with no replacement.",
                    }
                )

        mappings = self._build_mappings(rows, fixed, desc="Processing withdrawn")
        primary_ids = (
            self._extract_primary_ids(Path(complete_set_path)) if complete_set_path else None
        )
        return self.create_mapping_set(mappings, mapping_type="id", primary_ids=primary_ids)

    def parse_labels(self, complete_set_path: Path | str | None) -> BaseMappingSet:
        """Parse the complete file into old -> current symbol mappings.

        Each current entry lists the symbols it used to have (``prev_symbols``)
        and the ones it also goes by (``alias_symbols``). ``_label_type`` picks
        the predicate: ``previous`` renames, ``alias`` synonyms.

        Args:
            complete_set_path: Path to ``complete.tsv``.

        Returns:
            A ``LabelMappingSet``.

        Raises:
            ValueError: If *complete_set_path* is missing, or lacks its columns.
        """
        if complete_set_path is None:
            raise ValueError("complete_set_path must not be None")
        complete_set_path = Path(complete_set_path)
        self._resolve_version(complete_set_path)

        df = self._read_tsv(complete_set_path)
        id_col = self._find_column(df.columns, EXAMPLE_ID)
        symbol_col = self._find_column(df.columns, SYMBOL)
        if id_col is None or symbol_col is None:
            raise ValueError(f"No {EXAMPLE_ID}/{SYMBOL} column in {complete_set_path}")
        prev_col = self._find_column(df.columns, PREV_SYMBOLS)
        alias_col = self._find_column(df.columns, ALIAS_SYMBOLS)

        fixed = self._fixed_mapping_fields()

        rows: list[dict[str, str | None]] = []
        for row in df.iter_rows(named=True):
            current_id, current = row.get(id_col), row.get(symbol_col)
            if not current_id or not current:
                continue
            for col, label_type in ((prev_col, "previous"), (alias_col, "alias")):
                raw = row.get(col) if col else None
                for old in self._split_labels(labels_str=str(raw)) if raw else []:
                    rows.append(
                        {
                            "subject_label": old,
                            "object_id": current_id,
                            "object_label": current,
                            "_label_type": label_type,
                        }
                    )

        mappings = self._build_mappings(rows, fixed, desc="Processing symbols")
        return self.create_mapping_set(
            mappings,
            mapping_type="label",
            primary_labels=self._extract_primary_labels(complete_set_path),
        )

    def parse_primary_ids(self, complete_set_path: Path | str | None) -> BaseMappingSet:
        """Return a mapping set holding only the full current-ID list.

        Args:
            complete_set_path: Path to ``complete.tsv``.

        Returns:
            A mapping set with no mappings and ``_primary_ids`` populated.

        Raises:
            ValueError: If *complete_set_path* is missing.
        """
        if complete_set_path is None:
            raise ValueError("complete_set_path must not be None")
        complete_set_path = Path(complete_set_path)
        self._resolve_version(complete_set_path)
        return self.create_mapping_set(
            [], mapping_type="id", primary_ids=self._extract_primary_ids(complete_set_path)
        )

    def _extract_primary_ids(self, file_path: Path) -> set[str]:
        """Return every current ID in the complete file.

        Raises:
            ValueError: If the file lacks its ID column.
        """
        df = self._read_tsv(file_path)
        id_col = self._find_column(df.columns, EXAMPLE_ID)
        if id_col is None:
            raise ValueError(f"No {EXAMPLE_ID} column in {file_path}")
        return {str(v) for v in df[id_col].drop_nulls().to_list()}

    def _extract_primary_labels(self, file_path: Path) -> dict[str, set[str]]:
        """Return ``{current symbol: {ID, ...}}`` for the complete file.

        Raises:
            ValueError: If the file lacks its ID or symbol column.
        """
        df = self._read_tsv(file_path)
        id_col = self._find_column(df.columns, EXAMPLE_ID)
        symbol_col = self._find_column(df.columns, SYMBOL)
        if id_col is None or symbol_col is None:
            raise ValueError(f"No {EXAMPLE_ID}/{SYMBOL} column in {file_path}")
        result: dict[str, set[str]] = {}
        for id_, label in df.select([id_col, symbol_col]).drop_nulls().rows():
            result.setdefault(str(label), set()).add(str(id_))
        return result


__all__ = ["ExampleParser"]
