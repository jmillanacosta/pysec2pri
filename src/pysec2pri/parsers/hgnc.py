"""HGNC TSV file parser for secondary-to-primary identifier mappings."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import polars as pl

from pysec2pri.models import BaseMapping, HGNCMapping, MappingSet
from pysec2pri.parsers.base import BaseParser

# HGNC column names (case-insensitive matching used)
HGNC_ID = "hgnc_id"
SYMBOL = "symbol"
ALIAS_SYMBOL = "alias_symbol"
PREV_SYMBOL = "prev_symbol"
STATUS = "status"

# Merged info column has different naming variants across HGNC file versions
# Format evolved over time - need to handle multiple patterns
MERGED_INFO_PATTERNS = [
    "merged_into_report(i.e. hgnc_id/symbol/status)",
    "merged_into_report(i.e hgnc_id/symbol/status)",
    "merged_into_report(s) (i.e hgnc_id|symbol|status)",
]


class HGNCParser(BaseParser):
    """Parser for HGNC TSV files using Polars for memory efficiency.

    Extracts secondary-to-primary HGNC identifier mappings and
    symbol mappings from HGNC withdrawn and complete set files.
    """

    datasource_name = "HGNC"
    withdrawn_source_url = (
        "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/withdrawn.txt"
    )
    complete_set_source_url = (
        "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/"
        "tsv/hgnc_complete_set.txt"
    )

    def __init__(
        self,
        version: str | None = None,
        show_progress: bool = True,
        include_unmapped_genes: bool = False,
    ):
        """Initialize the HGNC parser.

        Args:
            version: Version/release identifier for the datasource.
            show_progress: Whether to show progress bars during parsing.
            include_unmapped_genes: If True, include entries for genes that
                have no alias or previous symbols (just primary ID mapping
                to itself). This matches the legacy R script behavior.
        """
        super().__init__(version=version, show_progress=show_progress)
        self.include_unmapped_genes = include_unmapped_genes

    def _read_withdrawn_file(self, file_path: Path) -> pl.DataFrame:
        """Read the withdrawn HGNC file into a Polars DataFrame."""
        return pl.read_csv(
            file_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=[""],
        )

    from collections.abc import Iterable

    def _process_withdrawn_rows(
        self,
        rows: Iterable[dict[str, str]],
        merged_col: str,
        hgnc_id_col: str,
        status_col: str | None,
        symbol_col: str | None,
    ) -> tuple[
        list[tuple[str, str | None]],
        list[tuple[str, str, str]],
        dict[str, int],
        dict[str, int],
    ]:
        primary_counts: dict[str, int] = defaultdict(int)
        secondary_counts: dict[str, int] = defaultdict(int)
        mappings_data: list[tuple[str, str, str]] = []
        withdrawn_ids: list[tuple[str, str | None]] = []

        rows = list(rows)
        rows = list(self._progress(rows, desc="Processing withdrawn"))

        for row in rows:
            hgnc_id = self.normalize_subject_id(row[hgnc_id_col])
            merged_info = row.get(merged_col)
            status = row.get(status_col) if status_col else None
            symbol = row.get(symbol_col) if symbol_col else None

            if (not merged_info or merged_info == "") and status:
                if "Entry Withdrawn" in str(status):
                    withdrawn_ids.append((hgnc_id, symbol))
                    continue

            if merged_info:
                parsed = self._parse_merged_info(merged_info)
                if parsed:
                    target_id, target_symbol = parsed
                    target_id = self.normalize_subject_id(target_id)
                    primary_counts[target_id] += 1
                    secondary_counts[hgnc_id] += 1
                    mappings_data.append((hgnc_id, target_id, target_symbol))

        return withdrawn_ids, mappings_data, primary_counts, secondary_counts

    def _create_withdrawn_mappings(
        self,
        withdrawn_ids: list[tuple[str, str | None]],
    ) -> list[HGNCMapping]:
        """Create mappings for withdrawn entries."""
        mappings = []
        for withdrawn_id, withdrawn_symbol in withdrawn_ids:
            mapping = HGNCMapping(
                subject_id=BaseMapping.withdrawn_entry,
                predicate_id="oboInOwl:consider",
                object_id=withdrawn_id,
                subject_label="",  # No label for withdrawn entry
                object_label=withdrawn_symbol or "",
                comment="ID (subject) withdrawn/deprecated.",
                source_url=self.withdrawn_source_url,
            )
            mappings.append(mapping)
        return mappings

    def _create_merged_mappings(
        self,
        mappings_data: list[tuple[str, str, str]],
        primary_counts: dict[str, int],
        secondary_counts: dict[str, int],
    ) -> list[HGNCMapping]:
        pairs = [(source_id, target_id) for source_id, target_id, _ in mappings_data]
        self.compute_cardinality(pairs)
        data_iter = self._progress(mappings_data, desc="Creating ID mappings")
        mappings = []
        for source_id, target_id, target_symbol in data_iter:
            mapping = HGNCMapping(
                subject_id=target_id,  # CURIE
                predicate_id="IAO:0100001",
                object_id=source_id,  # CURIE
                subject_label=target_symbol or "",  # symbol/label for subject
                object_label="",  # no label for object (ID mapping)
                source_url=self.withdrawn_source_url,
                comment=self._build_comment(
                    "object_label represents an alternative symbol by which "
                    "the subject_label is also known, though unrecommended."
                ),
            )
            mappings.append(mapping)
        return mappings

    def _parse_withdrawn(self, file_path: Path) -> list[HGNCMapping]:
        """Parse withdrawn HGNC file for ID-to-ID mappings using Polars."""
        df = self._read_withdrawn_file(file_path)

        merged_col = self._find_merged_column(df.columns, MERGED_INFO_PATTERNS)
        if merged_col is None:
            raise ValueError(f"Could not find merged_into_report column in {file_path}")

        hgnc_id_col = self._find_column(df.columns, HGNC_ID)
        if hgnc_id_col is None:
            raise ValueError(f"Could not find hgnc_id column in {file_path}")

        status_col = self._find_column(df.columns, STATUS)
        symbol_col = self._find_column(df.columns, "symbol")

        rows = df.iter_rows(named=True)
        (withdrawn_ids, mappings_data, primary_counts, secondary_counts) = (
            self._process_withdrawn_rows(rows, merged_col, hgnc_id_col, status_col, symbol_col)
        )

        mappings = []
        mappings.extend(self._create_withdrawn_mappings(withdrawn_ids))
        mappings.extend(
            self._create_merged_mappings(mappings_data, primary_counts, secondary_counts)
        )

        return mappings

    def _parse_complete_set(self, file_path: Path) -> list[HGNCMapping]:
        """Parse complete HGNC set for symbol mappings using Polars."""
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
        if status_col is None:
            raise ValueError(f"Could not find status column in {file_path}")
        if hgnc_id_col is None:
            raise ValueError(f"Could not find hgnc_id column in {file_path}")
        if symbol_col is None:
            raise ValueError(f"Could not find symbol column in {file_path}")
        df_approved = df.filter(pl.col(status_col) == "Approved")
        mappings: list[HGNCMapping] = []
        rows_iter = df_approved.iter_rows(named=True)
        rows = self._progress(rows_iter, desc="Processing symbols", total=len(df_approved))
        for row in rows:
            hgnc_id = self.normalize_subject_id(row[hgnc_id_col])
            symbol = row[symbol_col]
            alias_str = row.get(alias_col) if alias_col else None
            prev_str = row.get(prev_col) if prev_col else None
            aliases = self._split_symbols(alias_str) if alias_str else []
            prev_symbols = self._split_symbols(prev_str) if prev_str else []
            has_mappings = bool(aliases or prev_symbols)
            if self.include_unmapped_genes and not has_mappings:
                mapping = HGNCMapping(
                    subject_id=hgnc_id,  # CURIE
                    predicate_id="owl:sameAs",
                    object_id=hgnc_id,  # CURIE
                    subject_label=symbol or "",  # label
                    object_label="",  # no label for object
                    source_url=self.complete_set_source_url,
                    comment=self._build_comment("No secondaries for gene entry."),
                )
                mappings.append(mapping)
                continue
            for alias in aliases:
                mapping = HGNCMapping(
                    subject_id=hgnc_id,  # CURIE
                    predicate_id="oboInOwl:consider",
                    object_id=hgnc_id,  # CURIE (not alias!)
                    subject_label=symbol or "",  # label
                    object_label=alias or "",  # alias as label
                    source_url=self.complete_set_source_url,
                    comment=self._build_comment("Alias symbol mapping."),
                )
                mappings.append(mapping)
            for prev in prev_symbols:
                mapping = HGNCMapping(
                    subject_id=hgnc_id,  # CURIE
                    predicate_id="oboInOwl:consider",
                    object_id=hgnc_id,  # CURIE (not prev symbol!)
                    subject_label=symbol or "",  # label
                    object_label=prev or "",  # prev symbol as label
                    source_url=self.complete_set_source_url,
                    comment=self._build_comment("Previous symbol mapping."),
                )
                mappings.append(mapping)
        return mappings

    def parse(
        self,
        input_path: Path | str | None,
        complete_set_path: Path | str | None = None,
    ) -> MappingSet:
        """
        Parse HGNC TSV files and return a MappingSet.

        Args:
            input_path: Path to the withdrawn HGNC TSV file.
            complete_set_path: Optional path to the complete HGNC set TSV file.

        Returns:
            MappingSet containing HGNC identifier and symbol mappings.
        """
        return self._parse_with_optional_complete_set(input_path, complete_set_path)

    def _parse_with_optional_complete_set(
        self,
        input_path: Path | str | None,
        complete_set_path: Path | str | None = None,
    ) -> MappingSet:
        if input_path is None:
            raise ValueError("input_path must not be None")
        input_path = Path(input_path)
        mappings: list[HGNCMapping] = []
        withdrawn_mappings = self._parse_withdrawn(input_path)
        mappings.extend(withdrawn_mappings)
        if complete_set_path:
            complete_set_path = Path(complete_set_path)
            symbol_mappings = self._parse_complete_set(complete_set_path)
            mappings.extend(symbol_mappings)
        mapping_set = MappingSet(
            mappings=mappings,
            datasource_name=self.datasource_name,
            version=self.version,
            mapping_set_id="omicsfixid_hgnc_01",
            mapping_set_description="Secondary to primary ID mappings for HGNC"
            " database, generated for the omicsFixID project.",
            comment=self._build_comment(
                "object_label represents an alternative symbol by which the"
                "subject_label is also known."
            ),
            curie_map={
                "HGNC:": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",
                "IAO:": "http://purl.obolibrary.org/obo/IAO_",
                "oboInOwl:": "http://www.geneontology.org/formats/oboInOwl#",
            },
        )
        return mapping_set
