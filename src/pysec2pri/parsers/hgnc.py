"""HGNC TSV file parser for secondary-to-primary identifier mappings."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import polars as pl
from tqdm import tqdm

from pysec2pri.models import (
    HGNCMapping,
    MappingCardinality,
    MappingSet,
)
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


def _normalize_column_name(col: str) -> str:
    """Normalize column name for case-insensitive matching."""
    return col.lower().strip()


def _find_merged_column(columns: list[str]) -> str | None:
    """Find the merged info column regardless of naming variant."""
    normalized_patterns = [p.lower() for p in MERGED_INFO_PATTERNS]
    for col in columns:
        normalized = _normalize_column_name(col)
        if normalized in normalized_patterns:
            return col
        # Also check for partial match on key identifying part
        if "merged_into_report" in normalized:
            return col
    return None


def _find_column(columns: list[str], name: str) -> str | None:
    """Find column by case-insensitive name."""
    lower_name = name.lower()
    for col in columns:
        if col.lower() == lower_name:
            return col
    return None


def _parse_merged_info(merged_str: str) -> tuple[str, str] | None:
    """Parse merged_into_report to extract hgnc_id and symbol.

    Format varies by HGNC file version:
    - Old: "HGNC:12345/SYMBOL/Approved" (slash-separated)
    - New: "HGNC:12345|SYMBOL|Approved" (pipe-separated)

    Returns (hgnc_id, symbol) or None if parsing fails.
    """
    if not merged_str or merged_str == "":
        return None

    # Try pipe separator first (newer format), then slash (older format)
    if "|" in merged_str:
        parts = merged_str.split("|")
    else:
        parts = merged_str.split("/")

    if len(parts) >= 2:
        return (parts[0].strip(), parts[1].strip())
    return None


def _split_symbols(symbols_str: str) -> list[str]:
    """Split a pipe-separated string of symbols."""
    if not symbols_str:
        return []
    return [s.strip() for s in symbols_str.split("|") if s.strip()]


class HGNCParser(BaseParser):
    """Parser for HGNC TSV files using Polars for memory efficiency.

    Extracts secondary-to-primary HGNC identifier mappings and
    symbol mappings from HGNC withdrawn and complete set files.
    """

    datasource_name = "HGNC"
    withdrawn_source_url = (
        "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/"
        "tsv/withdrawn.txt"
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

    def parse(
        self,
        input_path: Path | str,
        complete_set_path: Path | str | None = None,
    ) -> MappingSet:
        """Parse HGNC TSV files.

        Args:
            input_path: Path to the withdrawn HGNC TSV file.
            complete_set_path: Path to the complete HGNC set TSV file.

        Returns:
            MappingSet with HGNC identifier and symbol mappings.
        """
        input_path = Path(input_path)

        mappings: list[HGNCMapping] = []

        # Parse withdrawn file for ID-to-ID mappings
        withdrawn_mappings = self._parse_withdrawn(input_path)
        mappings.extend(withdrawn_mappings)

        # Parse complete set for symbol mappings if provided
        if complete_set_path:
            complete_set_path = Path(complete_set_path)
            symbol_mappings = self._parse_complete_set(complete_set_path)
            mappings.extend(symbol_mappings)

        mapping_set = MappingSet(
            mappings=mappings,
            datasource_name=self.datasource_name,
            version=self.version,
            mapping_set_id="omicsfixid_hgnc_01",
            mapping_set_description=(
                "Secondary to primary ID mappings for HGNC database, "
                "generated for the omicsFixID project."
            ),
            comment=(
                "object_label represents an alternative symbol "
                "by which the subject_label is also known."
            ),
            curie_map={
                "HGNC:": "https://www.genenames.org/data/gene-symbol-report/"
                "#!/hgnc_id/",
                "IAO:": "http://purl.obolibrary.org/obo/IAO_",
                "oboInOwl:": "http://www.geneontology.org/formats/oboInOwl#",
            },
        )

        return mapping_set

    def _parse_withdrawn(self, file_path: Path) -> list[HGNCMapping]:
        """Parse withdrawn HGNC file for ID-to-ID mappings using Polars."""
        # Read the TSV file with Polars
        df = pl.read_csv(
            file_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=[""],
        )

        # Find columns with case-insensitive matching
        merged_col = _find_merged_column(df.columns)
        if merged_col is None:
            raise ValueError(
                f"Could not find merged_into_report column in {file_path}"
            )

        hgnc_id_col = _find_column(df.columns, HGNC_ID)
        if hgnc_id_col is None:
            raise ValueError(f"Could not find hgnc_id column in {file_path}")

        status_col = _find_column(df.columns, STATUS)
        symbol_col = _find_column(df.columns, "symbol")

        # Build mapping dictionaries for cardinality computation
        primary_counts: dict[str, int] = defaultdict(int)
        secondary_counts: dict[str, int] = defaultdict(int)
        mappings_data: list[tuple[str, str, str]] = []
        withdrawn_ids: list[tuple[str, str | None]] = []  # (id, symbol)

        # Process all rows to separate merged vs withdrawn
        rows = df.iter_rows(named=True)
        if self.show_progress:
            rows = tqdm(rows, total=len(df), desc="Processing withdrawn")

        for row in rows:
            hgnc_id = row[hgnc_id_col]
            merged_info = row.get(merged_col)
            status = row.get(status_col) if status_col else None
            symbol = row.get(symbol_col) if symbol_col else None

            # Check if this is a pure withdrawal (no replacement)
            # HGNC uses "Entry Withdrawn" in status field
            if (not merged_info or merged_info == "") and status:
                if "Entry Withdrawn" in str(status):
                    withdrawn_ids.append((hgnc_id, symbol))
                    continue

            # Handle merged entries
            if merged_info:
                parsed = _parse_merged_info(merged_info)
                if parsed:
                    target_id, target_symbol = parsed
                    primary_counts[target_id] += 1
                    secondary_counts[hgnc_id] += 1
                    mappings_data.append((hgnc_id, target_id, target_symbol))

        # Create mappings with computed cardinality
        mappings: list[HGNCMapping] = []

        # Add withdrawn entries (1:0 cardinality)
        for withdrawn_id, withdrawn_symbol in withdrawn_ids:
            mapping = HGNCMapping(
                primary_id=self.normalize_primary_id(None),
                secondary_id=withdrawn_id,
                primary_symbol=None,
                secondary_symbol=withdrawn_symbol,
                mapping_cardinality=MappingCardinality.ONE_TO_ZERO,
                comment="ID (subject) withdrawn/deprecated.",
                source_url=self.withdrawn_source_url,
            )
            mappings.append(mapping)

        # Add merged entries
        data_iter = mappings_data
        if self.show_progress:
            data_iter = tqdm(mappings_data, desc="Creating ID mappings")

        for source_id, target_id, target_symbol in data_iter:
            # Compute cardinality
            source_count = secondary_counts[source_id]
            target_count = primary_counts[target_id]

            if source_count == 1 and target_count == 1:
                cardinality = MappingCardinality.ONE_TO_ONE
            elif source_count == 1:
                cardinality = MappingCardinality.ONE_TO_MANY
            elif target_count == 1:
                cardinality = MappingCardinality.MANY_TO_ONE
            else:
                cardinality = MappingCardinality.MANY_TO_MANY

            mapping = HGNCMapping(
                primary_id=target_id,
                secondary_id=source_id,
                primary_symbol=target_symbol,
                mapping_cardinality=cardinality,
                secondary_symbol=target_id,  # TODO decide
                source_url=self.withdrawn_source_url,
                comment="object_label represents an alternative symbol by "
                "which the subject_label is also known, though unrecommended.",
            )
            mappings.append(mapping)

        return mappings

    def _parse_complete_set(self, file_path: Path) -> list[HGNCMapping]:
        """Parse complete HGNC set for symbol mappings using Polars."""
        # Read the TSV file with Polars
        df = pl.read_csv(
            file_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=[""],
        )

        # Find columns with case-insensitive matching
        status_col = _find_column(df.columns, STATUS)
        hgnc_id_col = _find_column(df.columns, HGNC_ID)
        symbol_col = _find_column(df.columns, SYMBOL)
        alias_col = _find_column(df.columns, ALIAS_SYMBOL)
        prev_col = _find_column(df.columns, PREV_SYMBOL)

        if status_col is None:
            raise ValueError(f"Could not find status column in {file_path}")
        if hgnc_id_col is None:
            raise ValueError(f"Could not find hgnc_id column in {file_path}")
        if symbol_col is None:
            raise ValueError(f"Could not find symbol column in {file_path}")

        # Filter to approved genes only
        df_approved = df.filter(pl.col(status_col) == "Approved")

        mappings: list[HGNCMapping] = []

        rows = df_approved.iter_rows(named=True)
        if self.show_progress:
            rows = tqdm(
                rows, total=len(df_approved), desc="Processing symbols"
            )

        for row in rows:
            hgnc_id = row[hgnc_id_col]
            symbol = row[symbol_col]

            # Get alias and previous symbols
            alias_str = row.get(alias_col) if alias_col else None
            prev_str = row.get(prev_col) if prev_col else None

            aliases = _split_symbols(alias_str) if alias_str else []
            prev_symbols = _split_symbols(prev_str) if prev_str else []

            has_mappings = bool(aliases or prev_symbols)

            # If include_unmapped_genes and gene has no aliases/prev symbols,
            # still add an entry (matching legacy behavior)
            if self.include_unmapped_genes and not has_mappings:
                mapping = HGNCMapping(
                    primary_id=hgnc_id,
                    secondary_id=None,
                    primary_symbol=symbol,
                    secondary_symbol=symbol,
                    source_url=self.complete_set_source_url,
                )
                mappings.append(mapping)
                continue

            # Add alias symbol mappings
            for alias in aliases:
                mapping = HGNCMapping(
                    primary_id=hgnc_id,
                    secondary_id=None,
                    primary_symbol=symbol,
                    secondary_symbol=alias,
                    source_url=self.complete_set_source_url,
                )
                mappings.append(mapping)

            # Add previous symbol mappings
            for prev in prev_symbols:
                mapping = HGNCMapping(
                    primary_id=hgnc_id,
                    secondary_id=None,
                    primary_symbol=symbol,
                    secondary_symbol=prev,
                    source_url=self.complete_set_source_url,
                )
                mappings.append(mapping)

        return mappings
