"""Base parser class for all datasource parsers."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

from tqdm import tqdm

if TYPE_CHECKING:
    from pysec2pri.models import MappingCardinality, MappingSet

# Import at runtime to avoid circular imports
from pysec2pri.models import BaseMapping


class BaseParser(ABC):
    """Abstract base class for all datasource parsers.

    Each parser is responsible for reading files from a specific datasource
    and extracting secondary-to-primary identifier mappings.
    """

    @staticmethod
    def compute_cardinality(
        pairs: list[tuple[str, str]],
    ) -> dict[tuple[str, str], MappingCardinality]:
        """
        Compute mapping cardinality for a list of (primary, secondary) pairs.

        Args:
            pairs: List of (primary, secondary) tuples.
        Returns:
            Dictionary mapping (primary, secondary) to MappingCardinality.
        """
        from collections import Counter

        from pysec2pri.models import MappingCardinality

        primary_counts: Counter[str] = Counter()
        secondary_counts: Counter[str] = Counter()
        for primary, secondary in pairs:
            primary_counts[primary] += 1
            secondary_counts[secondary] += 1
        result: dict[tuple[str, str], MappingCardinality] = {}
        for primary, secondary in pairs:
            p_count = primary_counts[primary]
            s_count = secondary_counts[secondary]
            if p_count == 1 and s_count == 1:
                cardinality = MappingCardinality.ONE_TO_ONE
            elif p_count == 1:
                cardinality = MappingCardinality.ONE_TO_MANY
            elif s_count == 1:
                cardinality = MappingCardinality.MANY_TO_ONE
            else:
                cardinality = MappingCardinality.MANY_TO_MANY
            result[(primary, secondary)] = cardinality
        return result

    @staticmethod
    def infer_cardinality(primary_count: int, secondary_count: int) -> MappingCardinality:
        """
        Infer mapping cardinality from counts.

        Args:
            primary_count: Number of times the primary appears.
            secondary_count: Number of times the secondary appears.
        Returns:
            MappingCardinality value.
        """
        from pysec2pri.models import MappingCardinality

        if primary_count == 1 and secondary_count == 1:
            return MappingCardinality.ONE_TO_ONE
        elif primary_count == 1:
            return MappingCardinality.ONE_TO_MANY
        elif secondary_count == 1:
            return MappingCardinality.MANY_TO_ONE
        else:
            return MappingCardinality.MANY_TO_MANY

    """Abstract base class for all datasource parsers.

    Each parser is responsible for reading files from a specific datasource
    and extracting secondary-to-primary identifier mappings.
    """

    # To be overridden by subclasses
    datasource_name: str = ""
    default_source_url: str = ""

    def __init__(
        self,
        version: str | None = None,
        show_progress: bool = True,
    ):
        """Initialize the parser.

        Args:
            version: Version/release identifier for the datasource.
            show_progress: Whether to show progress bars during parsing.
        """
        self.version = version
        self.show_progress = show_progress

    @abstractmethod
    def parse(self, input_path: Path | str | None) -> MappingSet:
        """Parse the input file(s) and return a MappingSet.

        Args:
            input_path: Path to the input file or directory.

        Returns:
            A MappingSet containing all extracted mappings.
        """

    def _progress(
        self, iterable: Iterable[Any], desc: str | None = None, total: int | None = None
    ) -> Iterable[Any]:
        """
        Wrap an iterable with a progress bar if enabled.

        Args:
            iterable: The iterable to wrap.
            desc: Description for the progress bar.
            total: Total number of items (if known).
        Returns:
            The iterable, optionally wrapped in tqdm.
        """
        if self.show_progress:
            return cast(Iterable[Any], tqdm(iterable, desc=desc, total=total))
        return iterable

    def _build_comment(
        self,
        base_comment: str,
        additional: str | None = None,
    ) -> str:
        """
        Build a comment string with version information.

        Args:
            base_comment: The base comment text.
            additional: Additional text to append.
        Returns:
            The complete comment string.
        """
        parts = [base_comment] if base_comment else []
        if additional:
            parts.append(additional)
        if self.version:
            parts.append(f"Release: {self.version}.")
        return " ".join(parts)

    def _find_merged_column(
        self,
        columns: list[str],
        merged_info_patterns: list[str],
    ) -> str | None:
        """Find the merged info column regardless of naming variant."""
        normalized_patterns = [p.lower() for p in merged_info_patterns]
        for col in columns:
            normalized = self._normalize_column_name(col)
            if normalized in normalized_patterns:
                return col
            # Also check for partial match on key identifying part
            if "merged_into_report" in normalized:
                return col
        return None

    @staticmethod
    def _normalize_column_name(col: str) -> str:
        """Normalize column name for case-insensitive matching."""
        return col.lower().strip()

    @staticmethod
    def _find_column(columns: list[str], name: str) -> str | None:
        """Find column by case-insensitive name."""
        lower_name = name.lower()
        for col in columns:
            if col.lower() == lower_name:
                return col
        return None

    @staticmethod
    def normalize_subject_id(subject_id: str | None) -> str:
        """
        Normalize a primary ID, converting empty/null values to withdrawn sentinel.

        Args:
            subject_id: The raw primary identifier from the source file.
        Returns:
            The normalized primary ID, or WITHDRAWN_SENTINEL for empty values.
        """
        if not subject_id or subject_id in ("-", ""):
            return BaseMapping.withdrawn_entry
        return subject_id

    @staticmethod
    def _split_symbols(symbols_str: str) -> list[str]:
        """Split a pipe-separated string of symbols."""
        if not symbols_str:
            return []
        return [s.strip() for s in symbols_str.split("|") if s.strip()]

    @staticmethod
    def is_withdrawn_primary(subject_id: str) -> bool:
        """Check if a primary ID represents a withdrawn/deleted entry.

        Args:
            subject_id: The primary identifier to check.

        Returns:
            True if the primary ID indicates a withdrawn entry.
        """
        return subject_id == BaseMapping.withdrawn_entry

    @staticmethod
    def _parse_merged_info(merged_str: str) -> tuple[str, str] | None:
        """Parse merged_into_report to extract hgnc_id and symbol.

        Returns (hgnc_id, symbol) or None if parsing fails.
        """
        if not merged_str or merged_str == "":
            return None
        # Try pipe separator first then slash
        if "|" in merged_str:
            parts = merged_str.split("|")
        else:
            parts = merged_str.split("/")
        if len(parts) >= 2:
            return (parts[0].strip(), parts[1].strip())
        return None


__all__ = ["BaseParser"]
