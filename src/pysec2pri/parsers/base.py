"""Base parser class for all datasource parsers."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING

from tqdm import tqdm

if TYPE_CHECKING:
    from pysec2pri.models import MappingSet

# Import at runtime to avoid circular imports
from pysec2pri.models import BaseMapping


class BaseParser(ABC):
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
    def parse(self, input_path: Path | str) -> MappingSet:
        """Parse the input file(s) and return a MappingSet.

        Args:
            input_path: Path to the input file or directory.

        Returns:
            A MappingSet containing all extracted mappings.
        """

    def _progress(self, iterable, desc: str | None = None, total: int | None = None):
        """Wrap an iterable with a progress bar if enabled.

        Args:
            iterable: The iterable to wrap.
            desc: Description for the progress bar.
            total: Total number of items (if known).

        Returns:
            The iterable, optionally wrapped in tqdm.
        """
        if self.show_progress:
            return tqdm(iterable, desc=desc, total=total)
        return iterable

    def _build_comment(
        self,
        base_comment: str,
        additional: str | None = None,
    ) -> str:
        """Build a comment string with version information.

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

    @staticmethod
    def normalize_primary_id(primary_id: str | None) -> str:
        """Normalize a primary ID, converting empty/null values to withdrawn sentinel.

        This method ensures consistent handling of withdrawn/deleted entries
        across all parsers. Empty strings, "-", and None values are all
        converted to the standard WITHDRAWN_SENTINEL value.

        Args:
            primary_id: The raw primary identifier from the source file.

        Returns:
            The normalized primary ID, or WITHDRAWN_SENTINEL for empty values.
        """
        if not primary_id or primary_id in ("-", ""):
            return BaseMapping.WITHDRAWN_SENTINEL
        return primary_id

    @staticmethod
    def is_withdrawn_primary(primary_id: str) -> bool:
        """Check if a primary ID represents a withdrawn/deleted entry.

        Args:
            primary_id: The primary identifier to check.

        Returns:
            True if the primary ID indicates a withdrawn entry.
        """
        return primary_id == BaseMapping.WITHDRAWN_SENTINEL


__all__ = ["BaseParser"]
