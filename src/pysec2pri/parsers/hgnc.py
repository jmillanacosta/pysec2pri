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
    BaseDownloader,
    BaseParser,
    LabelMappingSet,
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

    def parse(
        self,
        input_path: Path | str | None,
        complete_set_path: Path | str | None = None,
    ) -> Sec2PriMappingSet:
        """Parse HGNC withdrawn TSV file into an IdMappingSet.

        Args:
            input_path: Path to the withdrawn HGNC TSV file.
            complete_set_path: Optional path to the HGNC complete set TSV.
                When supplied, ``all_primary_ids|symbols`` on the returned mapping set
                is populated with every current HGNC ID, not just those that
                appear as ``object_id`` in a withdrawn to primary mapping.

        Returns:
            IdMappingSet with computed cardinalities based on IDs.
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        input_path = Path(input_path)
        self._resolve_version(input_path)

        # Parse withdrawn file for ID mappings
        mappings = self._parse_withdrawn(input_path)

        # Create IdMappingSet and compute cardinalities
        mapping_set = self._create_mapping_set(mappings, mapping_type="id")

        # Populate the full primary ID set when the complete set is available
        if complete_set_path is not None:
            object.__setattr__(
                mapping_set,
                "_primary_ids",
                self._extract_primary_ids(Path(complete_set_path)),
            )
        return mapping_set

    def parse_primary_symbols(
        self,
        complete_set_path: Path | str | None,
    ) -> Sec2PriMappingSet:
        """Return a mapping set whose only content is the full primary Symbol list.

        Reads the HGNC complete set to extract every current HGNC Symbol and
        stores it in ``_primary_symbols``.  The ``mappings`` list is intentionally
        left empty, this mapping set exists only to drive ``to_pri_symbols()``.

        Args:
            complete_set_path: Path to the HGNC complete set TSV file.

        Returns:
            :class:`~pysec2pri.parsers.base.LabelMappingSet` with no mappings and
            ``_primary_symbols`` populated with all current HGNC symbols.
        """
        if complete_set_path is None:
            raise ValueError("complete_set_path must not be None")
        complete_set_path = Path(complete_set_path)
        self._resolve_version(complete_set_path)

        mapping_set = self._create_mapping_set([], mapping_type="label")
        complete_pri_sym = self._extract_primary_symbols(complete_set_path)
        object.__setattr__(
            mapping_set,
            "_primary_symbols",
            complete_pri_sym,
        )
        return mapping_set

    def parse_primary_ids(
        self,
        complete_set_path: Path | str | None,
    ) -> Sec2PriMappingSet:
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

        mapping_set = self._create_mapping_set([], mapping_type="id")
        object.__setattr__(
            mapping_set,
            "_primary_ids",
            self._extract_primary_ids(complete_set_path),
        )
        return mapping_set

    def parse_symbols(
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

        # Parse complete set for symbol mappings
        mappings = self._parse_complete_set(complete_set_path, statuses=statuses)

        # Create LabelMappingSet and compute cardinalities
        mapping_set = self._create_mapping_set(mappings, mapping_type="label")
        # Set of all primary symbols
        object.__setattr__(
            mapping_set,
            "_primary_symbols",
            self._extract_primary_symbols(Path(complete_set_path)),
        )
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

    def _extract_primary_ids(self, file_path: Path) -> set[str]:
        """Extract all current HGNC IDs from the complete set file.

        Args:
            file_path: Path to the HGNC complete set TSV file.

        Returns:
            Set of all HGNC IDs present in the complete set.
        """
        df = pl.read_csv(
            file_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=[""],
        )
        hgnc_id_col = self._find_column(df.columns, HGNC_ID)
        if hgnc_id_col is None:
            raise ValueError(f"Could not find hgnc_id column in {file_path}")
        return {str(val) for val in df[hgnc_id_col].drop_nulls().to_list()}

    def _extract_primary_symbols(self, file_path: Path) -> dict[str, set[str]]:
        """Extract all current HGNC Symbols from the complete set file.

        Returns a ``dict`` mapping each symbol text to the set of primary HGNC
        IDs that carry that symbol.

        Args:
            file_path: Path to the HGNC complete set TSV file.

        Returns:
            ``dict[symbol, set[hgnc_id]]``
        """
        df = pl.read_csv(
            file_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=[""],
        )
        hgnc_sym_col = self._find_column(df.columns, "symbol")
        hgnc_id_col = self._find_column(df.columns, HGNC_ID)
        if hgnc_sym_col is None:
            raise ValueError(f"Could not find hgnc_symbol column in {file_path}")
        result: dict[str, set[str]] = {}
        for id_, symbol in (
            df.filter(pl.col("status") == "Approved")
            .select([hgnc_id_col, hgnc_sym_col])
            .drop_nulls()
            .rows()
        ):
            result.setdefault(str(symbol), set()).add(str(id_))
        return result

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
        fixed = {
            "mapping_justification": m_meta["mapping_justification"],
            "subject_source": m_meta.get("subject_source"),
            "object_source": m_meta.get("object_source"),
            "mapping_tool": m_meta.get("mapping_tool"),
            "license": m_meta.get("license"),
        }

        rows_data: list[dict[str, str | None]] = []
        for row in df.iter_rows(named=True):
            hgnc_id = row.get(hgnc_id_col)
            if not hgnc_id:
                continue

            merged_info = row.get(merged_col)
            status = row.get(status_col) if status_col else None
            symbol = row.get(symbol_col) if symbol_col else None

            # Case 1: Withdrawn with no replacement
            if not merged_info and status and "Entry Withdrawn" in str(status):
                rows_data.append(
                    {
                        "subject_id": hgnc_id,
                        "object_id": WITHDRAWN_ENTRY,
                        "subject_label": symbol or "",
                        "object_label": WITHDRAWN_ENTRY_LABEL,
                        "predicate_id": "oboInOwl:consider",
                        "comment": "Withdrawn entry with no replacement.",
                    }
                )
                continue

            # Case 2: Merged into another entry
            if merged_info:
                parsed = self._parse_merged_info(merged_info)
                if parsed:
                    target_id, target_symbol = parsed
                    rows_data.append(
                        {
                            "subject_id": hgnc_id,
                            "object_id": target_id,
                            "subject_label": symbol or "",
                            "object_label": target_symbol or "",
                            "predicate_id": m_meta["predicate_id"],
                            "predicate_label": m_meta.get("predicate_label"),
                        }
                    )

        return self._build_mappings(
            rows_data, fixed, desc="Processing withdrawn", total=len(rows_data)
        )

    def _parse_complete_set(
        self, file_path: Path, statuses: list[str] | None = None
    ) -> list[Mapping]:
        """Parse complete HGNC set for symbol (label) mappings.

        Args:
            file_path: Path to the complete HGNC set TSV file.
            statuses: Entry statuses to include (e.g. ``["Approved"]``).
                If ``None`` (default), all entries are included.

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
        assert hgnc_id_col is not None
        assert symbol_col is not None

        # Optionally filter by status
        if statuses is not None and status_col:
            df_approved = df.filter(pl.col(status_col).is_in(statuses))
        else:
            df_approved = df

        m_meta = self.get_mapping_metadata()
        fixed = {
            "mapping_justification": m_meta["mapping_justification"],
            "subject_source": m_meta.get("subject_source"),
            "object_source": m_meta.get("object_source"),
            "mapping_tool": m_meta.get("mapping_tool"),
            "license": m_meta.get("license"),
        }

        rows_data: list[dict[str, str | None]] = []
        for row in df_approved.iter_rows(named=True):
            hgnc_id = row.get(hgnc_id_col)
            symbol = row.get(symbol_col)
            if not hgnc_id or not symbol:
                continue

            alias_str = row.get(alias_col) if alias_col else None
            prev_str = row.get(prev_col) if prev_col else None
            aliases = self._split_symbols(symbols_str=alias_str) if alias_str else []
            prev_symbols = self._split_symbols(symbols_str=prev_str) if prev_str else []

            for alias in aliases:
                rows_data.append(
                    {
                        "subject_id": hgnc_id,
                        "subject_label": alias,
                        "object_id": hgnc_id,
                        "object_label": symbol,
                        "_label_type": "alias",
                        "comment": "Alias symbol mapping.",
                    }
                )

            for prev in prev_symbols:
                rows_data.append(
                    {
                        "subject_id": hgnc_id,
                        "subject_label": prev,
                        "object_id": hgnc_id,
                        "object_label": symbol,
                        "_label_type": "previous",
                        "comment": "Previous symbol mapping.",
                    }
                )

        return self._build_mappings(
            rows_data, fixed, desc="Processing symbols", total=len(rows_data)
        )

    def _create_mapping_set(
        self, mappings: list[Mapping], mapping_type: str = "id"
    ) -> Sec2PriMappingSet:
        """Create an IdMappingSet or LabelMappingSet with config metadata.

        Delegates to BaseParser.create_mapping_set().
        """
        return self.create_mapping_set(mappings, mapping_type)


class HGNCDownloader(BaseDownloader):
    """Downloader for HGNC data files from the quarterly archive."""

    datasource_name = "hgnc"

    def get_download_urls(
        self,
        version: str | None = None,
        **kwargs: object,
    ) -> dict[str, str]:
        """Get HGNC download URLs for *version* (``YYYY-MM-DD``), or latest."""
        from pysec2pri.download import _get_hgnc_urls_for_version, check_hgnc_release

        if version:
            return _get_hgnc_urls_for_version(version)
        return check_hgnc_release().files

    def download(
        self,
        output_dir: Path,
        version: str | None = None,
        decompress: bool = True,
        **kwargs: object,
    ) -> dict[str, Path]:
        """Download HGNC files into *output_dir*."""
        urls = self.get_download_urls(version)
        return self._download_urls(urls, output_dir, decompress)

    def list_versions(self) -> list[str]:
        """List all available HGNC quarterly archive versions.

        Queries the Google Cloud Storage API for all complete-set files and
        returns their dates in ascending order.

        Returns:
            Sorted list of version strings in ``YYYY-MM-DD`` format.
        """
        import re

        import httpx

        gcs_api_url = (
            "https://storage.googleapis.com/storage/v1/b/public-download-files/o"
            "?prefix=hgnc/archive/archive/quarterly/tsv/"
        )
        with httpx.Client(follow_redirects=True, timeout=30.0) as client:
            response = client.get(gcs_api_url)
            response.raise_for_status()
            data = response.json()

        items = data.get("items", [])
        versions: list[str] = []
        for item in items:
            name = item.get("name", "")
            match = re.search(r"hgnc_complete_set_(\d{4}-\d{2}-\d{2})\.txt$", name)
            if match:
                versions.append(match.group(1))
        return sorted(set(versions))


__all__ = ["HGNCDownloader", "HGNCParser"]
