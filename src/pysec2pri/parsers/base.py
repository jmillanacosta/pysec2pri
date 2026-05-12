"""Base parser class for all datasource parsers."""

from __future__ import annotations

import re
from abc import ABC, abstractmethod
from collections.abc import Iterable
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any, cast

import yaml
from sssom_schema import Mapping, MappingCardinalityEnum, MappingSet
from tqdm import tqdm

from pysec2pri.logging import logger
from pysec2pri.version import VERSION

# Values for withdrawn entries

WITHDRAWN_ENTRY = "sssom:NoTermFound"
WITHDRAWN_ENTRY_LABEL = "Withdrawn Entry"

# Config directory path

CONFIG_DIR = Path(__file__).parent.parent.parent.parent / "config"


@dataclass
class DatasourceConfig:
    """Configuration for a biological database datasource loaded from YAML."""

    name: str
    prefix: str
    curie_base_url: str
    default_output_filename: str = ""
    available_outputs: list[str] = field(default_factory=list)
    download_urls: dict[str, Any] = field(default_factory=dict)
    primary_file_key: str = ""
    id_pattern: str = ""
    archive_url: str = ""
    input_file_types: list[str] = field(default_factory=list)
    source: str = ""
    homepage: str = ""
    data_license: str = ""
    # SPARQL-based datasources (e.g., Wikidata)
    sparql_endpoint: str = ""
    queries: dict[str, str] = field(default_factory=dict)
    # For now, only ChEBI: version threshold for new TSV format.
    # Use if release files change location or serialization.
    new_format_version: int | None = None
    # Full metadata from YAML
    mappingset_metadata: dict[str, Any] = field(default_factory=dict)
    mapping_metadata: dict[str, Any] = field(default_factory=dict)


def load_config(datasource_name: str) -> dict[str, Any]:
    """Load configuration from a YAML file for a datasource.

    Args:
        datasource_name: Name of the datasource (e.g., 'chebi', 'hgnc').

    Returns:
        Dictionary with the full YAML configuration.

    Raises:
        FileNotFoundError: If the config file does not exist.
    """
    config_path = CONFIG_DIR / f"{datasource_name.lower()}.yaml"
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with config_path.open("r", encoding="utf-8") as f:
        result: dict[str, Any] = yaml.safe_load(f)
        return result


def get_datasource_config(datasource_name: str) -> DatasourceConfig:
    """Load and parse a DatasourceConfig from YAML.

    Args:
        datasource_name: Name of the datasource (e.g., 'chebi', 'hgnc').

    Returns:
        DatasourceConfig object populated from YAML.
    """
    raw = load_config(datasource_name)

    # Extract curie_map to get prefix and curie_base_url
    curie_map = raw.get("mappingset", {}).get("curie_map", {})
    # The first entry that's not a standard prefix is the datasource prefix
    standard_prefixes = {"IAO", "oboInOwl", "semapv", "skos", "sssom", "sec2pri"}
    prefix = ""
    curie_base_url = ""
    for k, v in curie_map.items():
        if k not in standard_prefixes:
            prefix = k
            curie_base_url = v
            break

    # Determine name - special case for wikidata
    if datasource_name.lower() == "wikidata":
        name = "Wikidata"
    else:
        name = datasource_name.upper()

    return DatasourceConfig(
        name=name,
        prefix=prefix,
        curie_base_url=curie_base_url,
        default_output_filename=raw.get("default_output_filename", ""),
        available_outputs=raw.get("available_outputs", []),
        download_urls=raw.get("download_urls", {}),
        primary_file_key=raw.get("primary_file_key", ""),
        id_pattern=raw.get("id_pattern", ""),
        archive_url=raw.get("archive_url", ""),
        input_file_types=raw.get("input_file_types", []),
        source=raw.get("source", ""),
        homepage=raw.get("homepage", ""),
        data_license=raw.get("data_license", ""),
        sparql_endpoint=raw.get("sparql_endpoint", ""),
        queries=raw.get("queries", {}),
        new_format_version=raw.get("new_format_version"),
        mappingset_metadata=raw.get("mappingset", {}),
        mapping_metadata=raw.get("mapping", {}),
    )


# Base Downloader Class


class BaseDownloader(ABC):
    """Abstract base class for datasource downloaders.

    Provides shared download logic that can be inherited by datasource-specific
    downloaders. Handles file downloads, URL construction, and version detection.
    """

    datasource_name: str = ""
    _config: DatasourceConfig | None = None

    def __init__(
        self,
        version: str | None = None,
        show_progress: bool = True,
    ) -> None:
        """Initialize the downloader.

        Args:
            version: Version/release identifier for the datasource.
            show_progress: Whether to show progress bars during downloads.
        """
        self.version = version
        self.show_progress = show_progress

        # Load config from YAML
        if self.datasource_name:
            try:
                self._config = get_datasource_config(self.datasource_name.lower())
            except FileNotFoundError:
                self._config = None

    @property
    def config(self) -> DatasourceConfig | None:
        """Get the loaded configuration."""
        return self._config

    @property
    def new_format_version(self) -> int | None:
        """Get the version threshold for new format (if any)."""
        if self._config:
            return self._config.new_format_version
        return None

    def is_new_format(self, version: str | None = None) -> bool:
        """Check if a version uses the new format.

        Args:
            version: Version to check. If None, uses self.version.

        Returns:
            True if version >= new_format_version threshold.
        """
        v = version or self.version
        threshold = self.new_format_version

        if threshold is None:
            return True  # No threshold means always "new" format

        if v is None:
            return True  # Default to new format for latest

        try:
            return int(v) >= threshold
        except ValueError:
            return True  # Default to new if version is not numeric

    @abstractmethod
    def get_download_urls(
        self,
        version: str | None = None,
        **kwargs: Any,
    ) -> dict[str, str]:
        """Get download URLs for the datasource.

        Args:
            version: Specific version to get URLs for.
            **kwargs: Additional options (e.g., subset, force_format).

        Returns:
            Dictionary mapping file keys to URLs.
        """

    @abstractmethod
    def download(
        self,
        output_dir: Path,
        version: str | None = None,
        decompress: bool = True,
        **kwargs: Any,
    ) -> dict[str, Path]:
        """Download files for the datasource.

        Args:
            output_dir: Directory to save downloaded files.
            version: Specific version to download.
            decompress: Whether to decompress .gz files.
            **kwargs: Additional options.

        Returns:
            Dictionary mapping file keys to downloaded paths.
        """

    def _download_file(
        self,
        url: str,
        output_path: Path,
        decompress_gz: bool = True,
        timeout: float | None = None,
        description: str | None = None,
    ) -> Path:
        """Download a file from URL to the specified path.

        Args:
            url: URL to download from.
            output_path: Where to save the file.
            decompress_gz: Whether to decompress .gz files automatically.
            timeout: Request timeout in seconds.
            description: Description for the progress bar.

        Returns:
            Path to the downloaded file.
        """
        # Import here to avoid circular imports
        from pysec2pri.download import download_file

        return download_file(
            url,
            output_path,
            decompress_gz=decompress_gz,
            timeout=timeout,
            show_progress=self.show_progress,
            description=description,
        )

    def _download_urls(
        self,
        urls: dict[str, str],
        output_dir: Path,
        decompress: bool = True,
    ) -> dict[str, Path]:
        """Download files from URLs to output directory.

        Args:
            urls: Dictionary mapping file keys to URLs.
            output_dir: Directory to save files.
            decompress: Whether to decompress .gz files.

        Returns:
            Dictionary mapping file keys to downloaded paths.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        downloaded: dict[str, Path] = {}

        for key, url in urls.items():
            filename = url.split("/")[-1]

            if decompress and filename.endswith(".gz"):
                filename = filename[:-3]

            output_path = output_dir / filename
            logger.info("Downloading %s: %s", key, url)
            self._download_file(url, output_path, decompress_gz=decompress)
            downloaded[key] = output_path
            logger.info("Saved to: %s", output_path)

        return downloaded


class Sec2PriMappingSet(MappingSet):  # type: ignore[misc]
    """A MappingSet for Sec2Pri, with helpers for cardinality."""

    def _compute_cardinalities(self, on: str = "id") -> None:
        """Compute and set mapping_cardinality on all mappings.

        'on' can be 'id' (uses subject_id/object_id) or 'label'.
        """
        if not self.mappings:  # type: ignore[has-type]
            return

        mappings = self._normalize_mappings()

        if on == "label":
            sec_field, pri_field = "subject_label", "object_label"
        else:
            sec_field, pri_field = "subject_id", "object_id"

        import polars as pl

        sec_vals = [str(getattr(m, sec_field, None) or "") for m in mappings]
        pri_vals = [str(getattr(m, pri_field, None) or "") for m in mappings]

        cardinalities: list[str] = (
            pl.DataFrame({"sec": sec_vals, "pri": pri_vals})
            .with_columns(
                pl.col("sec").count().over("sec").alias("s_count"),
                pl.col("pri").count().over("pri").alias("p_count"),
            )
            .select(
                pl.when(pl.col("s_count") == 1, pl.col("p_count") == 1)
                .then(pl.lit("1:1"))
                .when(pl.col("s_count") > 1, pl.col("p_count") == 1)
                .then(pl.lit("n:1"))
                .when(pl.col("s_count") == 1, pl.col("p_count") > 1)
                .then(pl.lit("1:n"))
                .when(pl.col("s_count") > 1, pl.col("p_count") > 1)
                .then(pl.lit("n:n"))
                .otherwise(pl.lit("1:0"))
                .alias("cardinality")
            )
            .get_column("cardinality")
            .to_list()
        )

        for m, card in zip(mappings, cardinalities, strict=False):
            m.mapping_cardinality = MappingCardinalityEnum(card)

        self.mappings = mappings

    def _normalize_mappings(self) -> list[Mapping]:
        """Normalize mappings to a list of Mapping objects.

        Returns:
            List of Mapping objects.
        """
        mappings = self.mappings
        if not isinstance(mappings, list):
            mappings = [mappings]
        for i, m in enumerate(mappings):
            if isinstance(m, dict):
                mappings[i] = Mapping(**m)
        return mappings


class IdMappingSet(Sec2PriMappingSet):
    """Child class for ID-based mapping sets."""

    def compute_cardinalities(self) -> None:
        """Compute cardinalities using subject_id and object_id fields."""
        self._compute_cardinalities(on="id")


class LabelMappingSet(Sec2PriMappingSet):
    """Child class for label-based mapping sets."""

    def compute_cardinalities(self) -> None:
        """Compute cardinalities using subject_label and object_label."""
        self._compute_cardinalities(on="label")


class BaseParser(ABC):
    """Abstract base class for all datasource parsers.

    Each parser is responsible for reading files from a specific datasource
    and extracting secondary-to-primary identifier Mapping Sets.
    """

    # To be overridden by subclasses
    datasource_name: str = ""
    default_source_url: str = ""
    _config: DatasourceConfig | None = None

    def __init__(
        self,
        version: str | None = None,
        show_progress: bool = True,
        config_name: str | None = None,
    ):
        """Initialize the parser.

        Args:
            version: Version/release identifier for the datasource.
            show_progress: Whether to show progress bars during parsing.
            config_name: Name of config file to load (defaults to class name).
        """
        self.version = version
        self.show_progress = show_progress

        # Load config from YAML if available
        if config_name:
            self._config = get_datasource_config(config_name)
        elif self.datasource_name:
            try:
                self._config = get_datasource_config(self.datasource_name.lower())
            except FileNotFoundError:
                self._config = None

    @property
    def config(self) -> DatasourceConfig | None:
        """Get the loaded configuration."""
        return self._config

    def get_download_url(self, key: str) -> str | None:
        """Get a download URL from config by key."""
        if self._config:
            return self._config.download_urls.get(key)
        return None

    def get_curie_map(self) -> dict[str, str]:
        """Get the CURIE map from config."""
        if self._config and self._config.mappingset_metadata:
            result: dict[str, str] = self._config.mappingset_metadata.get("curie_map", {})
            return result
        return {}

    def get_mappingset_metadata(self) -> dict[str, Any]:
        """Get mapping set metadata from config."""
        if self._config:
            result: dict[str, Any] = self._config.mappingset_metadata
            return result
        return {}

    def get_mapping_metadata(self) -> dict[str, Any]:
        """Get mapping metadata from config."""
        if self._config:
            result: dict[str, Any] = self._config.mapping_metadata
            return result
        return {}

    def load_metadata(self, yaml_path: str) -> dict[str, Any]:
        """Load metadata from a YAML config file."""
        with open(yaml_path, encoding="utf-8") as f:
            result: dict[str, Any] = yaml.safe_load(f)
            return result

    def apply_metadata_to_mappingset(
        self,
        mappingset: MappingSet,
        metadata: dict[str, Any],
    ) -> None:
        """Apply metadata to a MappingSet and its Mappings."""
        # Set MappingSet fields
        for key, value in metadata.get("mappingset", {}).items():
            if hasattr(mappingset, key) and value is not None:
                setattr(mappingset, key, value)
        # Set Mapping fields
        if hasattr(mappingset, "mappings") and mappingset.mappings:
            for mapping in mappingset.mappings:
                for key, value in metadata.get("mapping", {}).items():
                    if hasattr(mapping, key) and value is not None:
                        setattr(mapping, key, value)

    def _extract_version_from_file(self, file_path: Path) -> str | None:
        """Extract a version string embedded in a data file's header.

        Override in subclasses where the source file contains release
        metadata (e.g. ``Release: 2026_01`` in UniProt flat files).

        Args:
            file_path: Path to the data file to inspect.

        Returns:
            Version string, or ``None`` if not found.
        """
        return None

    def _resolve_version(self, file_path: Path | None = None) -> str:
        """Resolve the dataset version to use for source version fields.

        Resolution order:
        1. ``self.version`` if already set explicitly.
        2. Version extracted from file header via ``_extract_version_from_file``.
        3. ISO date or release token found in the filename stem
           (e.g. ``withdrawn_2026-04-07.txt`` -> ``2026-04-07``,
           ``chebi_245.sdf`` -> ``245``).
        4. File modification date (ISO-8601) when a path is provided.
        5. Today's date as a last resort.

        Sets ``self.version`` to the resolved value so that
        ``create_mapping_set`` picks it up for ``subject_source_version`` /
        ``object_source_version`` automatically.

        Args:
            file_path: Optional path to the primary input file.

        Returns:
            Resolved version string.
        """
        if self.version:
            return self.version

        if file_path is not None:
            file_path = Path(file_path)
            # 1. Try header-embedded version (parser-specific override)
            extracted = self._extract_version_from_file(file_path)
            if extracted:
                self.version = extracted
                return self.version
            # 2. Try ISO date (YYYY-MM-DD) in the filename stem
            iso_match = re.search(r"\d{4}-\d{2}-\d{2}", file_path.stem)
            if iso_match:
                self.version = iso_match.group(0)
                return self.version
            # 3. Try a plain numeric/semver token in the filename stem
            #    e.g. "chebi_245" -> "245", "gene_history_v2" -> "2"
            num_match = re.search(r"(?<![.\d])(\d{3,})(?![.\d])", file_path.stem)
            if num_match:
                self.version = num_match.group(1)
                return self.version
            # 4. Fall back to file modification time
            try:
                mtime = file_path.stat().st_mtime
                self.version = date.fromtimestamp(mtime).isoformat()
                return self.version
            except OSError:
                pass

        self.version = date.today().isoformat()
        return self.version

    @abstractmethod
    def parse(self, input_path: Path | str | None) -> MappingSet:
        """Parse the input file(s) and return a MappingSet.

        Args:
            input_path: Path to the input file or directory.

        Returns:
            A MappingSet containing all extracted mappings.
        """

    def _progress(
        self,
        iterable: Iterable[Any],
        desc: str | None = None,
        total: int | None = None,
    ) -> Iterable[Any]:
        """Wrap an iterable with a progress bar if enabled.

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

    def _build_mappings(
        self,
        rows: Iterable[dict[str, Any]],
        fixed_fields: dict[str, Any] | None = None,
        *,
        desc: str = "Building mappings",
        total: int | None = None,
    ) -> list[Mapping]:
        """Build SSSOM Mapping objects from row dicts.

        Automatically injects mapping-level fields from the parser's
        config metadata (e.g. ``confidence``) unless the caller already
        provides them in ``fixed_fields`` or individual row dicts.

        Args:
            rows: Per-row fields as dicts (subject_id, object_id, etc.).
            fixed_fields: Fields shared by all rows (predicate_id, license, etc.).
            desc: Progress bar description.
            total: Total count for the progress bar.

        Returns:
            List of Mapping objects.
        """
        _auto_fields = ("confidence",)
        m_meta = self.get_mapping_metadata()
        auto: dict[str, Any] = {
            k: m_meta[k] for k in _auto_fields if k in m_meta and m_meta[k] is not None
        }

        # Build base: auto fields < fixed_fields (fixed wins over auto)
        base = {**auto, **(fixed_fields or {})}

        if base:
            merged: Iterable[dict[str, Any]] = ({**base, **row} for row in rows)
        else:
            merged = rows
        return [Mapping(**row) for row in self._progress(merged, desc=desc, total=total)]

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
    def normalize_withdrawn_id(subject_id: str | None) -> str:
        """Normalize a primary ID, converting empty/null to withdrawn.

        Args:
            subject_id: The raw primary identifier from the source file.

        Returns:
            The normalized primary ID, or WITHDRAWN_ENTRY for empty values.
        """
        if not subject_id or subject_id in ("-", ""):
            return WITHDRAWN_ENTRY
        return subject_id

    @staticmethod
    def is_withdrawn(identifier: str) -> bool:
        """Return if is withdrawn."""
        return WITHDRAWN_ENTRY == identifier

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
        return subject_id == WITHDRAWN_ENTRY

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

    def create_mapping_set(
        self,
        mappings: list[Mapping],
        mapping_type: str = "id",
    ) -> Sec2PriMappingSet:
        """Create an IdMappingSet or LabelMappingSet with config metadata.

        This is the common factory method for creating mapping sets with
        all SSSOM metadata populated from the YAML config. It also computes
        cardinalities automatically using Polars vectorised counting.

        Args:
            mappings: List of SSSOM Mapping objects.
            mapping_type: "id" for IdMappingSet (cardinality by ID),
                         "label" for LabelMappingSet (cardinality by label).

        Returns:
            Configured MappingSet with computed cardinalities.
        """
        ms_meta = self.get_mappingset_metadata()
        curie_map = self.get_curie_map()

        # Choose the appropriate MappingSet class
        if mapping_type == "label":
            mapping_set_class = LabelMappingSet
        else:
            mapping_set_class = IdMappingSet

        # Build description with version if available
        description = ms_meta.get("mapping_set_description", "")
        if self.version and description:
            description = f"{description} Version: {self.version}."

        # Create the mapping set with all SSSOM metadata
        mapping_set = mapping_set_class(
            mappings=mappings,
            curie_map=curie_map,
            mapping_set_id=ms_meta.get("mapping_set_id"),
            mapping_set_version=self.version,
            mapping_set_title=ms_meta.get("mapping_set_title"),
            mapping_set_description=description or None,
            creator_id=ms_meta.get("creator_id"),
            creator_label=ms_meta.get("creator_label"),
            comment=ms_meta.get("comment"),
            license=ms_meta.get("license"),
            subject_source=ms_meta.get("subject_source"),
            subject_source_version=self.version,
            object_source=ms_meta.get("object_source"),
            object_source_version=self.version,
            mapping_provider=ms_meta.get("mapping_provider"),
            mapping_tool=ms_meta.get("mapping_tool"),
            mapping_tool_version=VERSION,
            mapping_date=date.today().isoformat(),
            see_also=ms_meta.get("see_also"),
            issue_tracker=ms_meta.get("issue_tracker"),
            subject_preprocessing=ms_meta.get("subject_preprocessing"),
            object_preprocessing=ms_meta.get("object_preprocessing"),
        )

        # Compute cardinalities using Polars vectorised counting
        mapping_set.compute_cardinalities()

        return mapping_set


__all__ = [
    "WITHDRAWN_ENTRY",
    "WITHDRAWN_ENTRY_LABEL",
    "BaseDownloader",
    "BaseParser",
    "DatasourceConfig",
    "IdMappingSet",
    "LabelMappingSet",
    "Sec2PriMappingSet",
    "get_datasource_config",
    "load_config",
]
