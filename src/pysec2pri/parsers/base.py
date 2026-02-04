"""Base parser class for all datasource parsers."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections import Counter
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any, cast

import yaml
from sssom_schema import Mapping, MappingCardinalityEnum, MappingSet
from tqdm import tqdm

from pysec2pri.version import VERSION

# =============================================================================
# Values for withdrawn entries
# =============================================================================

WITHDRAWN_ENTRY = "sssom:NoTermFound"
WITHDRAWN_ENTRY_LABEL = "Withdrawn Entry"

# =============================================================================
# Config directory path
# =============================================================================

CONFIG_DIR = Path(__file__).parent.parent.parent.parent / "config"


@dataclass
class DatasourceConfig:
    """Configuration for a biological database datasource loaded from YAML."""

    name: str
    prefix: str
    curie_base_url: str
    default_output_filename: str = ""
    available_outputs: list[str] = field(default_factory=list)
    download_urls: dict[str, str] = field(default_factory=dict)
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
        mappingset_metadata=raw.get("mappingset", {}),
        mapping_metadata=raw.get("mapping", {}),
    )


def _determine_cardinality(p_count: int, s_count: int) -> str:
    """Determine SSSOM cardinality value from counts.

    Args:
        p_count: Number of mappings with the same primary.
        s_count: Number of mappings with the same secondary.

    Returns:
        SSSOM cardinality string.
    """
    if not p_count or not s_count:
        return "1:0"
    if s_count == 1 and p_count == 1:
        return "1:1"
    if s_count > 1 and p_count == 1:
        return "n:1"
    if s_count == 1 and p_count > 1:
        return "1:n"
    return "n:n"


def _get_field_accessors(
    on: str,
) -> tuple[Callable[[Mapping], str], Callable[[Mapping], str]]:
    """Get accessor functions for primary and secondary fields.

    Args:
        on: Either 'id' or 'label'.

    Returns:
        Tuple of (get_primary, get_secondary) functions.
    """
    if on == "label":
        return (
            lambda m: str(getattr(m, "subject_label", None)),
            lambda m: str(getattr(m, "object_label", None)),
        )
    return (
        lambda m: str(getattr(m, "subject_id", None)),
        lambda m: str(getattr(m, "object_id", None)),
    )


class Sec2PriMappingSet(MappingSet):  # type: ignore[misc]
    """A MappingSet for Sec2Pri, with helpers for cardinality."""

    def _compute_cardinalities(self, on: str = "id") -> None:
        """Compute mapping cardinalities and update each Mapping.

        Updates each Mapping's mapping_cardinality field using SSSOM-compliant
        values. 'on' can be 'id' or 'label'.
        """
        if not self.mappings:  # type: ignore[has-type]
            return

        mappings = self._normalize_mappings()
        get_primary, get_secondary = _get_field_accessors(on)

        primary_counts = Counter(get_primary(m) for m in mappings)
        secondary_counts = Counter(get_secondary(m) for m in mappings)

        for m in mappings:
            p_count = primary_counts.get(get_primary(m), 0)
            s_count = secondary_counts.get(get_secondary(m), 0)
            cardinality = _determine_cardinality(p_count, s_count)
            # Use constructor to get proper enum type
            m.mapping_cardinality = MappingCardinalityEnum(cardinality)

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
        cardinalities automatically.

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

        # Compute cardinalities using the appropriate method
        mapping_set.compute_cardinalities()

        return mapping_set


__all__ = [
    "WITHDRAWN_ENTRY",
    "WITHDRAWN_ENTRY_LABEL",
    "BaseParser",
    "DatasourceConfig",
    "IdMappingSet",
    "LabelMappingSet",
    "Sec2PriMappingSet",
    "get_datasource_config",
    "load_config",
]
