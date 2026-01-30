"""Pydantic data models for secondary to primary identifier mapping.

This module provides a simple, SSSOM-focused data model for representing
mappings between secondary (retired/withdrawn) identifiers and primary
(current) identifiers across various biological databases (ChEBI, HMDB,
HGNC, NCBI, UniProt, etc.).

SSSOM Semantics:
- subject_id = PRIMARY/CURRENT identifier (the one to use)
- object_id = SECONDARY/OLD identifier (the retired/merged one)
- subject_label = label/symbol for the primary
- object_label = label/symbol for the secondary

For legacy export formats:
- primaryID = subject_id
- secondaryID = object_id
- primarySymbol/primaryLabel = subject_label
- secondarySymbol/secondaryLabel = object_label
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from datetime import date
from enum import Enum
from typing import TYPE_CHECKING, ClassVar

from pydantic import BaseModel, ConfigDict, Field

from pysec2pri.constants import (
    CHEBI,
    HGNC,
    HMDB,
    MAPPING_JUSTIFICATION,
    NCBI,
    STANDARD_PREFIX_MAP,
    UNIPROT,
)

if TYPE_CHECKING:
    from sssom import MappingSetDocument
    from sssom.util import MappingSetDataFrame
    from sssom_schema import MappingSet as SSSOMMapingSet


# =============================================================================
# Enums for mapping metadata
# =============================================================================


class MappingCardinality(str, Enum):
    """Cardinality of mapping between secondary and primary identifiers."""

    ONE_TO_ONE = "1:1"  # Single secondary replaced by single primary
    MANY_TO_ONE = "n:1"  # Multiple secondaries merged into one primary
    ONE_TO_MANY = "1:n"  # Single secondary split into multiple primaries
    MANY_TO_MANY = "n:n"  # Complex merge/split
    ONE_TO_ZERO = "1:0"  # Withdrawn with no replacement


class PredicateID(str, Enum):
    """SSSOM predicate IDs."""

    TERM_REPLACED_BY = "IAO:0100001"  # Clear replacement (1:1, n:1)
    CONSIDER = "oboInOwl:consider"  # Split/complex/withdrawn (1:n, n:n, 1:0)
    SAMEAS = "owl:sameAs"  # Same entity


def get_predicate_for_cardinality(
    cardinality: MappingCardinality | None,
) -> PredicateID | None:
    """Get appropriate predicate for a cardinality type."""
    if cardinality is None:
        return None
    if cardinality in (
        MappingCardinality.ONE_TO_ONE,
        MappingCardinality.MANY_TO_ONE,
    ):
        return PredicateID.TERM_REPLACED_BY
    return PredicateID.CONSIDER


def get_comment_for_cardinality(cardinality: MappingCardinality | None) -> str:
    """Get human-readable comment for a cardinality type."""
    comments = {
        MappingCardinality.ONE_TO_ONE: "ID replaced.",
        MappingCardinality.MANY_TO_ONE: "IDs merged into one.",
        MappingCardinality.ONE_TO_MANY: "ID split into multiple.",
        MappingCardinality.MANY_TO_MANY: "IDs merged/split into multiple.",
        MappingCardinality.ONE_TO_ZERO: "ID withdrawn/deprecated.",
    }
    return comments.get(cardinality, "") if cardinality else ""


# =============================================================================
# Base Mapping Class
# =============================================================================


WITHDRAWN_ENTRY = "sec2pri:WithdrawnEntry"


class BaseMapping(BaseModel):
    """Base class for SSSOM-style subject-predicate-object mappings.

    Mapping semantics:
        subject_id = PRIMARY (current) identifier
        object_id = SECONDARY (old/retired) identifier
        subject_label = label for primary
        object_label = label for secondary

    For withdrawn entries with no replacement:
        subject_id = WITHDRAWN_ENTRY sentinel
        object_id = the withdrawn identifier
    """

    model_config = ConfigDict(use_enum_values=False, frozen=True)

    # Class-level constants (override in subclasses)
    datasource_name: ClassVar[str] = ""
    datasource_prefix: ClassVar[str] = ""
    curie_base_url: ClassVar[str] = ""

    # Instance fields: SymbolMapping and IdentifierMapping decide which fields
    # are needed
    subject_id: str = Field(..., description="Primary/current identifier")
    predicate_id: None | str = Field(..., description="SSSOM predicate")
    object_id: None | str = Field(..., description="Secondary/old identifier")
    subject_label: str | None = Field(default=None, description="Label for primary")
    object_label: str | None = Field(default=None, description="Label for secondary")
    comment: str | None = Field(default=None, description="Mapping context")
    source_url: str | None = Field(default=None, description="Source data URL")
    mapping_cardinality: MappingCardinality | None = Field(
        default=None, description="Cardinality of the mapping"
    )

    # Class constant for withdrawn entries
    withdrawn_entry: ClassVar[str] = WITHDRAWN_ENTRY

    @classmethod
    def curie(cls, identifier: str | None) -> str | None:
        """Convert identifier to CURIE format if needed."""
        if not identifier:
            return None
        if ":" in identifier:
            return identifier
        # Use getattr to get the value if it's a FieldInfo
        prefix = getattr(cls, "datasource_prefix", "")
        if isinstance(prefix, str) and prefix:
            return f"{prefix}:{identifier}"
        return identifier

    def to_sssom_dict(self) -> dict[str, str | None]:
        """Convert to SSSOM-compatible dictionary."""
        return {
            "subject_id": self.curie(self.subject_id),
            "subject_label": getattr(self, "subject_label", None),
            "predicate_id": self.predicate_id,
            "object_id": self.curie(self.object_id),
            "object_label": getattr(self, "object_label", None),
            "mapping_cardinality": (
                self.mapping_cardinality.name if self.mapping_cardinality else None
            ),
            "comment": self.comment,
            "source": self.source_url,
        }

    def to_legacy_dict(self) -> dict[str, str | None]:
        """Convert to legacy export format.

        Maps SSSOM fields to legacy column names:
        - subject_id -> primaryID
        - object_id -> secondaryID
        - subject_label -> primarySymbol/primaryLabel
        - object_label -> secondarySymbol/secondaryLabel
        """
        return {
            "primaryID": self.subject_id,
            "secondaryID": self.object_id,
            "primaryLabel": getattr(self, "subject_label", "") or "",
            "secondaryLabel": getattr(self, "object_label", "") or "",
            "primarySymbol": getattr(self, "subject_label", "") or "",
            "secondarySymbol": getattr(self, "object_label", "") or "",
            "predicateID": self.predicate_id,
            "mapping_cardinality_sec2pri": (
                self.mapping_cardinality.name if self.mapping_cardinality else ""
            ),
            "comment": self.comment or "",
            "source": self.source_url or "",
        }


# =============================================================================
# Datasource-specific mapping classes
# =============================================================================


class IdMapping(BaseMapping):
    """Mapping for ID-only resources (ChEBI, HMDB, UniProt).

    Requires subject_id, predicate_id, and object_id (all identifiers).
    Labels are optional.
    """

    # Override to make object_id required (not None)
    object_id: str | None = Field(default=None, description="Secondary/old ID")


class SymbolMapping(BaseMapping):
    """Mapping with symbol/label fields (HGNC, NCBI).

    Requires subject_id, predicate_id, and labels for subject and object.
    object_id may be None for symbol-only mappings.
    """

    # Override to make labels required
    subject_label: str | None = Field(default=None, description="Primary label/symbol")
    object_label: str | None = Field(default=None, description="Secondary label/symbol")


class ChEBIMapping(IdMapping):
    """ChEBI identifier mapping."""

    datasource_name: ClassVar[str] = CHEBI.name
    datasource_prefix: ClassVar[str] = CHEBI.prefix
    curie_base_url: ClassVar[str] = CHEBI.curie_base_url


class HMDBMapping(IdMapping):
    """HMDB identifier mapping."""

    datasource_name: ClassVar[str] = HMDB.name
    datasource_prefix: ClassVar[str] = HMDB.prefix
    curie_base_url: ClassVar[str] = HMDB.curie_base_url


class UniProtMapping(IdMapping):
    """UniProt identifier mapping."""

    datasource_name: ClassVar[str] = UNIPROT.name
    datasource_prefix: ClassVar[str] = UNIPROT.prefix
    curie_base_url: ClassVar[str] = UNIPROT.curie_base_url


class HGNCMapping(SymbolMapping):
    """HGNC gene identifier mapping."""

    datasource_name: ClassVar[str] = HGNC.name
    datasource_prefix: ClassVar[str] = HGNC.prefix
    curie_base_url: ClassVar[str] = HGNC.curie_base_url


class NCBIGeneMapping(SymbolMapping):
    """NCBI Gene (Entrez Gene) identifier mapping."""

    datasource_name: ClassVar[str] = NCBI.name
    datasource_prefix: ClassVar[str] = NCBI.prefix
    curie_base_url: ClassVar[str] = NCBI.curie_base_url


# =============================================================================
# Mapping set container
# =============================================================================


class MappingSet(BaseModel):
    """A collection of mappings with associated metadata.

    This class represents a complete set of mappings from a datasource,
    including metadata needed for SSSOM output.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    mappings: list[BaseMapping] = Field(default_factory=list)
    datasource_name: str = Field(..., description="Name of the datasource")
    version: str | None = Field(
        default=None,
        description="Version/release of the datasource",
    )
    mapping_set_id: str | None = Field(
        default=None,
        description="Unique identifier for this set",
    )
    mapping_set_description: str | None = Field(
        default=None,
        description="Description of the mapping set",
    )
    mapping_date: str | None = Field(
        default=None,
        description="Date the mappings were generated",
    )
    license_url: str = Field(
        default="https://creativecommons.org/publicdomain/zero/1.0/",
        description="License URL for the mapping set",
    )
    curie_map: dict[str, str] = Field(
        default_factory=dict,
        description="CURIE prefix to URL mappings",
    )
    source_metadata: str | None = Field(
        default=None, description="Source metadata (version info, etc.)"
    )
    comment: str | None = Field(
        default=None,
        description="General comment about the mapping set",
    )

    def add_mapping(self, mapping: BaseMapping) -> None:
        """Add a mapping to the set."""
        # Convert to list if needed
        if not isinstance(self.mappings, list):
            self.mappings = list(self.mappings)
        self.mappings.append(mapping)

    def add_mappings(self, mappings: Sequence[BaseMapping]) -> None:
        """Add multiple mappings to the set."""
        if not isinstance(self.mappings, list):
            self.mappings = list(self.mappings)
        self.mappings.extend(mappings)

    def __len__(self) -> int:
        """Return the number of mappings in the set."""
        return len(self.mappings)

    def iter_mappings(self) -> Iterator[BaseMapping]:
        """Iterate over mappings."""
        return iter(self.mappings)

    def to_sssom_mapping_set(
        self,
        mapping_date: str | None = None,
    ) -> SSSOMMapingSet:
        """Convert this MappingSet to an sssom_schema MappingSet.

        Args:
            mapping_date: Date string (YYYY-MM-DD). Defaults to today.

        Returns:
            An sssom_schema MappingSet object.
        """
        from sssom_schema import Mapping as SSSOMMapping
        from sssom_schema import MappingSet as SSSOMMapingSet

        if mapping_date is None:
            mapping_date = date.today().isoformat()

        # Convert individual mappings
        sssom_mappings = []
        for m in self.mappings:
            sssom_dict = m.to_sssom_dict()

            # Build SSSOM Mapping with required fields
            predicate_id = sssom_dict.get("predicate_id")
            if predicate_id is None:
                predicate_id = "skos:relatedMatch"  # fallback

            sssom_mapping = SSSOMMapping(
                subject_id=sssom_dict.get("subject_id"),
                subject_label=sssom_dict.get("subject_label"),
                predicate_id=predicate_id,
                object_id=sssom_dict.get("object_id"),
                object_label=sssom_dict.get("object_label"),
                mapping_justification=MAPPING_JUSTIFICATION,
                mapping_cardinality=sssom_dict.get("mapping_cardinality"),
                comment=sssom_dict.get("comment"),
            )
            sssom_mappings.append(sssom_mapping)

        # Build the MappingSet
        sssom_mapping_set = SSSOMMapingSet(
            mapping_set_id=(self.mapping_set_id or "https://w3id.org/sssom/mappings"),
            license=self.license_url,
            mappings=sssom_mappings,
            mapping_set_version=self.version,
            mapping_set_description=self.mapping_set_description,
            mapping_date=mapping_date,
            comment=self.comment,
        )

        return sssom_mapping_set

    def to_sssom_document(
        self,
        mapping_date: str | None = None,
    ) -> MappingSetDocument:
        """Convert this MappingSet to an SSSOM MappingSetDocument.

        Args:
            mapping_date: Date string (YYYY-MM-DD). Defaults to today.

        Returns:
            An SSSOM MappingSetDocument with converter.
        """
        from curies import Converter
        from sssom import MappingSetDocument

        sssom_ms = self.to_sssom_mapping_set(mapping_date)

        # Build prefix map for converter (start with standard prefixes)
        prefix_map = dict(STANDARD_PREFIX_MAP)

        # Add datasource-specific prefixes
        for prefix, url in dict(self.curie_map).items():
            # Remove trailing colon if present
            clean_prefix = prefix.rstrip(":")
            prefix_map[clean_prefix] = url

        converter = Converter.from_prefix_map(prefix_map)

        return MappingSetDocument(
            mapping_set=sssom_ms,
            converter=converter,
        )

    def to_sssom_dataframe(
        self,
        mapping_date: str | None = None,
    ) -> MappingSetDataFrame:
        """Convert this MappingSet to an SSSOM MappingSetDataFrame.

        Args:
            mapping_date: Date string (YYYY-MM-DD). Defaults to today.

        Returns:
            An SSSOM MappingSetDataFrame.

        Raises:
            ImportError: If pandas is not installed.
        """
        from sssom.parsers import to_mapping_set_dataframe  # type: ignore[attr-defined]

        doc = self.to_sssom_document(mapping_date)
        return to_mapping_set_dataframe(doc)


# =============================================================================
# Utility functions for computing cardinality
# =============================================================================


def compute_cardinality(
    mappings: list[tuple[str, str]],
) -> dict[tuple[str, str], MappingCardinality]:
    """Compute cardinality for a list of (subject_id, object_id) pairs.

    This function analyzes the relationships between primary and secondary IDs
    to determine the cardinality of each mapping.

    Args:
        mappings: List of (subject_id, object_id) tuples.

    Returns:
    Dictionary mapping each (subject_id, object_id) pair to its cardinality.
    """
    from collections import Counter

    # Count how many times each primary/secondary appears
    primary_counts: Counter[str] = Counter()
    secondary_counts: Counter[str] = Counter()

    for subject_id, object_id in mappings:
        if subject_id and object_id:
            primary_counts[subject_id] += 1
            secondary_counts[object_id] += 1

    result: dict[tuple[str, str], MappingCardinality] = {}

    for subject_id, object_id in mappings:
        if not object_id:
            continue

        # Check for withdrawn entries
        if subject_id in (BaseMapping.withdrawn_entry, "-", "", None):
            result[(subject_id, object_id)] = MappingCardinality.ONE_TO_ZERO
            continue

        pri_count = primary_counts.get(subject_id, 1)
        sec_count = secondary_counts.get(object_id, 1)

        if pri_count > 1 and sec_count == 1:
            # Multiple secondary IDs map to one primary ID (merge)
            cardinality = MappingCardinality.MANY_TO_ONE
        elif pri_count == 1 and sec_count == 1:
            # Simple 1:1 replacement
            cardinality = MappingCardinality.ONE_TO_ONE
        elif pri_count == 1 and sec_count > 1:
            # One secondary ID maps to multiple primary IDs (split)
            cardinality = MappingCardinality.ONE_TO_MANY
        elif pri_count > 1 and sec_count > 1:
            # Complex many-to-many relationship
            cardinality = MappingCardinality.MANY_TO_MANY
        else:
            cardinality = MappingCardinality.ONE_TO_ONE  # fallback

        result[(subject_id, object_id)] = cardinality

    return result


__all__ = [
    # Constants
    "WITHDRAWN_ENTRY",
    # Base classes
    "BaseMapping",
    # Datasource-specific classes
    "ChEBIMapping",
    "HGNCMapping",
    "HMDBMapping",
    "IdMapping",
    # Enums
    "MappingCardinality",
    # Container
    "MappingSet",
    "NCBIGeneMapping",
    "PredicateID",
    "SymbolMapping",
    "UniProtMapping",
    # Utility functions
    "compute_cardinality",
    "get_comment_for_cardinality",
    "get_predicate_for_cardinality",
]
