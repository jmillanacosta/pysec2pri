"""Pydantic data models for secondary to primary identifier mapping.

This module provides a class-based, DRY data model for representing mappings
between secondary (retired/withdrawn) identifiers and primary (current) identifiers
across various biological databases (ChEBI, HMDB, HGNC, NCBI, UniProt, etc.).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from enum import Enum
from typing import ClassVar, Optional

from pydantic import BaseModel, ConfigDict, Field, computed_field

from pysec2pri.constants import (
    CHEBI,
    HGNC,
    HMDB,
    NCBI,
    UNIPROT,
    WITHDRAWN_CURIE,
)


class MappingCardinality(str, Enum):
    """Cardinality of mapping between secondary and primary identifiers.

    - ONE_TO_ONE (1:1): A single secondary ID replaced by a single primary ID.
    - MANY_TO_ONE (n:1): Multiple secondary IDs merged into one primary ID.
    - ONE_TO_MANY (1:n): A single secondary ID split into multiple primary IDs.
    - MANY_TO_MANY (n:n): Multiple secondary IDs merged/split into multiple primary IDs.
    - ONE_TO_ZERO (1:0): A secondary ID that was withdrawn/deprecated with no replacement.
    """

    ONE_TO_ONE = "1:1"
    MANY_TO_ONE = "n:1"
    ONE_TO_MANY = "1:n"
    MANY_TO_MANY = "n:n"
    ONE_TO_ZERO = "1:0"


class PredicateID(str, Enum):
    """Predicate IDs for SSSOM mapping relationships.

    - TERM_REPLACED_BY: The subject term was replaced by the object term (IAO:0100001).
      Used for 1:1 and n:1 mappings (clear replacement or merge).
    - CONSIDER: The subject term should be considered in relation to the object
      (oboInOwl:consider). Used for 1:n, n:n, and 1:0 mappings (split, complex, withdrawn).
    """

    TERM_REPLACED_BY = "IAO:0100001"
    CONSIDER = "oboInOwl:consider"


def get_predicate_for_cardinality(cardinality: MappingCardinality | None) -> PredicateID | None:
    """Determine the appropriate predicate based on mapping cardinality.

    Args:
        cardinality: The mapping cardinality type.

    Returns:
        The appropriate predicate ID or None if cardinality is None.
    """
    if cardinality is None:
        return None
    if cardinality in (MappingCardinality.ONE_TO_ONE, MappingCardinality.MANY_TO_ONE):
        return PredicateID.TERM_REPLACED_BY
    return PredicateID.CONSIDER


def get_comment_for_cardinality(cardinality: MappingCardinality | None) -> str:
    """Generate a human-readable comment for a mapping cardinality.

    Args:
        cardinality: The mapping cardinality type.

    Returns:
        A descriptive comment string.
    """
    comments = {
        MappingCardinality.ONE_TO_ONE: "ID (subject) is replaced.",
        MappingCardinality.MANY_TO_ONE: "This ID (subject) and other ID(s) are merged into one ID.",
        MappingCardinality.ONE_TO_MANY: "ID (subject) is split into multiple.",
        MappingCardinality.MANY_TO_MANY: (
            "This ID (subject) and other ID(s) are merged/split into multiple ID(s)."
        ),
        MappingCardinality.ONE_TO_ZERO: "ID (subject) withdrawn/deprecated.",
    }
    return comments.get(cardinality, "") if cardinality else ""


class BaseMapping(BaseModel, ABC):
    """Abstract base class for all mapping types.

    Provides common functionality for secondary-to-primary identifier mappings.
    Subclasses implement specific mapping types (ID-only or ID+symbol).
    """

    model_config = ConfigDict(use_enum_values=False, frozen=True)

    # Class-level metadata to be overridden by subclasses
    datasource_name: ClassVar[str] = ""
    datasource_prefix: ClassVar[str] = ""
    curie_base_url: ClassVar[str] = ""

    primary_id: str = Field(..., description="The primary (current) identifier")
    secondary_id: str | None = Field(default=None, description="The secondary (retired) identifier")
    mapping_cardinality: MappingCardinality | None = Field(
        default=None, description="Cardinality of the mapping relationship"
    )
    comment: str | None = Field(default=None, description="Additional context about the mapping")
    source_url: str | None = Field(default=None, description="URL of the source data file")

    @computed_field
    @property
    def predicate_id(self) -> PredicateID | None:
        """Compute the predicate based on cardinality."""
        return get_predicate_for_cardinality(self.mapping_cardinality)

    @computed_field
    @property
    def is_withdrawn(self) -> bool:
        """Check if this mapping represents a withdrawn/deprecated entry."""
        return self.mapping_cardinality == MappingCardinality.ONE_TO_ZERO

    # CURIE for withdrawn/deleted entries with no replacement
    # Uses sec2pri namespace: sec2pri:WithdrawnEntry
    WITHDRAWN_SENTINEL: ClassVar[str] = WITHDRAWN_CURIE

    @classmethod
    def curie(cls, identifier: str | None) -> str | None:
        """Format an identifier as a CURIE (Compact URI).

        Args:
            identifier: The raw identifier.

        Returns:
            The identifier formatted with the datasource prefix,
            or None if identifier is None/empty.
        """
        if not identifier:
            return None
        # Withdrawn sentinel is already a valid CURIE
        if identifier == cls.WITHDRAWN_SENTINEL:
            return identifier
        # If already a CURIE (contains :), return as-is
        if ":" in identifier:
            return identifier
        if cls.datasource_prefix:
            return f"{cls.datasource_prefix}:{identifier}"
        return identifier

    @classmethod
    def format_primary_id(cls, primary_id: str) -> str | None:
        """Format primary ID for output, handling withdrawn entries.

        Args:
            primary_id: The primary identifier.

        Returns:
            Formatted ID or None for withdrawn entries.
        """
        if primary_id == cls.WITHDRAWN_SENTINEL:
            return None
        return cls.curie(primary_id)

    @abstractmethod
    def to_sssom_dict(self) -> dict:
        """Convert the mapping to a dictionary suitable for SSSOM output.

        Returns:
            Dictionary with SSSOM-compatible field names.
        """

    def to_legacy_dict(self) -> dict:
        """Convert the mapping to a dictionary for legacy output formats.

        Returns a consistent dictionary with all possible fields.
        Subclasses may override to add additional fields.

        Returns:
            Dictionary with legacy format field names.
        """
        return {
            "primaryID": self.primary_id,
            "secondaryID": self.secondary_id or "",
            "predicateID": self.predicate_id.value if self.predicate_id else "",
            "mapping_cardinality_sec2pri": (
                self.mapping_cardinality.value if self.mapping_cardinality else ""
            ),
            "comment": self.comment or "",
            "source": self.source_url or "",
            "objectID": self.primary_id,
        }


class IdMapping(BaseMapping):
    """Mapping between secondary and primary identifiers only.

    Used for databases like ChEBI, HMDB, and UniProt where the mapping
    is purely identifier-based without symbol/name information.
    """

    primary_label: str | None = Field(
        default=None, description="Name/label for the primary ID"
    )
    secondary_label: str | None = Field(
        default=None, description="Synonym/label for the secondary ID"
    )

    def to_sssom_dict(self) -> dict:
        """Convert to SSSOM dictionary format."""
        predicate = self.predicate_id.value if self.predicate_id else None
        object_id = self.curie(self.secondary_id) if self.secondary_id else None
        cardinality = (
            self.mapping_cardinality.value if self.mapping_cardinality else None
        )
        return {
            "subject_id": self.curie(self.primary_id),
            "subject_label": self.primary_label,
            "predicate_id": predicate,
            "object_id": object_id,
            "object_label": self.secondary_label,
            "mapping_cardinality": cardinality,
            "comment": self.comment,
            "source": self.source_url,
        }

    def to_legacy_dict(self) -> dict:
        """Convert to legacy format with label fields."""
        base = super().to_legacy_dict()
        base["primaryLabel"] = self.primary_label or ""
        base["secondaryLabel"] = self.secondary_label or ""
        return base


class SymbolMapping(BaseMapping):
    """Mapping that includes symbol/gene name information.

    Used for gene databases like HGNC and NCBI where both identifiers
    and symbols need to be tracked.
    """

    primary_symbol: str | None = Field(
        default=None, description="Primary symbol (e.g., gene symbol)"
    )
    secondary_symbol: str | None = Field(
        default=None, description="Secondary/previous symbol"
    )

    def to_sssom_dict(self) -> dict:
        """Convert to SSSOM dictionary format."""
        predicate = self.predicate_id.value if self.predicate_id else None
        object_id = self.curie(self.secondary_id) if self.secondary_id else None
        cardinality = (
            self.mapping_cardinality.value if self.mapping_cardinality else None
        )
        return {
            "subject_id": self.curie(self.primary_id),
            "subject_label": self.primary_symbol,
            "predicate_id": predicate,
            "object_id": object_id,
            "object_label": self.secondary_symbol,
            "mapping_cardinality": cardinality,
            "comment": self.comment,
            "source": self.source_url,
        }

    def to_legacy_dict(self) -> dict:
        """Convert to legacy format with symbol fields."""
        base = super().to_legacy_dict()
        base["primarySymbol"] = self.primary_symbol or ""
        base["secondarySymbol"] = self.secondary_symbol or ""
        return base


# =============================================================================
# Datasource-specific mapping classes
# =============================================================================


class ChEBIMapping(IdMapping):
    """ChEBI (Chemical Entities of Biological Interest) identifier mapping."""

    datasource_name: ClassVar[str] = CHEBI.name
    datasource_prefix: ClassVar[str] = CHEBI.prefix
    curie_base_url: ClassVar[str] = CHEBI.curie_base_url


class HMDBMapping(IdMapping):
    """HMDB (Human Metabolome Database) identifier mapping."""

    datasource_name: ClassVar[str] = HMDB.name
    datasource_prefix: ClassVar[str] = HMDB.prefix
    curie_base_url: ClassVar[str] = HMDB.curie_base_url


class UniProtMapping(IdMapping):
    """UniProt protein identifier mapping."""

    datasource_name: ClassVar[str] = UNIPROT.name
    datasource_prefix: ClassVar[str] = UNIPROT.prefix
    curie_base_url: ClassVar[str] = UNIPROT.curie_base_url


class HGNCMapping(SymbolMapping):
    """HGNC (HUGO Gene Nomenclature Committee) gene identifier mapping."""

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
    version: str | None = Field(default=None, description="Version/release of the datasource")
    mapping_set_id: str | None = Field(default=None, description="Unique identifier for this set")
    mapping_set_description: str | None = Field(
        default=None, description="Description of the mapping set"
    )
    mapping_date: str | None = Field(default=None, description="Date the mappings were generated")
    license_url: str = Field(
        default="https://creativecommons.org/publicdomain/zero/1.0/",
        description="License URL for the mapping set",
    )
    curie_map: dict[str, str] = Field(
        default_factory=dict, description="CURIE prefix to URL mappings"
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
        self.mappings.append(mapping)

    def add_mappings(self, mappings: list[BaseMapping]) -> None:
        """Add multiple mappings to the set."""
        self.mappings.extend(mappings)

    def __len__(self) -> int:
        """Return the number of mappings in the set."""
        return len(self.mappings)

    def iter_mappings(self):
        """Iterate over mappings."""
        return iter(self.mappings)


# =============================================================================
# Utility functions for computing cardinality
# =============================================================================


def compute_cardinality(
    mappings: list[tuple[str, str]],
) -> dict[tuple[str, str], MappingCardinality]:
    """Compute cardinality for a list of (primary_id, secondary_id) pairs.

    This function analyzes the relationships between primary and secondary IDs
    to determine the cardinality of each mapping.

    Args:
        mappings: List of (primary_id, secondary_id) tuples.

    Returns:
        Dictionary mapping each (primary_id, secondary_id) pair to its cardinality.
    """
    from collections import Counter

    # Count how many times each primary/secondary appears
    primary_counts: Counter[str] = Counter()
    secondary_counts: Counter[str] = Counter()

    for primary_id, secondary_id in mappings:
        if primary_id and secondary_id:
            primary_counts[primary_id] += 1
            secondary_counts[secondary_id] += 1

    result: dict[tuple[str, str], MappingCardinality] = {}

    for primary_id, secondary_id in mappings:
        if not secondary_id:
            continue

        # Check for withdrawn entries
        if primary_id in (BaseMapping.WITHDRAWN_SENTINEL, "-", "", None):
            result[(primary_id, secondary_id)] = MappingCardinality.ONE_TO_ZERO
            continue

        pri_count = primary_counts.get(primary_id, 1)
        sec_count = secondary_counts.get(secondary_id, 1)

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

        result[(primary_id, secondary_id)] = cardinality

    return result


__all__ = [
    # Enums
    "MappingCardinality",
    "PredicateID",
    # Base classes
    "BaseMapping",
    "IdMapping",
    "SymbolMapping",
    # Datasource-specific classes
    "ChEBIMapping",
    "HMDBMapping",
    "UniProtMapping",
    "HGNCMapping",
    "NCBIGeneMapping",
    # Container
    "MappingSet",
    # Utility functions
    "compute_cardinality",
    "get_predicate_for_cardinality",
    "get_comment_for_cardinality",
]
