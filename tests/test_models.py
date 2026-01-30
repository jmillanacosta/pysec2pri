"""Tests for pysec2pri.models module.

Tests cover:
- Enums (MappingCardinality, PredicateID)
- Helper functions (get_predicate_for_cardinality, get_comment_for_cardinality)
- BaseMapping class and its methods
- Datasource-specific mapping classes
- MappingSet container class
- compute_cardinality function
"""

from pysec2pri.models import (
    WITHDRAWN_ENTRY,
    BaseMapping,
    ChEBIMapping,
    HGNCMapping,
    HMDBMapping,
    IdMapping,
    MappingCardinality,
    MappingSet,
    NCBIGeneMapping,
    PredicateID,
    SymbolMapping,
    UniProtMapping,
    compute_cardinality,
    get_comment_for_cardinality,
    get_predicate_for_cardinality,
)

# =============================================================================
# Tests for Enums
# =============================================================================


class TestMappingCardinality:
    """Tests for MappingCardinality enum."""

    def test_cardinality_values(self) -> None:
        """Test that cardinality enum has expected values."""
        assert MappingCardinality.ONE_TO_ONE.value == "1:1"
        assert MappingCardinality.MANY_TO_ONE.value == "n:1"
        assert MappingCardinality.ONE_TO_MANY.value == "1:n"
        assert MappingCardinality.MANY_TO_MANY.value == "n:n"
        assert MappingCardinality.ONE_TO_ZERO.value == "1:0"

    def test_cardinality_is_string_enum(self) -> None:
        """Test that cardinality can be used as string."""
        assert str(MappingCardinality.ONE_TO_ONE) == "MappingCardinality.ONE_TO_ONE"
        assert MappingCardinality.ONE_TO_ONE.value == "1:1"


class TestPredicateID:
    """Tests for PredicateID enum."""

    def test_predicate_values(self) -> None:
        """Test that predicate enum has expected values."""
        assert PredicateID.TERM_REPLACED_BY.value == "IAO:0100001"
        assert PredicateID.CONSIDER.value == "oboInOwl:consider"
        assert PredicateID.SAMEAS.value == "owl:sameAs"


# =============================================================================
# Tests for Helper Functions
# =============================================================================


class TestGetPredicateForCardinality:
    """Tests for get_predicate_for_cardinality function."""

    def test_one_to_one_returns_replaced_by(self) -> None:
        """1:1 cardinality should return TERM_REPLACED_BY."""
        result = get_predicate_for_cardinality(MappingCardinality.ONE_TO_ONE)
        assert result == PredicateID.TERM_REPLACED_BY

    def test_many_to_one_returns_replaced_by(self) -> None:
        """n:1 cardinality should return TERM_REPLACED_BY."""
        result = get_predicate_for_cardinality(MappingCardinality.MANY_TO_ONE)
        assert result == PredicateID.TERM_REPLACED_BY

    def test_one_to_many_returns_consider(self) -> None:
        """1:n cardinality should return CONSIDER."""
        result = get_predicate_for_cardinality(MappingCardinality.ONE_TO_MANY)
        assert result == PredicateID.CONSIDER

    def test_many_to_many_returns_consider(self) -> None:
        """n:n cardinality should return CONSIDER."""
        result = get_predicate_for_cardinality(MappingCardinality.MANY_TO_MANY)
        assert result == PredicateID.CONSIDER

    def test_one_to_zero_returns_consider(self) -> None:
        """1:0 cardinality should return CONSIDER."""
        result = get_predicate_for_cardinality(MappingCardinality.ONE_TO_ZERO)
        assert result == PredicateID.CONSIDER

    def test_none_returns_none(self) -> None:
        """None cardinality should return None."""
        result = get_predicate_for_cardinality(None)
        assert result is None


class TestGetCommentForCardinality:
    """Tests for get_comment_for_cardinality function."""

    def test_one_to_one_comment(self) -> None:
        """1:1 cardinality should have appropriate comment."""
        result = get_comment_for_cardinality(MappingCardinality.ONE_TO_ONE)
        assert result == "ID replaced."

    def test_many_to_one_comment(self) -> None:
        """n:1 cardinality should have appropriate comment."""
        result = get_comment_for_cardinality(MappingCardinality.MANY_TO_ONE)
        assert result == "IDs merged into one."

    def test_one_to_many_comment(self) -> None:
        """1:n cardinality should have appropriate comment."""
        result = get_comment_for_cardinality(MappingCardinality.ONE_TO_MANY)
        assert result == "ID split into multiple."

    def test_many_to_many_comment(self) -> None:
        """n:n cardinality should have appropriate comment."""
        result = get_comment_for_cardinality(MappingCardinality.MANY_TO_MANY)
        assert result == "IDs merged/split into multiple."

    def test_one_to_zero_comment(self) -> None:
        """1:0 cardinality should have appropriate comment."""
        result = get_comment_for_cardinality(MappingCardinality.ONE_TO_ZERO)
        assert result == "ID withdrawn/deprecated."

    def test_none_returns_empty_string(self) -> None:
        """None cardinality should return empty string."""
        result = get_comment_for_cardinality(None)
        assert result == ""


# =============================================================================
# Tests for BaseMapping
# =============================================================================


class TestBaseMapping:
    """Tests for BaseMapping class."""

    def test_create_basic_mapping(self) -> None:
        """Test creating a basic mapping with required fields."""
        mapping = BaseMapping(
            subject_id="PRIMARY123",
            predicate_id="IAO:0100001",
            object_id="SECONDARY456",
        )
        assert mapping.subject_id == "PRIMARY123"
        assert mapping.predicate_id == "IAO:0100001"
        assert mapping.object_id == "SECONDARY456"
        assert mapping.subject_label is None
        assert mapping.object_label is None

    def test_create_mapping_with_all_fields(self) -> None:
        """Test creating a mapping with all fields."""
        mapping = BaseMapping(
            subject_id="PRIMARY123",
            predicate_id="IAO:0100001",
            object_id="SECONDARY456",
            subject_label="Primary Label",
            object_label="Secondary Label",
            comment="Test comment",
            source_url="https://example.com",
            mapping_cardinality=MappingCardinality.ONE_TO_ONE,
        )
        assert mapping.subject_label == "Primary Label"
        assert mapping.object_label == "Secondary Label"
        assert mapping.comment == "Test comment"
        assert mapping.source_url == "https://example.com"
        assert mapping.mapping_cardinality == MappingCardinality.ONE_TO_ONE

    def test_withdrawn_entry_constant(self) -> None:
        """Test the withdrawn entry constant."""
        assert BaseMapping.withdrawn_entry == WITHDRAWN_ENTRY
        assert WITHDRAWN_ENTRY == "sec2pri:WithdrawnEntry"

    def test_curie_with_prefix(self) -> None:
        """Test curie method when identifier has no prefix."""
        # BaseMapping has no datasource_prefix set
        result = BaseMapping.curie("12345")
        assert result == "12345"  # No prefix added

    def test_curie_with_existing_curie(self) -> None:
        """Test curie method when identifier is already a CURIE."""
        result = BaseMapping.curie("CHEBI:12345")
        assert result == "CHEBI:12345"

    def test_curie_with_none(self) -> None:
        """Test curie method with None input."""
        result = BaseMapping.curie(None)
        assert result is None

    def test_curie_with_empty_string(self) -> None:
        """Test curie method with empty string."""
        result = BaseMapping.curie("")
        assert result is None


class TestBaseMappingToSssomDict:
    """Tests for BaseMapping.to_sssom_dict method."""

    def test_to_sssom_dict_basic(self) -> None:
        """Test converting mapping to SSSOM dictionary."""
        mapping = BaseMapping(
            subject_id="PRIMARY123",
            predicate_id="IAO:0100001",
            object_id="SECONDARY456",
        )
        result = mapping.to_sssom_dict()
        assert result["subject_id"] == "PRIMARY123"
        assert result["predicate_id"] == "IAO:0100001"
        assert result["object_id"] == "SECONDARY456"
        assert result["subject_label"] is None
        assert result["object_label"] is None
        assert result["mapping_cardinality"] is None

    def test_to_sssom_dict_with_cardinality(self) -> None:
        """Test that cardinality value is extracted correctly."""
        mapping = BaseMapping(
            subject_id="PRIMARY123",
            predicate_id="IAO:0100001",
            object_id="SECONDARY456",
            mapping_cardinality=MappingCardinality.ONE_TO_ONE,
        )
        result = mapping.to_sssom_dict()
        assert result["mapping_cardinality"] == "ONE_TO_ONE"


class TestBaseMappingToLegacyDict:
    """Tests for BaseMapping.to_legacy_dict method."""

    def test_to_legacy_dict_basic(self) -> None:
        """Test converting mapping to legacy dictionary."""
        mapping = BaseMapping(
            subject_id="PRIMARY123",
            predicate_id="IAO:0100001",
            object_id="SECONDARY456",
        )
        result = mapping.to_legacy_dict()
        assert result["primaryID"] == "PRIMARY123"
        assert result["secondaryID"] == "SECONDARY456"
        assert result["primaryLabel"] == ""
        assert result["secondaryLabel"] == ""
        assert result["predicateID"] == "IAO:0100001"

    def test_to_legacy_dict_with_labels(self) -> None:
        """Test that labels are mapped correctly."""
        mapping = BaseMapping(
            subject_id="PRIMARY123",
            predicate_id="IAO:0100001",
            object_id="SECONDARY456",
            subject_label="Primary Label",
            object_label="Secondary Label",
        )
        result = mapping.to_legacy_dict()
        assert result["primaryLabel"] == "Primary Label"
        assert result["secondaryLabel"] == "Secondary Label"
        assert result["primarySymbol"] == "Primary Label"
        assert result["secondarySymbol"] == "Secondary Label"


# =============================================================================
# Tests for Datasource-specific Mapping Classes
# =============================================================================


class TestChEBIMapping:
    """Tests for ChEBIMapping class."""

    def test_class_constants(self) -> None:
        """Test ChEBI class constants."""
        assert ChEBIMapping.datasource_name == "ChEBI"
        assert ChEBIMapping.datasource_prefix == "CHEBI"
        assert "obolibrary" in ChEBIMapping.curie_base_url

    def test_curie_method(self) -> None:
        """Test that curie adds CHEBI prefix."""
        result = ChEBIMapping.curie("12345")
        assert result == "CHEBI:12345"

    def test_curie_preserves_existing_prefix(self) -> None:
        """Test that existing prefix is preserved."""
        result = ChEBIMapping.curie("CHEBI:12345")
        assert result == "CHEBI:12345"


class TestHMDBMapping:
    """Tests for HMDBMapping class."""

    def test_class_constants(self) -> None:
        """Test HMDB class constants."""
        assert HMDBMapping.datasource_name == "HMDB"
        assert HMDBMapping.datasource_prefix == "HMDB"

    def test_curie_method(self) -> None:
        """Test that curie adds HMDB prefix."""
        result = HMDBMapping.curie("0000001")
        assert result == "HMDB:0000001"


class TestUniProtMapping:
    """Tests for UniProtMapping class."""

    def test_class_constants(self) -> None:
        """Test UniProt class constants."""
        assert UniProtMapping.datasource_name == "UniProt"
        assert UniProtMapping.datasource_prefix == "UniProt"


class TestHGNCMapping:
    """Tests for HGNCMapping class."""

    def test_class_constants(self) -> None:
        """Test HGNC class constants."""
        assert HGNCMapping.datasource_name == "HGNC"
        assert HGNCMapping.datasource_prefix == "HGNC"

    def test_curie_method(self) -> None:
        """Test that curie adds HGNC prefix."""
        result = HGNCMapping.curie("1234")
        assert result == "HGNC:1234"


class TestNCBIGeneMapping:
    """Tests for NCBIGeneMapping class."""

    def test_class_constants(self) -> None:
        """Test NCBI class constants."""
        assert NCBIGeneMapping.datasource_name == "NCBI"
        assert NCBIGeneMapping.datasource_prefix == "NCBIGene"


class TestIdMappingAndSymbolMapping:
    """Tests for IdMapping and SymbolMapping base classes."""

    def test_id_mapping_is_base_mapping(self) -> None:
        """Test IdMapping inherits from BaseMapping."""
        assert issubclass(IdMapping, BaseMapping)

    def test_symbol_mapping_is_base_mapping(self) -> None:
        """Test SymbolMapping inherits from BaseMapping."""
        assert issubclass(SymbolMapping, BaseMapping)

    def test_chebi_is_id_mapping(self) -> None:
        """Test ChEBIMapping is IdMapping."""
        assert issubclass(ChEBIMapping, IdMapping)

    def test_hgnc_is_symbol_mapping(self) -> None:
        """Test HGNCMapping is SymbolMapping."""
        assert issubclass(HGNCMapping, SymbolMapping)


# =============================================================================
# Tests for MappingSet
# =============================================================================


class TestMappingSet:
    """Tests for MappingSet container class."""

    def test_create_empty_mapping_set(self) -> None:
        """Test creating an empty mapping set."""
        ms = MappingSet(datasource_name="Test")
        assert len(ms) == 0
        assert ms.datasource_name == "Test"

    def test_create_mapping_set_with_mappings(self) -> None:
        """Test creating a mapping set with initial mappings."""
        mapping = ChEBIMapping(
            subject_id="12345",
            predicate_id="IAO:0100001",
            object_id="67890",
        )
        ms = MappingSet(datasource_name="ChEBI", mappings=[mapping])
        assert len(ms) == 1

    def test_add_mapping(self) -> None:
        """Test adding a single mapping."""
        ms = MappingSet(datasource_name="Test")
        mapping = BaseMapping(
            subject_id="A",
            predicate_id="P",
            object_id="B",
        )
        ms.add_mapping(mapping)
        assert len(ms) == 1

    def test_add_mappings(self) -> None:
        """Test adding multiple mappings."""
        ms = MappingSet(datasource_name="Test")
        mappings = [
            BaseMapping(subject_id="A", predicate_id="P", object_id="B"),
            BaseMapping(subject_id="C", predicate_id="P", object_id="D"),
        ]
        ms.add_mappings(mappings)
        assert len(ms) == 2

    def test_iter_mappings(self) -> None:
        """Test iterating over mappings."""
        ms = MappingSet(datasource_name="Test")
        mappings = [
            BaseMapping(subject_id="A", predicate_id="P", object_id="B"),
            BaseMapping(subject_id="C", predicate_id="P", object_id="D"),
        ]
        ms.add_mappings(mappings)

        result = list(ms.iter_mappings())
        assert len(result) == 2
        assert result[0].subject_id == "A"
        assert result[1].subject_id == "C"

    def test_mapping_set_metadata(self) -> None:
        """Test mapping set metadata fields."""
        ms = MappingSet(
            datasource_name="ChEBI",
            version="2024-01",
            mapping_set_id="https://example.com/mappings",
            mapping_set_description="Test mappings",
            mapping_date="2024-01-15",
            comment="Generated for testing",
        )
        assert ms.version == "2024-01"
        assert ms.mapping_set_id == "https://example.com/mappings"
        assert ms.mapping_set_description == "Test mappings"
        assert ms.mapping_date == "2024-01-15"
        assert ms.comment == "Generated for testing"

    def test_default_license(self) -> None:
        """Test default license URL."""
        ms = MappingSet(datasource_name="Test")
        assert "creativecommons.org" in ms.license_url

    def test_curie_map(self) -> None:
        """Test CURIE map handling."""
        ms = MappingSet(
            datasource_name="Test",
            curie_map={"TEST": "http://example.com/"},
        )
        assert ms.curie_map == {"TEST": "http://example.com/"}


# =============================================================================
# Tests for compute_cardinality
# =============================================================================


class TestComputeCardinality:
    """Tests for compute_cardinality function."""

    def test_one_to_one_mapping(self) -> None:
        """Test 1:1 mapping detection."""
        mappings = [("A", "1")]
        result = compute_cardinality(mappings)
        assert result[("A", "1")] == MappingCardinality.ONE_TO_ONE

    def test_many_to_one_mapping(self) -> None:
        """Test n:1 mapping detection (merge)."""
        # Multiple secondary IDs map to one primary
        mappings = [
            ("A", "1"),
            ("A", "2"),
        ]
        result = compute_cardinality(mappings)
        assert result[("A", "1")] == MappingCardinality.MANY_TO_ONE
        assert result[("A", "2")] == MappingCardinality.MANY_TO_ONE

    def test_one_to_many_mapping(self) -> None:
        """Test 1:n mapping detection (split)."""
        # One secondary ID maps to multiple primaries
        mappings = [
            ("A", "1"),
            ("B", "1"),
        ]
        result = compute_cardinality(mappings)
        assert result[("A", "1")] == MappingCardinality.ONE_TO_MANY
        assert result[("B", "1")] == MappingCardinality.ONE_TO_MANY

    def test_many_to_many_mapping(self) -> None:
        """Test n:n mapping detection."""
        # Complex relationship
        mappings = [
            ("A", "1"),
            ("A", "2"),
            ("B", "1"),
        ]
        result = compute_cardinality(mappings)
        # A appears twice, 1 appears twice -> complex
        assert result[("A", "1")] == MappingCardinality.MANY_TO_MANY
        assert result[("B", "1")] == MappingCardinality.ONE_TO_MANY
        assert result[("A", "2")] == MappingCardinality.MANY_TO_ONE

    def test_withdrawn_entry(self) -> None:
        """Test withdrawn entry detection."""
        mappings = [(WITHDRAWN_ENTRY, "OLD123")]
        result = compute_cardinality(mappings)
        assert result[(WITHDRAWN_ENTRY, "OLD123")] == MappingCardinality.ONE_TO_ZERO

    def test_empty_mappings(self) -> None:
        """Test empty mapping list."""
        result = compute_cardinality([])
        assert result == {}

    def test_none_object_id_skipped(self) -> None:
        """Test that None object_id is skipped."""
        mappings = [("A", None)]
        result = compute_cardinality(mappings)  # type: ignore[arg-type]
        assert len(result) == 0

    def test_empty_subject_id_is_withdrawn(self) -> None:
        """Test that empty subject_id is treated as withdrawn."""
        mappings = [("", "OLD123")]
        result = compute_cardinality(mappings)
        assert result[("", "OLD123")] == MappingCardinality.ONE_TO_ZERO

    def test_dash_subject_id_is_withdrawn(self) -> None:
        """Test that '-' subject_id is treated as withdrawn."""
        mappings = [("-", "OLD123")]  # ignore
        result = compute_cardinality(mappings)
        assert result[("-", "OLD123")] == MappingCardinality.ONE_TO_ZERO


# =============================================================================
# Tests for edge cases
# =============================================================================


class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    def test_mapping_with_special_characters(self) -> None:
        """Test mapping with special characters in IDs."""
        mapping = BaseMapping(
            subject_id="ID:with/special_chars",
            predicate_id="P",
            object_id="Another:ID",
        )
        assert mapping.subject_id == "ID:with/special_chars"

    def test_mapping_with_unicode(self) -> None:
        """Test mapping with unicode in labels."""
        mapping = BaseMapping(
            subject_id="A",
            predicate_id="P",
            object_id="B",
            subject_label="α-glucose",
            object_label="Ω-3 fatty acid",
        )
        assert mapping.subject_label == "α-glucose"
        assert mapping.object_label == "Ω-3 fatty acid"

    def test_large_mapping_set(self) -> None:
        """Test MappingSet with many mappings."""
        ms = MappingSet(datasource_name="Test")
        for i in range(1000):
            mapping = BaseMapping(
                subject_id=f"PRIMARY{i}",
                predicate_id="P",
                object_id=f"SECONDARY{i}",
            )
            ms.add_mapping(mapping)
        assert len(ms) == 1000
