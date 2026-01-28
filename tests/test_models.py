"""Tests for pysec2pri data models."""

from __future__ import annotations

import pytest

from pysec2pri.models import (
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


class TestMappingCardinality:
    """Tests for MappingCardinality enum."""

    def test_cardinality_values(self) -> None:
        """Test cardinality enum values."""
        assert MappingCardinality.ONE_TO_ONE.value == "1:1"
        assert MappingCardinality.MANY_TO_ONE.value == "n:1"
        assert MappingCardinality.ONE_TO_MANY.value == "1:n"
        assert MappingCardinality.MANY_TO_MANY.value == "n:n"
        assert MappingCardinality.ONE_TO_ZERO.value == "1:0"


class TestPredicateID:
    """Tests for PredicateID enum."""

    def test_predicate_values(self) -> None:
        """Test predicate enum values."""
        assert PredicateID.TERM_REPLACED_BY.value == "IAO:0100001"
        assert PredicateID.CONSIDER.value == "oboInOwl:consider"


class TestGetPredicateForCardinality:
    """Tests for get_predicate_for_cardinality function."""

    def test_one_to_one_returns_replaced_by(self) -> None:
        """Test 1:1 returns TERM_REPLACED_BY."""
        assert (
            get_predicate_for_cardinality(MappingCardinality.ONE_TO_ONE)
            == PredicateID.TERM_REPLACED_BY
        )

    def test_many_to_one_returns_replaced_by(self) -> None:
        """Test n:1 returns TERM_REPLACED_BY."""
        assert (
            get_predicate_for_cardinality(MappingCardinality.MANY_TO_ONE)
            == PredicateID.TERM_REPLACED_BY
        )

    def test_one_to_many_returns_consider(self) -> None:
        """Test 1:n returns CONSIDER."""
        assert (
            get_predicate_for_cardinality(MappingCardinality.ONE_TO_MANY)
            == PredicateID.CONSIDER
        )

    def test_many_to_many_returns_consider(self) -> None:
        """Test n:n returns CONSIDER."""
        assert (
            get_predicate_for_cardinality(MappingCardinality.MANY_TO_MANY)
            == PredicateID.CONSIDER
        )

    def test_one_to_zero_returns_consider(self) -> None:
        """Test 1:0 returns CONSIDER."""
        assert (
            get_predicate_for_cardinality(MappingCardinality.ONE_TO_ZERO)
            == PredicateID.CONSIDER
        )

    def test_none_returns_none(self) -> None:
        """Test None returns None."""
        assert get_predicate_for_cardinality(None) is None


class TestGetCommentForCardinality:
    """Tests for get_comment_for_cardinality function."""

    def test_one_to_one_comment(self) -> None:
        """Test 1:1 comment."""
        comment = get_comment_for_cardinality(MappingCardinality.ONE_TO_ONE)
        assert "replaced" in comment.lower()

    def test_one_to_zero_comment(self) -> None:
        """Test 1:0 comment."""
        comment = get_comment_for_cardinality(MappingCardinality.ONE_TO_ZERO)
        assert "withdrawn" in comment.lower() or "deprecated" in comment.lower()

    def test_none_returns_empty_string(self) -> None:
        """Test None returns empty string."""
        assert get_comment_for_cardinality(None) == ""


class TestIdMapping:
    """Tests for IdMapping class."""

    def test_create_basic_mapping(self) -> None:
        """Test creating a basic IdMapping."""
        mapping = ChEBIMapping(
            primary_id="12345",
            secondary_id="67890",
            mapping_cardinality=MappingCardinality.ONE_TO_ONE,
        )
        assert mapping.primary_id == "12345"
        assert mapping.secondary_id == "67890"
        assert mapping.mapping_cardinality == MappingCardinality.ONE_TO_ONE

    def test_computed_predicate_id(self) -> None:
        """Test computed predicate_id property."""
        mapping = ChEBIMapping(
            primary_id="12345",
            secondary_id="67890",
            mapping_cardinality=MappingCardinality.ONE_TO_ONE,
        )
        assert mapping.predicate_id == PredicateID.TERM_REPLACED_BY

    def test_is_withdrawn_property(self) -> None:
        """Test is_withdrawn computed property."""
        withdrawn = ChEBIMapping(
            primary_id="12345",
            mapping_cardinality=MappingCardinality.ONE_TO_ZERO,
        )
        active = ChEBIMapping(
            primary_id="12345",
            secondary_id="67890",
            mapping_cardinality=MappingCardinality.ONE_TO_ONE,
        )
        assert withdrawn.is_withdrawn is True
        assert active.is_withdrawn is False

    def test_curie_method(self) -> None:
        """Test CURIE formatting."""
        mapping = ChEBIMapping(primary_id="12345")
        curie = mapping.curie("12345")
        assert curie == "CHEBI:12345"

    def test_curie_already_prefixed(self) -> None:
        """Test CURIE when already prefixed."""
        mapping = ChEBIMapping(primary_id="12345")
        curie = mapping.curie("CHEBI:12345")
        assert curie == "CHEBI:12345"

    def test_to_sssom_dict(self) -> None:
        """Test SSSOM dict conversion."""
        mapping = ChEBIMapping(
            primary_id="12345",
            secondary_id="67890",
            mapping_cardinality=MappingCardinality.ONE_TO_ONE,
            primary_label="Test Compound",
        )
        sssom = mapping.to_sssom_dict()
        assert sssom["subject_id"] == "CHEBI:12345"
        assert sssom["object_id"] == "CHEBI:67890"
        assert sssom["predicate_id"] == "IAO:0100001"
        assert sssom["subject_label"] == "Test Compound"


class TestSymbolMapping:
    """Tests for SymbolMapping class."""

    def test_create_symbol_mapping(self) -> None:
        """Test creating a SymbolMapping with symbols."""
        mapping = HGNCMapping(
            primary_id="HGNC:1234",
            secondary_id="HGNC:5678",
            primary_symbol="TP53",
            secondary_symbol="P53",
            mapping_cardinality=MappingCardinality.ONE_TO_ONE,
        )
        assert mapping.primary_symbol == "TP53"
        assert mapping.secondary_symbol == "P53"

    def test_to_sssom_dict_with_symbols(self) -> None:
        """Test SSSOM dict includes symbols as labels."""
        mapping = HGNCMapping(
            primary_id="HGNC:1234",
            secondary_id="HGNC:5678",
            primary_symbol="TP53",
            secondary_symbol="P53",
            mapping_cardinality=MappingCardinality.ONE_TO_ONE,
        )
        sssom = mapping.to_sssom_dict()
        assert sssom["subject_label"] == "TP53"
        assert sssom["object_label"] == "P53"


class TestDatasourceSpecificMappings:
    """Tests for datasource-specific mapping classes."""

    def test_chebi_mapping_metadata(self) -> None:
        """Test ChEBIMapping class metadata."""
        assert ChEBIMapping.datasource_name == "ChEBI"
        assert ChEBIMapping.datasource_prefix == "CHEBI"
        assert "obolibrary" in ChEBIMapping.curie_base_url

    def test_hmdb_mapping_metadata(self) -> None:
        """Test HMDBMapping class metadata."""
        assert HMDBMapping.datasource_name == "HMDB"
        assert HMDBMapping.datasource_prefix == "HMDB"

    def test_uniprot_mapping_metadata(self) -> None:
        """Test UniProtMapping class metadata."""
        assert UniProtMapping.datasource_name == "UniProt"
        assert UniProtMapping.datasource_prefix == "UniProt"

    def test_hgnc_mapping_metadata(self) -> None:
        """Test HGNCMapping class metadata."""
        assert HGNCMapping.datasource_name == "HGNC"
        assert HGNCMapping.datasource_prefix == "HGNC"

    def test_ncbi_mapping_metadata(self) -> None:
        """Test NCBIGeneMapping class metadata."""
        assert NCBIGeneMapping.datasource_name == "NCBI"
        assert NCBIGeneMapping.datasource_prefix == "NCBIGene"


class TestMappingSet:
    """Tests for MappingSet container class."""

    def test_create_empty_mapping_set(self) -> None:
        """Test creating an empty MappingSet."""
        ms = MappingSet(datasource_name="ChEBI")
        assert len(ms) == 0
        assert ms.datasource_name == "ChEBI"

    def test_add_mapping(self) -> None:
        """Test adding a single mapping."""
        ms = MappingSet(datasource_name="ChEBI")
        mapping = ChEBIMapping(primary_id="12345")
        ms.add_mapping(mapping)
        assert len(ms) == 1

    def test_add_mappings(self) -> None:
        """Test adding multiple mappings."""
        ms = MappingSet(datasource_name="ChEBI")
        mappings = [
            ChEBIMapping(primary_id="12345"),
            ChEBIMapping(primary_id="67890"),
        ]
        ms.add_mappings(mappings)
        assert len(ms) == 2

    def test_iter_mappings(self) -> None:
        """Test iterating over mappings."""
        ms = MappingSet(datasource_name="ChEBI")
        mappings = [
            ChEBIMapping(primary_id="12345"),
            ChEBIMapping(primary_id="67890"),
        ]
        ms.add_mappings(mappings)

        count = 0
        for mapping in ms.iter_mappings():
            count += 1
            assert isinstance(mapping, ChEBIMapping)
        assert count == 2

    def test_mapping_set_with_metadata(self) -> None:
        """Test MappingSet with full metadata."""
        ms = MappingSet(
            datasource_name="ChEBI",
            version="v220",
            mapping_set_id="https://example.com/chebi_sec2pri",
            mapping_set_description="ChEBI secondary to primary mappings",
            mapping_date="2024-01-01",
        )
        assert ms.version == "v220"
        assert ms.mapping_set_id == "https://example.com/chebi_sec2pri"
        assert ms.mapping_date == "2024-01-01"


class TestComputeCardinality:
    """Tests for compute_cardinality function."""

    def test_one_to_one(self) -> None:
        """Test 1:1 cardinality detection."""
        mappings = [
            ChEBIMapping(primary_id="100", secondary_id="1"),
            ChEBIMapping(primary_id="200", secondary_id="2"),
        ]
        cardinality = compute_cardinality("1", "100", mappings)
        assert cardinality == MappingCardinality.ONE_TO_ONE

    def test_many_to_one(self) -> None:
        """Test n:1 cardinality detection (multiple secondaries to one primary)."""
        mappings = [
            ChEBIMapping(primary_id="100", secondary_id="1"),
            ChEBIMapping(primary_id="100", secondary_id="2"),
        ]
        cardinality = compute_cardinality("1", "100", mappings)
        assert cardinality == MappingCardinality.MANY_TO_ONE

    def test_one_to_many(self) -> None:
        """Test 1:n cardinality detection (one secondary to multiple primaries)."""
        mappings = [
            ChEBIMapping(primary_id="100", secondary_id="1"),
            ChEBIMapping(primary_id="200", secondary_id="1"),
        ]
        cardinality = compute_cardinality("1", "100", mappings)
        assert cardinality == MappingCardinality.ONE_TO_MANY

    def test_one_to_zero(self) -> None:
        """Test 1:0 cardinality detection (withdrawn)."""
        mappings = [
            ChEBIMapping(primary_id="100", secondary_id=None),
        ]
        cardinality = compute_cardinality("100", "100", mappings)
        assert cardinality == MappingCardinality.ONE_TO_ZERO
