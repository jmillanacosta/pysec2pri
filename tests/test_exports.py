"""Tests for pysec2pri.exports module."""

import tempfile
from pathlib import Path

import pytest
from sssom_schema import Mapping

from pysec2pri.exports import (
    write_name2synonym,
    write_sec2pri,
    write_symbol2prev,
)
from pysec2pri.parsers.base import IdMappingSet, LabelMappingSet


@pytest.fixture
def id_mapping_set() -> IdMappingSet:
    """Create a simple IdMappingSet for testing."""
    mappings = [
        Mapping(
            subject_id="CHEBI:10001",
            object_id="CHEBI:99901",
            predicate_id="IAO:0100001",
            mapping_justification="semapv:BackgroundKnowledgeBasedMatching",
        ),
        Mapping(
            subject_id="CHEBI:10002",
            object_id="CHEBI:99902",
            predicate_id="IAO:0100001",
            mapping_justification="semapv:BackgroundKnowledgeBasedMatching",
        ),
    ]
    ms = IdMappingSet(
        mapping_set_id="https://example.com/test",
        license="https://creativecommons.org/publicdomain/zero/1.0/",
        mappings=mappings,
    )
    ms.compute_cardinalities()
    return ms


@pytest.fixture
def label_mapping_set() -> LabelMappingSet:
    """Create a simple LabelMappingSet for testing."""
    mappings = [
        Mapping(
            subject_id="HGNC:1",
            subject_label="BRCA1",
            object_label="BRCC1",
            predicate_id="oboInOwl:hasExactSynonym",
            mapping_justification="semapv:BackgroundKnowledgeBasedMatching",
        ),
        Mapping(
            subject_id="HGNC:2",
            subject_label="TP53",
            object_label="P53",
            predicate_id="oboInOwl:hasExactSynonym",
            mapping_justification="semapv:BackgroundKnowledgeBasedMatching",
        ),
    ]
    ms = LabelMappingSet(
        mapping_set_id="https://example.com/test",
        license="https://creativecommons.org/publicdomain/zero/1.0/",
        mappings=mappings,
    )
    ms.compute_cardinalities()
    return ms


class TestWriteSec2Pri:
    """Tests for write_sec2pri function."""

    def test_write_sec2pri(self, id_mapping_set: IdMappingSet) -> None:
        """Test writing sec2pri mappings to TSV file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.tsv"
            result = write_sec2pri(id_mapping_set, output_path)

            assert result.exists()
            content = result.read_text()
            lines = content.strip().split("\n")

            # Header + 2 data rows
            assert len(lines) == 3
            assert "subject_id" in lines[0]
            assert "CHEBI:10001" in lines[1]


class TestWriteName2Synonym:
    """Tests for write_name2synonym function."""

    def test_write_name2synonym(self, label_mapping_set: LabelMappingSet) -> None:
        """Test writing name2synonym mappings to TSV file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.tsv"
            result = write_name2synonym(label_mapping_set, output_path)

            assert result.exists()
            content = result.read_text()
            lines = content.strip().split("\n")

            # Header + 2 data rows
            assert len(lines) == 3
            assert "subject_label" in lines[0]


class TestWriteSymbol2Prev:
    """Tests for write_symbol2prev function."""

    def test_write_symbol2prev(self, label_mapping_set: LabelMappingSet) -> None:
        """Test writing symbol2prev mappings to TSV file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.tsv"
            result = write_symbol2prev(label_mapping_set, output_path)

            assert result.exists()
            content = result.read_text()
            lines = content.strip().split("\n")

            assert len(lines) == 3
            assert "mapping_cardinality" in lines[0]
