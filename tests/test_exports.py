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


@pytest.fixture
def prev_label_mapping_set() -> LabelMappingSet:
    """LabelMappingSet with IAO:0100001 (previous-symbol / deprecation) rows."""
    mappings = [
        Mapping(
            subject_id="HGNC:1",
            subject_label="PSCP",
            object_id="HGNC:1",
            object_label="BRCA1",
            predicate_id="IAO:0100001",
            mapping_justification="semapv:BackgroundKnowledgeBasedMatching",
        ),
        Mapping(
            subject_id="HGNC:2",
            subject_label="tumor_p53",
            object_id="HGNC:2",
            object_label="TP53",
            predicate_id="IAO:0100001",
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
            # header now uses primary_id and secondary_id
            assert "primary_id" in lines[0]
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
            # header now uses name/synonym with primary_id
            assert "name" in lines[0]


class TestWriteSymbol2Prev:
    """Tests for write_symbol2prev function."""

    def test_write_symbol2prev(self, prev_label_mapping_set: LabelMappingSet) -> None:
        """Test writing symbol2prev mappings to TSV file.

        Uses a mapping set with IAO:0100001 rows only; write_symbol2prev
        must skip hasExactSynonym rows and emit only deprecation rows.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.tsv"
            result = write_symbol2prev(prev_label_mapping_set, output_path)

            assert result.exists()
            content = result.read_text()
            lines = content.strip().split("\n")

            assert len(lines) == 3
            # header includes mapping_cardinality and uses primary_symbol
            assert "mapping_cardinality" in lines[0]
            assert "primary_symbol" in lines[0]

    def test_write_symbol2prev_skips_synonyms(self, label_mapping_set: LabelMappingSet) -> None:
        """write_symbol2prev emits only the header when all rows are hasExactSynonym."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.tsv"
            write_symbol2prev(label_mapping_set, output_path)
            lines = output_path.read_text().strip().split("\n")
            # Only the header row — no synonym rows should appear
            assert len(lines) == 1
            assert "primary_symbol" in lines[0]
