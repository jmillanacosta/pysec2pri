"""Tests for pysec2pri.parsers module."""

from pathlib import Path

import pytest
from sssom_schema import Mapping

from pysec2pri.parsers import (
    ChEBIParser,
    HGNCParser,
    HMDBParser,
    NCBIParser,
    UniProtParser,
)
from pysec2pri.parsers.base import Sec2PriMappingSet

TEST_DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def chebi_sdf_path() -> Path:
    """Test data directory."""
    return TEST_DATA_DIR / "fake_chebi.sdf"


@pytest.fixture
def hmdb_xml_path() -> Path:
    """Test data directory."""
    return TEST_DATA_DIR / "fake_hmdb.xml"


@pytest.fixture
def hgnc_withdrawn_path() -> Path:
    """Test data directory."""
    return TEST_DATA_DIR / "fake_hgnc_withdrawn.tsv"


@pytest.fixture
def hgnc_complete_path() -> Path:
    """Test data directory."""
    return TEST_DATA_DIR / "fake_hgnc_complete.tsv"


@pytest.fixture
def ncbi_history_path() -> Path:
    """Test data directory."""
    return TEST_DATA_DIR / "fake_gene_history.tsv"


@pytest.fixture
def ncbi_info_path() -> Path:
    """Test data directory."""
    return TEST_DATA_DIR / "fake_gene_info.tsv"


@pytest.fixture
def uniprot_sec_ac_path() -> Path:
    """Test data directory."""
    return TEST_DATA_DIR / "fake_sec_ac.txt"


class TestChEBIParser:
    """Tests for ChEBI parser."""

    def test_parse(self, chebi_sdf_path: Path) -> None:
        """Test ChEBI parsing output."""
        parser = ChEBIParser(show_progress=False)
        result = parser.parse(chebi_sdf_path)
        assert isinstance(result, Sec2PriMappingSet)
        assert len(result.mappings) > 0

    def test_parse_synonyms(self, chebi_sdf_path: Path) -> None:
        """Test parsing output for ChEBI synonyms."""
        parser = ChEBIParser(show_progress=False)
        result = parser.parse_synonyms(chebi_sdf_path)
        assert isinstance(result, Sec2PriMappingSet)


class TestHMDBParser:
    """Tests for HMDB parser."""

    def test_parse(self, hmdb_xml_path: Path) -> None:
        """Test HMDB parsing output."""
        parser = HMDBParser(show_progress=False)
        result = parser.parse(hmdb_xml_path)
        assert isinstance(result, Sec2PriMappingSet)


class TestHGNCParser:
    """Tests for HGNC parser."""

    def test_parse(self, hgnc_withdrawn_path: Path) -> None:
        """Test parsing output for HGNC IDs."""
        parser = HGNCParser(show_progress=False)
        result = parser.parse(hgnc_withdrawn_path)
        assert isinstance(result, Sec2PriMappingSet)
        assert len(result.mappings) > 0

    def test_parse_symbols(self, hgnc_complete_path: Path) -> None:
        """Test parsing output for HGNC symbols."""
        parser = HGNCParser(show_progress=False)
        result = parser.parse_symbols(hgnc_complete_path)
        assert isinstance(result, Sec2PriMappingSet)


class TestNCBIParser:
    """Tests for NCBI parser."""

    def test_parse(self, ncbi_history_path: Path) -> None:
        """Test parsing output for NCBI IDs."""
        parser = NCBIParser(show_progress=False)
        result = parser.parse(ncbi_history_path, tax_id="9606")
        assert isinstance(result, Sec2PriMappingSet)
        assert len(result.mappings) > 0

    def test_parse_symbols(self, ncbi_info_path: Path) -> None:
        """Test parsing output for NCBI symbols."""
        parser = NCBIParser(show_progress=False)
        result = parser.parse_symbols(ncbi_info_path, tax_id="9606")
        assert isinstance(result, Sec2PriMappingSet)


class TestUniProtParser:
    """Tests for UniProt parser."""

    def test_parse(self, uniprot_sec_ac_path: Path) -> None:
        """Test parsing output for UniProt."""
        parser = UniProtParser(show_progress=False)
        result = parser.parse(uniprot_sec_ac_path)
        assert isinstance(result, Sec2PriMappingSet)
        assert len(result.mappings) >= 4


class TestParserIntegration:
    """Integration tests for all parsers."""

    def test_all_mappings_are_mapping_instances(
        self,
        chebi_sdf_path: Path,
        hmdb_xml_path: Path,
        hgnc_withdrawn_path: Path,
        ncbi_history_path: Path,
        uniprot_sec_ac_path: Path,
    ) -> None:
        """Test that all mappings are mapping instances."""
        results = [
            ChEBIParser(show_progress=False).parse(chebi_sdf_path),
            HMDBParser(show_progress=False).parse(hmdb_xml_path),
            HGNCParser(show_progress=False).parse(hgnc_withdrawn_path),
            NCBIParser(show_progress=False).parse(ncbi_history_path, tax_id="9606"),
            UniProtParser(show_progress=False).parse(uniprot_sec_ac_path),
        ]

        for result in results:
            assert isinstance(result, Sec2PriMappingSet)
            for m in result.mappings:
                assert isinstance(m, Mapping)
                assert m.subject_id is not None
                assert m.predicate_id is not None
