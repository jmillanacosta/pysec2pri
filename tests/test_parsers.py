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
def hmdb_proteins_xml_path() -> Path:
    """Fake HMDB proteins XML test data."""
    return TEST_DATA_DIR / "fake_hmdb_proteins.xml"


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

    def test_parse_returns_mapping_set(self, hmdb_xml_path: Path) -> None:
        """parse() returns a Sec2PriMappingSet."""
        parser = HMDBParser(show_progress=False)
        result = parser.parse(hmdb_xml_path)
        assert isinstance(result, Sec2PriMappingSet)

    def test_parse_produces_mappings(self, hmdb_xml_path: Path) -> None:
        """parse() extracts sec->pri mappings (fake file has 3)."""
        parser = HMDBParser(show_progress=False)
        result = parser.parse(hmdb_xml_path)
        assert len(result.mappings or []) == 3

    def test_parse_subject_object_ids(self, hmdb_xml_path: Path) -> None:
        """Secondary IDs are subjects; primary IDs are objects."""
        parser = HMDBParser(show_progress=False)
        result = parser.parse(hmdb_xml_path)
        mappings = result.mappings or []
        subjects = {m.subject_id for m in mappings}
        objects = {m.object_id for m in mappings}
        assert "HMDB:HMDB00001" in subjects
        assert "HMDB:HMDB0001001" in subjects
        assert "HMDB:HMDB0000001" in objects
        assert "HMDB:HMDB00002" in subjects
        assert "HMDB:HMDB0000002" in objects

    def test_parse_no_secondary_skipped(self, hmdb_xml_path: Path) -> None:
        """Records with empty secondary_accessions produce no mappings."""
        parser = HMDBParser(show_progress=False)
        result = parser.parse(hmdb_xml_path)
        subjects = {m.subject_id for m in (result.mappings or [])}
        # HMDB0000003 has no secondary accessions
        assert "HMDB:HMDB0000003" not in subjects

    def test_parse_proteins_returns_mapping_set(self, hmdb_proteins_xml_path: Path) -> None:
        """parse_proteins() returns a Sec2PriMappingSet."""
        parser = HMDBParser(show_progress=False)
        result = parser.parse_proteins(hmdb_proteins_xml_path)
        assert isinstance(result, Sec2PriMappingSet)

    def test_parse_proteins_produces_mappings(self, hmdb_proteins_xml_path: Path) -> None:
        """parse_proteins() extracts sec->pri mappings (fake has 3)."""
        parser = HMDBParser(show_progress=False)
        result = parser.parse_proteins(hmdb_proteins_xml_path)
        assert len(result.mappings or []) == 3

    def test_parse_proteins_bare_number_normalised(self, hmdb_proteins_xml_path: Path) -> None:
        """Bare numeric secondary accessions are normalised to HMDBP prefix."""
        parser = HMDBParser(show_progress=False)
        result = parser.parse_proteins(hmdb_proteins_xml_path)
        subjects = {m.subject_id for m in (result.mappings or [])}
        # "5229" should become "HMDB:HMDBP05229" (secondary to subject)
        assert "HMDB:HMDBP05229" in subjects
        # full HMDBP accession unchanged
        assert "HMDB:HMDBP05261" in subjects

    def test_parse_proteins_no_secondary_skipped(self, hmdb_proteins_xml_path: Path) -> None:
        """Protein records with empty secondary_accessions produce no mappings."""
        parser = HMDBParser(show_progress=False)
        result = parser.parse_proteins(hmdb_proteins_xml_path)
        subjects = {m.subject_id for m in (result.mappings or [])}
        # HMDBP00003 has no secondary accessions
        assert "HMDB:HMDBP00003" not in subjects


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

    def test_predicate_label_on_sec2pri_mappings(
        self,
        chebi_sdf_path: Path,
        hmdb_xml_path: Path,
        hgnc_withdrawn_path: Path,
        ncbi_history_path: Path,
        uniprot_sec_ac_path: Path,
    ) -> None:
        """Sec2pri mappings use IAO:0100001 with label 'term replaced by'."""
        results = [
            ChEBIParser(show_progress=False).parse(chebi_sdf_path),
            HMDBParser(show_progress=False).parse(hmdb_xml_path),
            HGNCParser(show_progress=False).parse(hgnc_withdrawn_path),
            NCBIParser(show_progress=False).parse(ncbi_history_path, tax_id="9606"),
            UniProtParser(show_progress=False).parse(uniprot_sec_ac_path),
        ]
        for result in results:
            for m in result.mappings:
                if m.predicate_id == "IAO:0100001":
                    assert m.predicate_label == "term replaced by", (
                        f"predicate_label should be 'term replaced by' "
                        f"for IAO:0100001, got {m.predicate_label!r}"
                    )
