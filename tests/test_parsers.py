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
    return TEST_DATA_DIR / "mock_chebi.sdf"


@pytest.fixture
def chebi_names_tsv_path() -> Path:
    """Fake ChEBI names.tsv test data."""
    return TEST_DATA_DIR / "mock_chebi_names.tsv"


@pytest.fixture
def chebi_compounds_tsv_path() -> Path:
    """Fake ChEBI compounds.tsv test data."""
    return TEST_DATA_DIR / "mock_chebi_compounds.tsv"


@pytest.fixture
def hmdb_xml_path() -> Path:
    """Test data directory."""
    return TEST_DATA_DIR / "mock_hmdb.xml"


@pytest.fixture
def hmdb_proteins_xml_path() -> Path:
    """Fake HMDB proteins XML test data."""
    return TEST_DATA_DIR / "mock_hmdb_proteins.xml"


@pytest.fixture
def hgnc_withdrawn_path() -> Path:
    """Test data directory."""
    return TEST_DATA_DIR / "mock_hgnc_withdrawn.tsv"


@pytest.fixture
def hgnc_complete_path() -> Path:
    """Test data directory."""
    return TEST_DATA_DIR / "mock_hgnc_complete.tsv"


@pytest.fixture
def ncbi_history_path() -> Path:
    """Test data directory."""
    return TEST_DATA_DIR / "mock_gene_history.tsv"


@pytest.fixture
def ncbi_info_path() -> Path:
    """Test data directory."""
    return TEST_DATA_DIR / "mock_gene_info.tsv"


@pytest.fixture
def uniprot_sec_ac_path() -> Path:
    """Test data directory."""
    return TEST_DATA_DIR / "mock_sec_ac.txt"


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

    def test_parse_synonyms_direction(self, chebi_sdf_path: Path) -> None:
        """Synonym mappings must be sec:pri: synonym is subject, primary name is object."""
        parser = ChEBIParser(show_progress=False)
        result = parser.parse_synonyms(chebi_sdf_path)
        mappings = result.mappings or []
        assert len(mappings) > 0
        subject_labels = {m.subject_label for m in mappings}
        object_labels = {m.object_label for m in mappings}
        # mock_chebi.sdf: CHEBI:10001 has primary name "alpha-glucose"
        # and synonyms "glucose", "D-glucose": synonyms must be subjects
        assert "glucose" in subject_labels
        assert "D-glucose" in subject_labels
        assert "alpha-glucose" in object_labels
        # the primary/canonical name must never appear as a subject_label
        assert "alpha-glucose" not in subject_labels

    def test_parse_synonyms_tsv_direction(
        self,
        chebi_names_tsv_path: Path,
        chebi_compounds_tsv_path: Path,
    ) -> None:
        """TSV path: primary name comes from compounds.tsv, NOT from a names.tsv heuristic."""
        parser = ChEBIParser(show_progress=False, subset="3star")
        result = parser.parse_synonyms(
            names_path=chebi_names_tsv_path,
            compounds_path=chebi_compounds_tsv_path,
        )
        mappings = result.mappings or []
        assert len(mappings) > 0
        subject_labels = {m.subject_label for m in mappings}
        object_labels = {m.object_label for m in mappings}
        # "alpha-glucose" is the compounds.tsv canonical name for CHEBI:10001
        # : must be object_label only; its names.tsv entry is excluded as self-mapping
        assert "alpha-glucose" in object_labels
        assert "alpha-glucose" not in subject_labels
        # All other names.tsv entries are synonyms : subject_label
        assert "glucose" in subject_labels
        assert "D-glucose" in subject_labels
        # Same logic for CHEBI:10002 and CHEBI:10003
        assert "beta-glucose" in object_labels
        assert "B-glucose" in subject_labels
        assert "water" in object_labels
        assert "H2O" in subject_labels
        assert "dihydrogen oxide" in subject_labels


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

    def test_parse_populates_primary_ids_including_no_secondary(self, hmdb_xml_path: Path) -> None:
        """_primary_ids includes ALL metabolites, even those with no secondaries."""
        parser = HMDBParser(show_progress=False)
        result = parser.parse(hmdb_xml_path)
        pri_ids = result.to_pri_ids()
        # fake file has 3 metabolites; HMDB0000003 has no secondaries but must appear
        assert "HMDB:HMDB0000001" in pri_ids
        assert "HMDB:HMDB0000002" in pri_ids
        assert "HMDB:HMDB0000003" in pri_ids


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

    def test_parse_with_gene_info_populates_primary_ids(
        self, ncbi_history_path: Path, ncbi_info_path: Path
    ) -> None:
        """Passing gene_info_path should populate _primary_ids with all current IDs."""
        parser = NCBIParser(show_progress=False)
        result = parser.parse(ncbi_history_path, tax_id="9606", gene_info_path=ncbi_info_path)
        assert isinstance(result, Sec2PriMappingSet)
        pri_ids = result.to_pri_ids()
        # _primary_ids comes from gene_info, so it should include IDs even for
        # genes that have no secondary in gene_history
        assert len(pri_ids) > 0
        assert all(id_.startswith("NCBIGene:") for id_ in pri_ids)

        """Test parsing output for NCBI symbols."""
        parser = NCBIParser(show_progress=False)
        result = parser.parse_symbols(ncbi_info_path, tax_id="9606")
        assert isinstance(result, Sec2PriMappingSet)

    def test_parse_symbols_direction(self, ncbi_info_path: Path) -> None:
        """Symbol mappings must be sec:pri: synonym is subject, primary symbol is object."""
        parser = NCBIParser(show_progress=False)
        result = parser.parse_symbols(ncbi_info_path, tax_id="9606")
        mappings = result.mappings or []
        assert len(mappings) > 0
        subject_labels = {m.subject_label for m in mappings}
        object_labels = {m.object_label for m in mappings}
        # mock_gene_info.tsv: GENE1 has synonyms ALT_SYM1 and ALT_SYM2
        # synonyms must be subjects; current Symbol must be the object
        assert "ALT_SYM1" in subject_labels
        assert "ALT_SYM2" in subject_labels
        assert "GENE1" in object_labels
        # the current primary symbol must never appear as a subject_label
        assert "GENE1" not in subject_labels


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
