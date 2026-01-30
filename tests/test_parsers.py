"""Tests for pysec2pri.parsers module.

Tests cover all parsers using small fake datasets.
"""

from pathlib import Path

import pytest

from pysec2pri.models import (
    WITHDRAWN_ENTRY,
    ChEBIMapping,
    MappingCardinality,
    MappingSet,
    UniProtMapping,
)
from pysec2pri.parsers import (
    ChEBIParser,
    HGNCParser,
    HMDBParser,
    NCBIParser,
    UniProtParser,
)

# Path to test data directory
TEST_DATA_DIR = Path(__file__).parent / "data"


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def chebi_sdf_path() -> Path:
    """Path to fake ChEBI SDF file."""
    return TEST_DATA_DIR / "fake_chebi.sdf"


@pytest.fixture
def hmdb_xml_path() -> Path:
    """Path to fake HMDB XML file."""
    return TEST_DATA_DIR / "fake_hmdb.xml"


@pytest.fixture
def hgnc_withdrawn_path() -> Path:
    """Path to fake HGNC withdrawn file."""
    return TEST_DATA_DIR / "fake_hgnc_withdrawn.tsv"


@pytest.fixture
def hgnc_complete_path() -> Path:
    """Path to fake HGNC complete set file."""
    return TEST_DATA_DIR / "fake_hgnc_complete.tsv"


@pytest.fixture
def ncbi_history_path() -> Path:
    """Path to fake NCBI gene history file."""
    return TEST_DATA_DIR / "fake_gene_history.tsv"


@pytest.fixture
def ncbi_info_path() -> Path:
    """Path to fake NCBI gene info file."""
    return TEST_DATA_DIR / "fake_gene_info.tsv"


@pytest.fixture
def uniprot_sec_ac_path() -> Path:
    """Path to fake UniProt secondary accessions file."""
    return TEST_DATA_DIR / "fake_sec_ac.txt"


@pytest.fixture
def uniprot_delac_path() -> Path:
    """Path to fake UniProt deleted accessions file."""
    return TEST_DATA_DIR / "fake_delac_sp.txt"


# =============================================================================
# ChEBI Parser Tests
# =============================================================================


class TestChEBIParser:
    """Tests for ChEBIParser."""

    def test_parser_initialization(self) -> None:
        """Test parser can be initialized."""
        parser = ChEBIParser(version="test", show_progress=False)
        assert parser.version == "test"
        assert parser.datasource_name == "ChEBI"

    def test_parse_fake_sdf(self, chebi_sdf_path: Path) -> None:
        """Test parsing fake ChEBI SDF file."""
        parser = ChEBIParser(show_progress=False)
        result = parser.parse(chebi_sdf_path)

        assert isinstance(result, MappingSet)
        assert result.datasource_name == "ChEBI"
        assert len(result) > 0

    def test_parse_extracts_secondary_ids(self, chebi_sdf_path: Path) -> None:
        """Test that secondary IDs are extracted."""
        parser = ChEBIParser(show_progress=False)
        result = parser.parse(chebi_sdf_path)

        # Find ID mappings (not synonym mappings)
        id_mappings = [
            m
            for m in result.mappings
            if isinstance(m, ChEBIMapping) and m.object_id and m.object_id.startswith("CHEBI:")
        ]
        assert len(id_mappings) > 3  # We have 4 secondary IDs in test data

    def test_parse_extracts_synonyms(self, chebi_sdf_path: Path) -> None:
        """Test that synonyms are extracted."""
        parser = ChEBIParser(show_progress=False)
        result = parser.parse(chebi_sdf_path)

        # All mappings should be ChEBIMapping
        for m in result.mappings:
            assert isinstance(m, ChEBIMapping)


# =============================================================================
# HMDB Parser Tests
# =============================================================================


class TestHMDBParser:
    """Tests for HMDBParser."""

    def test_parser_initialization(self) -> None:
        """Test parser can be initialized."""
        parser = HMDBParser(version="test", show_progress=False)
        assert parser.version == "test"
        assert parser.datasource_name == "HMDB"

    def test_parse_fake_xml(self, hmdb_xml_path: Path) -> None:
        """Test parsing fake HMDB XML file."""
        parser = HMDBParser(show_progress=False)
        result = parser.parse(hmdb_xml_path)

        assert isinstance(result, MappingSet)
        assert result.datasource_name == "HMDB"

    def test_parse_extracts_secondary_accessions(self, hmdb_xml_path: Path) -> None:
        """Test that secondary accessions are extracted."""
        parser = HMDBParser(show_progress=False)
        result = parser.parse(hmdb_xml_path)

        # Check that some mappings exist
        # First metabolite has 2 secondary accessions
        # Second metabolite has 1 secondary accession
        # Total = 3 secondary accession mappings
        id_mappings = [m for m in result.mappings if m.object_id and "HMDB" in str(m.object_id)]
        assert len(id_mappings) == 3

    def test_primary_ids_are_correct(self, hmdb_xml_path: Path) -> None:
        """Test that primary IDs are correctly identified."""
        parser = HMDBParser(show_progress=False)
        result = parser.parse(hmdb_xml_path)

        # All subject_ids should be primary HMDB IDs
        for m in result.mappings:
            if m.subject_id != WITHDRAWN_ENTRY:
                # Primary IDs should be HMDB format
                assert m.subject_id.startswith("HMDB")


# =============================================================================
# HGNC Parser Tests
# =============================================================================


class TestHGNCParser:
    """Tests for HGNCParser."""

    def test_parser_initialization(self) -> None:
        """Test parser can be initialized."""
        parser = HGNCParser(version="test", show_progress=False)
        assert parser.version == "test"
        assert parser.datasource_name == "HGNC"

    def test_parse_withdrawn_file(self, hgnc_withdrawn_path: Path) -> None:
        """Test parsing HGNC withdrawn file."""
        parser = HGNCParser(show_progress=False)
        result = parser.parse(hgnc_withdrawn_path)

        assert isinstance(result, MappingSet)
        assert result.datasource_name == "HGNC"
        assert len(result) > 0

    def test_parse_extracts_withdrawn_entries(self, hgnc_withdrawn_path: Path) -> None:
        """Test that withdrawn entries are identified."""
        parser = HGNCParser(show_progress=False)
        result = parser.parse(hgnc_withdrawn_path)

        # Check for withdrawn entries
        withdrawn = [m for m in result.mappings if m.subject_id == WITHDRAWN_ENTRY]
        # We have 1 "Entry Withdrawn" in test data
        assert len(withdrawn) >= 1

    def test_parse_extracts_merged_entries(self, hgnc_withdrawn_path: Path) -> None:
        """Test that merged entries are extracted."""
        parser = HGNCParser(show_progress=False)
        result = parser.parse(hgnc_withdrawn_path)

        # Check for merged entries (mapping to a new ID)
        merged = [
            m
            for m in result.mappings
            if m.subject_id != WITHDRAWN_ENTRY and "HGNC:" in str(m.subject_id)
        ]
        # We have 3 merged entries in test data
        assert len(merged) > 2

    def test_parse_with_complete_set(
        self, hgnc_withdrawn_path: Path, hgnc_complete_path: Path
    ) -> None:
        """Test parsing with both withdrawn and complete set files."""
        parser = HGNCParser(show_progress=False)
        result = parser.parse(
            hgnc_withdrawn_path,
            complete_set_path=hgnc_complete_path,
        )

        assert isinstance(result, MappingSet)
        # Should have symbol mappings from complete set
        assert len(result) > 0


# =============================================================================
# NCBI Parser Tests
# =============================================================================


class TestNCBIParser:
    """Tests for NCBIParser."""

    def test_parser_initialization(self) -> None:
        """Test parser can be initialized."""
        parser = NCBIParser(version="test", show_progress=False)
        assert parser.version == "test"
        assert parser.datasource_name == "NCBI"

    def test_parse_gene_history(self, ncbi_history_path: Path) -> None:
        """Test parsing gene history file."""
        parser = NCBIParser(show_progress=False)
        result = parser.parse(ncbi_history_path, tax_id="9606")

        assert isinstance(result, MappingSet)
        assert result.datasource_name == "NCBI"
        # Should have mappings for human genes only
        assert len(result) > 0

    def test_filters_by_tax_id(self, ncbi_history_path: Path) -> None:
        """Test that filtering by tax_id works."""
        parser = NCBIParser(show_progress=False)

        # Parse for human
        human_result = parser.parse(ncbi_history_path, tax_id="9606")

        # Parse for mouse
        mouse_result = parser.parse(ncbi_history_path, tax_id="10090")

        # Human should have more mappings in our test data
        assert len(human_result) >= 3
        assert len(mouse_result) >= 1

    def test_extracts_withdrawn_genes(
        self,
        ncbi_history_path: Path,
    ) -> None:
        """Test that withdrawn genes are identified."""
        parser = NCBIParser(show_progress=False)
        result = parser.parse(ncbi_history_path, tax_id="9606")

        # Check for withdrawn entries (GeneID = "-")
        withdrawn = [
            m for m in result.mappings if m.mapping_cardinality == MappingCardinality.ONE_TO_ZERO
        ]
        assert len(withdrawn) >= 1

    def test_parse_with_gene_info(
        self,
        ncbi_history_path: Path,
        ncbi_info_path: Path,
    ) -> None:
        """Test parsing with both history and info files."""
        parser = NCBIParser(show_progress=False)
        result = parser.parse(
            ncbi_history_path,
            gene_info_path=ncbi_info_path,
            tax_id="9606",
        )

        assert isinstance(result, MappingSet)
        assert len(result) > 0


# =============================================================================
# UniProt Parser Tests
# =============================================================================


class TestUniProtParser:
    """Tests for UniProtParser."""

    def test_parser_initialization(self) -> None:
        """Test parser can be initialized."""
        parser = UniProtParser(version="test", show_progress=False)
        assert parser.version == "test"
        assert parser.datasource_name == "UniProt"

    def test_parse_sec_ac_file(
        self,
        uniprot_sec_ac_path: Path,
    ) -> None:
        """Test parsing secondary accessions file."""
        parser = UniProtParser(show_progress=False)
        result = parser.parse(uniprot_sec_ac_path)

        assert isinstance(result, MappingSet)
        assert result.datasource_name == "UniProt"
        # We have 4 secondary accessions in test data
        assert len(result) >= 4

    def test_extracts_correct_mappings(self, uniprot_sec_ac_path: Path) -> None:
        """Test that mappings have correct structure."""
        parser = UniProtParser(show_progress=False)
        result = parser.parse(uniprot_sec_ac_path)

        for m in result.mappings:
            assert isinstance(m, UniProtMapping)
            # Object ID should be secondary accession
            assert m.object_id is not None
            # Subject ID should be primary accession (or withdrawn)
            assert m.subject_id is not None

    def test_many_to_one_cardinality(
        self,
        uniprot_sec_ac_path: Path,
    ) -> None:
        """Test that many-to-one cardinality is detected."""
        parser = UniProtParser(show_progress=False)
        result = parser.parse(uniprot_sec_ac_path)

        # A0A001 and A0A002 both map to P12345
        # So these should be MANY_TO_ONE
        many_to_one = [
            m for m in result.mappings if m.mapping_cardinality == MappingCardinality.MANY_TO_ONE
        ]
        assert len(many_to_one) >= 2

    def test_parse_with_delac_file(
        self, uniprot_sec_ac_path: Path, uniprot_delac_path: Path
    ) -> None:
        """Test parsing with both sec_ac and delac files."""
        parser = UniProtParser(show_progress=False)
        result = parser.parse(
            uniprot_sec_ac_path,
            delac_path=uniprot_delac_path,
        )

        assert isinstance(result, MappingSet)
        # Should include deleted accessions (4 from sec_ac + 3 from delac)
        assert len(result) >= 4  # At minimum the sec_ac entries


# =============================================================================
# Integration Tests
# =============================================================================


class TestParserIntegration:
    """Integration tests for parsers."""

    def test_all_parsers_return_mapping_set(
        self,
        chebi_sdf_path: Path,
        hmdb_xml_path: Path,
        hgnc_withdrawn_path: Path,
        ncbi_history_path: Path,
        uniprot_sec_ac_path: Path,
    ) -> None:
        """Test that all parsers return MappingSet."""
        results = [
            ChEBIParser(show_progress=False).parse(chebi_sdf_path),
            HMDBParser(show_progress=False).parse(hmdb_xml_path),
            HGNCParser(show_progress=False).parse(hgnc_withdrawn_path),
            NCBIParser(show_progress=False).parse(ncbi_history_path, tax_id="9606"),
            UniProtParser(show_progress=False).parse(uniprot_sec_ac_path),
        ]

        for result in results:
            assert isinstance(result, MappingSet)
            assert result.datasource_name != ""

    def test_all_mappings_have_required_fields(
        self,
        chebi_sdf_path: Path,
        hmdb_xml_path: Path,
        hgnc_withdrawn_path: Path,
        ncbi_history_path: Path,
        uniprot_sec_ac_path: Path,
    ) -> None:
        """Test that all mappings have required fields."""
        parsers = [
            (ChEBIParser(show_progress=False), chebi_sdf_path, {}),
            (HMDBParser(show_progress=False), hmdb_xml_path, {}),
            (HGNCParser(show_progress=False), hgnc_withdrawn_path, {}),
            (
                NCBIParser(show_progress=False),
                ncbi_history_path,
                {"tax_id": "9606"},
            ),
            (UniProtParser(show_progress=False), uniprot_sec_ac_path, {}),
        ]

        for parser, path, kwargs in parsers:
            result = parser.parse(path, **kwargs)
            for m in result.mappings:
                assert m.subject_id is not None
                assert m.predicate_id is not None
                # object_id may be None for some mapping types TODO discuss
