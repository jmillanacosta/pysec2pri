"""Tests for pysec2pri.parsers module."""

from pathlib import Path

import pytest
from sssom_schema import Mapping

from pysec2pri.parsers import (
    ChEBIParser,
    HGNCParser,
    HMDBMetaboliteParser,
    HMDBProteinParser,
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


class TestChEBIDistributionEras:
    """ChEBI's SDF/TSV format switch is resolved via DatasourceConfig.era_for."""

    def test_legacy_version_resolves_sdf_url(self) -> None:
        """A pre-245 release resolves to the legacy SDF era, not a missing config."""
        from pysec2pri.parsers.chebi import ChEBIDownloader

        downloader = ChEBIDownloader(subset="3star")
        urls = downloader.get_download_urls("100")
        assert urls.keys() == {"sdf"}
        assert "chebi_legacy" in urls["sdf"]
        assert "rel100" in urls["sdf"]
        assert downloader.get_format("100") == "sdf"

    def test_legacy_version_complete_subset_resolves_sdf_url(self) -> None:
        """The 'complete' subset picks the complete (non-3star) SDF file."""
        from pysec2pri.parsers.chebi import ChEBIDownloader

        urls = ChEBIDownloader(subset="complete").get_download_urls("100")
        assert urls["sdf"].endswith("ChEBI_complete.sdf.gz")

    def test_new_version_resolves_tsv_urls(self) -> None:
        """A >=245 release still resolves to the TSV flat-file URLs."""
        from pysec2pri.parsers.chebi import ChEBIDownloader

        downloader = ChEBIDownloader(subset="3star")
        urls = downloader.get_download_urls("245")
        assert urls.keys() == {"secondary_ids", "names", "compounds"}
        assert "rel245" in urls["secondary_ids"]
        assert downloader.get_format("245") == "tsv"

    def test_era_for_boundary_versions(self) -> None:
        """244 is the last SDF release, 245 is the first TSV release."""
        from pysec2pri.parsers.base import get_datasource_config

        config = get_datasource_config("chebi")
        sdf_era = config.era_for("244")
        tsv_era = config.era_for("245")
        assert sdf_era is not None
        assert tsv_era is not None
        assert sdf_era.id == "sdf"
        assert tsv_era.id == "tsv"

    def test_era_for_no_version_returns_none(self) -> None:
        """No version given means no era match; callers fall back to top-level config."""
        from pysec2pri.parsers.base import get_datasource_config

        config = get_datasource_config("chebi")
        assert config.era_for(None) is None


class TestHMDBParsers:
    """Tests for HMDB parsers."""

    def test_parse_met_returns_mapping_set(self, hmdb_xml_path: Path) -> None:
        """parse() returns a Sec2PriMappingSet."""
        parser = HMDBMetaboliteParser(show_progress=False)
        result = parser.parse(hmdb_xml_path)
        assert isinstance(result, Sec2PriMappingSet)

    def test_parse_produces_mappings(self, hmdb_xml_path: Path) -> None:
        """parse() extracts sec->pri mappings (fake file has 3)."""
        parser = HMDBMetaboliteParser(show_progress=False)
        result = parser.parse(hmdb_xml_path)
        assert len(result.mappings or []) == 3

    def test_parse_subject_object_ids(self, hmdb_xml_path: Path) -> None:
        """Secondary IDs are subjects; primary IDs are objects."""
        parser = HMDBMetaboliteParser(show_progress=False)
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
        parser = HMDBMetaboliteParser(show_progress=False)
        result = parser.parse(hmdb_xml_path)
        subjects = {m.subject_id for m in (result.mappings or [])}
        # HMDB0000003 has no secondary accessions
        assert "HMDB:HMDB0000003" not in subjects

    def test_parse_prot_returns_mapping_set(self, hmdb_proteins_xml_path: Path) -> None:
        """parse() returns a Sec2PriMappingSet."""
        parser = HMDBProteinParser(show_progress=False)
        result = parser.parse(hmdb_proteins_xml_path)
        assert isinstance(result, Sec2PriMappingSet)

    def test_parse_prot_produces_mappings(self, hmdb_proteins_xml_path: Path) -> None:
        """parse() extracts sec->pri mappings (fake has 3)."""
        parser = HMDBProteinParser(show_progress=False)
        result = parser.parse(hmdb_proteins_xml_path)
        assert len(result.mappings or []) == 3

    def test_parse_bare_number_normalised(self, hmdb_proteins_xml_path: Path) -> None:
        """Bare numeric secondary accessions are normalised to HMDBP prefix."""
        parser = HMDBProteinParser(show_progress=False)
        result = parser.parse(hmdb_proteins_xml_path)
        subjects = {m.subject_id for m in (result.mappings or [])}
        # "5229" should become "HMDB:HMDBP05229" (secondary to subject)
        assert "HMDBP:HMDBP05229" in subjects
        # full HMDBP accession unchanged
        assert "HMDBP:HMDBP05261" in subjects

    def test_parse_prot_no_secondary_skipped(self, hmdb_proteins_xml_path: Path) -> None:
        """Protein records with empty secondary_accessions produce no mappings."""
        parser = HMDBProteinParser(show_progress=False)
        result = parser.parse(hmdb_proteins_xml_path)
        subjects = {m.subject_id for m in (result.mappings or [])}
        # HMDBP00003 has no secondary accessions
        assert "HMDBP:HMDBP00003" not in subjects

    def test_parse_populates_primary_ids_including_no_secondary(self, hmdb_xml_path: Path) -> None:
        """_primary_ids includes ALL metabolites, even those with no secondaries."""
        parser = HMDBMetaboliteParser(show_progress=False)
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

    def test_parse_labels(self, hgnc_complete_path: Path) -> None:
        """Test parsing output for HGNC labels."""
        parser = HGNCParser(show_progress=False)
        result = parser.parse_labels(hgnc_complete_path)
        assert isinstance(result, Sec2PriMappingSet)

    def test_name2synonym_excludes_previous_labels(self, hgnc_complete_path: Path) -> None:
        """to_name2synonym must only contain hasExactSynonym rows (alias labels).

        Previous labels carry IAO:0100001 and are deprecation mappings, not
        synonyms.  They must appear in to_label_sec2pri but never in
        to_name2synonym.

        mock_hgnc_complete.tsv:
          BRCA1  alias: BRCC1, RNF53   prev: PSCP
          TP53   alias: P53, LFS1      prev: tumor_p53
          MYC    alias: c-Myc          prev: v-myc
        """
        parser = HGNCParser(show_progress=False)
        result = parser.parse_labels(hgnc_complete_path)

        df = result.to_name2synonym()
        synonyms = set(df["synonym"].tolist())
        # Aliases (hasExactSynonym) must appear in name2synonym
        assert "BRCC1" in synonyms
        assert "RNF53" in synonyms
        assert "P53" in synonyms
        assert "LFS1" in synonyms
        assert "c-Myc" in synonyms
        # Previous labels (IAO:0100001) must NOT appear in name2synonym
        assert "PSCP" not in synonyms
        assert "tumor_p53" not in synonyms
        assert "v-myc" not in synonyms

    def test_label_sec2pri_contains_previous_labels(self, hgnc_complete_path: Path) -> None:
        """to_label_sec2pri contains all label mappings including previous labels."""
        parser = HGNCParser(show_progress=False)
        result = parser.parse_labels(hgnc_complete_path)

        df = result.to_label_sec2pri()
        # Column is now "secondary_label" (was subject_label)
        all_secondary = set(df["secondary_label"].tolist())
        # Both aliases and previous labels must appear in the full label mapping
        assert "BRCC1" in all_secondary
        assert "PSCP" in all_secondary
        assert "tumor_p53" in all_secondary
        assert "v-myc" in all_secondary


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

        """Test parsing output for NCBI labels."""
        parser = NCBIParser(show_progress=False)
        result = parser.parse_labels(ncbi_info_path, tax_id="9606")
        assert isinstance(result, Sec2PriMappingSet)

    def test_parse_labels_direction(self, ncbi_info_path: Path) -> None:
        """Label mappings must be sec:pri: synonym is subject, primary label is object."""
        parser = NCBIParser(show_progress=False)
        result = parser.parse_labels(ncbi_info_path, tax_id="9606")
        mappings = result.mappings or []
        assert len(mappings) > 0
        subject_labels = {m.subject_label for m in mappings}
        object_labels = {m.object_label for m in mappings}
        # mock_gene_info.tsv: GENE1 has synonyms ALT_SYM1 and ALT_SYM2
        # synonyms must be subjects; current label must be the object
        assert "ALT_SYM1" in subject_labels
        assert "ALT_SYM2" in subject_labels
        assert "GENE1" in object_labels
        # the current primary label must never appear as a subject_label
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
            HMDBMetaboliteParser(show_progress=False).parse(hmdb_xml_path),
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
            HMDBMetaboliteParser(show_progress=False).parse(hmdb_xml_path),
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


class TestMappingDate:
    """Tests for release-date driven SSSOM ``mapping_date`` resolution."""

    def test_resolve_mapping_date_prefers_release_datetime(self) -> None:
        """A datetime release date is rendered as an ISO date string."""
        from datetime import datetime

        parser = HGNCParser(show_progress=False)
        parser.release_date = datetime(2021, 4, 1, 13, 30, 0)
        assert parser._resolve_mapping_date() == "2021-04-01"

    def test_resolve_mapping_date_accepts_date_and_str(self) -> None:
        """``date`` objects and ISO strings pass through unchanged."""
        from datetime import date

        parser = HGNCParser(show_progress=False)
        parser.release_date = date(2022, 7, 15)
        assert parser._resolve_mapping_date() == "2022-07-15"

        parser.release_date = "2019-12-31"
        assert parser._resolve_mapping_date() == "2019-12-31"

    def test_resolve_mapping_date_falls_back_to_iso_version(self) -> None:
        """When no release date is set, an ISO version (e.g. HGNC) is used."""
        parser = HGNCParser(version="2020-10-01", show_progress=False)
        assert parser.release_date is None
        assert parser._resolve_mapping_date() == "2020-10-01"

    def test_resolve_mapping_date_defaults_to_today(self) -> None:
        """A non-date version with no release date falls back to today."""
        from datetime import date

        parser = ChEBIParser(version="245", show_progress=False)
        assert parser._resolve_mapping_date() == date.today().isoformat()

    def test_mapping_set_uses_release_date(self, hgnc_withdrawn_path: Path) -> None:
        """The generated mapping set's ``mapping_date`` is the release date."""
        from datetime import datetime

        parser = HGNCParser(show_progress=False)
        parser.release_date = datetime(2021, 4, 1)
        result = parser.parse(hgnc_withdrawn_path)
        assert result.mapping_date == "2021-04-01"


class TestPerMappingDate:
    """Tests for per-Mapping ``mapping_date`` sourced from upstream record dates."""

    def test_hgnc_prev_symbol_uses_date_symbol_changed(self, hgnc_complete_path: Path) -> None:
        """Previous-symbol mappings carry HGNC's own ``date_symbol_changed``.

        mock_hgnc_complete.tsv: BRCA1 prev_symbol=PSCP, date_symbol_changed=1996-03-01.
        """
        parser = HGNCParser(show_progress=False)
        result = parser.parse_labels(hgnc_complete_path)
        prev_mappings = [m for m in result.mappings if m.subject_label == "PSCP"]
        assert len(prev_mappings) == 1
        assert prev_mappings[0].mapping_date == "1996-03-01"

    def test_hgnc_alias_symbol_has_no_mapping_date(self, hgnc_complete_path: Path) -> None:
        """Alias-symbol mappings have no associated change date in HGNC's data."""
        parser = HGNCParser(show_progress=False)
        result = parser.parse_labels(hgnc_complete_path)
        alias_mappings = [m for m in result.mappings if m.subject_label == "BRCC1"]
        assert len(alias_mappings) == 1
        assert alias_mappings[0].mapping_date is None

    def test_ncbi_discontinued_uses_discontinue_date(self, ncbi_history_path: Path) -> None:
        """Discontinued-ID mappings carry NCBI's own ``Discontinue_Date``, as ISO.

        mock_gene_history.tsv: GeneID=1001 -> Discontinued_GeneID=100001,
        Discontinue_Date=20230115.
        """
        parser = NCBIParser(show_progress=False)
        result = parser.parse(ncbi_history_path, tax_id="9606")
        matches = [m for m in result.mappings if m.subject_id == "NCBIGene:100001"]
        assert len(matches) == 1
        assert matches[0].mapping_date == "2023-01-15"

    def test_chebi_id_mapping_uses_consolidated_date(
        self, chebi_sdf_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """ChEBI ID mappings carry the consolidated first-seen date, when known.

        mock_chebi.sdf: CHEBI:10001 has secondary IDs CHEBI:99901, CHEBI:99902.
        Only the first is faked into the consolidated cache; the second must
        fall back to no per-row date (set-level pinned date applies instead).
        """
        parser = ChEBIParser(show_progress=False)
        m_meta = parser.get_mapping_metadata()
        record_id = parser._record_id(str(m_meta["record_id"]), "CHEBI:10001", "CHEBI:99901")
        monkeypatch.setattr(
            "pysec2pri.consolidate.load_mapping_dates",
            lambda *args, **kwargs: {record_id: "2013-02-15"},
        )

        result = parser.parse(chebi_sdf_path)

        known = [m for m in result.mappings if m.subject_id == "CHEBI:99901"]
        assert len(known) == 1
        assert known[0].mapping_date == "2013-02-15"

        unknown = [m for m in result.mappings if m.subject_id == "CHEBI:99902"]
        assert len(unknown) == 1
        assert unknown[0].mapping_date is None

    def test_uniprot_id_mapping_uses_consolidated_date(
        self, uniprot_sec_ac_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """UniProt ID mappings carry the consolidated first-seen date, when known.

        mock_sec_ac.txt: secondary A0A009QSN8 -> primary A0ACD6B9T2.
        """
        parser = UniProtParser(show_progress=False)
        m_meta = parser.get_mapping_metadata()
        record_id = parser._record_id(
            str(m_meta["record_id"]), "UniProtKB:A0ACD6B9T2", "UniProtKB:A0A009QSN8"
        )
        monkeypatch.setattr(
            "pysec2pri.consolidate.load_mapping_dates",
            lambda *args, **kwargs: {record_id: "2015-07-01"},
        )

        result = parser.parse(uniprot_sec_ac_path)

        known = [m for m in result.mappings if m.subject_id == "UniProtKB:A0A009QSN8"]
        assert len(known) == 1
        assert known[0].mapping_date == "2015-07-01"
