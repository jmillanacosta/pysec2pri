"""Tests for pysec2pri.parsers module."""

from pathlib import Path

import pytest
from sssom_schema import Mapping

from pysec2pri.parsers import (
    ChEBIParser,
    EnsemblParser,
    HGNCParser,
    HMDBMetaboliteParser,
    HMDBProteinParser,
    NCBIParser,
    UniProtParser,
    VGNCParser,
)
from pysec2pri.parsers.base import BaseMappingSet

TEST_DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def chebi_sdf_path() -> Path:
    """Return the fake ChEBI SDF test file."""
    return TEST_DATA_DIR / "mock_chebi.sdf"


@pytest.fixture
def chebi_names_tsv_path() -> Path:
    """Return the fake ChEBI names.tsv test file."""
    return TEST_DATA_DIR / "mock_chebi_names.tsv"


@pytest.fixture
def chebi_compounds_tsv_path() -> Path:
    """Return the fake ChEBI compounds.tsv test file."""
    return TEST_DATA_DIR / "mock_chebi_compounds.tsv"


@pytest.fixture
def hmdb_xml_path() -> Path:
    """Return the fake HMDB metabolites XML test file."""
    return TEST_DATA_DIR / "mock_hmdb.xml"


@pytest.fixture
def hmdb_proteins_xml_path() -> Path:
    """Return the fake HMDB proteins XML test file."""
    return TEST_DATA_DIR / "mock_hmdb_proteins.xml"


@pytest.fixture
def hgnc_withdrawn_path() -> Path:
    """Return the fake HGNC withdrawn.txt test file."""
    return TEST_DATA_DIR / "mock_hgnc_withdrawn.tsv"


@pytest.fixture
def hgnc_complete_path() -> Path:
    """Return the fake HGNC complete-set test file."""
    return TEST_DATA_DIR / "mock_hgnc_complete.tsv"


@pytest.fixture
def ncbi_history_path() -> Path:
    """Return the fake NCBI gene_history test file."""
    return TEST_DATA_DIR / "mock_gene_history.tsv"


@pytest.fixture
def ncbi_info_path() -> Path:
    """Return the fake NCBI gene_info test file."""
    return TEST_DATA_DIR / "mock_gene_info.tsv"


@pytest.fixture
def uniprot_sec_ac_path() -> Path:
    """Return the fake UniProt sec_ac.txt test file."""
    return TEST_DATA_DIR / "mock_sec_ac.txt"


@pytest.fixture
def vgnc_withdrawn_path() -> Path:
    """Return the fake VGNC withdrawn test file."""
    return TEST_DATA_DIR / "mock_vgnc_withdrawn.tsv"


@pytest.fixture
def vgnc_gene_set_path() -> Path:
    """Return the fake VGNC gene-set test file."""
    return TEST_DATA_DIR / "mock_vgnc_gene_set.tsv"


@pytest.fixture
def ensembl_stable_id_event_path() -> Path:
    """Return the fake Ensembl stable_id_event.txt test file (real-data slice)."""
    return TEST_DATA_DIR / "mock_ensembl_stable_id_event.txt"


@pytest.fixture
def ensembl_mapping_session_path() -> Path:
    """Return the fake Ensembl mapping_session.txt test file (real-data slice)."""
    return TEST_DATA_DIR / "mock_ensembl_mapping_session.txt"


@pytest.fixture
def ensembl_gene_path() -> Path:
    """Return the fake Ensembl gene.txt test file (real-data slice)."""
    return TEST_DATA_DIR / "mock_ensembl_gene.txt"


@pytest.fixture
def ensembl_xref_path() -> Path:
    """Return the fake Ensembl xref.txt test file (real-data slice)."""
    return TEST_DATA_DIR / "mock_ensembl_xref.txt"


@pytest.fixture
def ensembl_external_synonym_path() -> Path:
    """Return the fake Ensembl external_synonym.txt test file (real-data slice)."""
    return TEST_DATA_DIR / "mock_ensembl_external_synonym.txt"


class TestChEBIParser:
    """Tests for the ChEBI parser."""

    def test_parse(self, chebi_sdf_path: Path) -> None:
        """parse() returns a populated mapping set."""
        result = ChEBIParser(show_progress=False).parse(chebi_sdf_path)
        assert isinstance(result, BaseMappingSet)
        assert len(result.mappings) > 0

    def test_parse_synonyms_direction(self, chebi_sdf_path: Path) -> None:
        """Synonym mappings are sec:pri: synonym is subject, primary name is object."""
        result = ChEBIParser(show_progress=False).parse_synonyms(chebi_sdf_path)
        mappings = result.mappings or []
        assert len(mappings) > 0
        subject_labels = {m.subject_label for m in mappings}
        object_labels = {m.object_label for m in mappings}
        # mock_chebi.sdf: CHEBI:10001 has primary name "alpha-glucose" and
        # synonyms "glucose", "D-glucose": synonyms must be subjects, the
        # primary/canonical name must never appear as a subject_label.
        assert {"glucose", "D-glucose"} <= subject_labels
        assert "alpha-glucose" in object_labels
        assert "alpha-glucose" not in subject_labels

    def test_parse_synonyms_tsv_direction(
        self, chebi_names_tsv_path: Path, chebi_compounds_tsv_path: Path
    ) -> None:
        """TSV path: the primary name comes from compounds.tsv, not a names.tsv heuristic."""
        result = ChEBIParser(show_progress=False, subset="3star").parse_synonyms(
            names_path=chebi_names_tsv_path, compounds_path=chebi_compounds_tsv_path
        )
        mappings = result.mappings or []
        subject_labels = {m.subject_label for m in mappings}
        object_labels = {m.object_label for m in mappings}
        # "alpha-glucose" is the compounds.tsv canonical name for CHEBI:10001:
        # it must be object_label only; its names.tsv self-entry is excluded.
        assert "alpha-glucose" in object_labels
        assert "alpha-glucose" not in subject_labels
        assert {"glucose", "D-glucose"} <= subject_labels


class TestChEBIDistributionEras:
    """ChEBI's SDF/TSV format switch is resolved via DatasourceConfig.era_for."""

    def test_legacy_version_resolves_sdf_url(self) -> None:
        """A pre-245 release resolves to the legacy SDF era for both subsets."""
        from pysec2pri.downloads.chebi import ChEBIDownloader

        urls = ChEBIDownloader(subset="3star").get_download_urls("100")
        assert urls.keys() == {"sdf"}
        assert "chebi_legacy" in urls["sdf"] and "rel100" in urls["sdf"]
        assert ChEBIDownloader(subset="3star").get_format("100") == "sdf"
        assert (
            ChEBIDownloader(subset="complete")
            .get_download_urls("100")["sdf"]
            .endswith("ChEBI_complete.sdf.gz")
        )

    def test_new_version_resolves_tsv_urls(self) -> None:
        """A >=245 release resolves to the TSV flat-file URLs."""
        from pysec2pri.downloads.chebi import ChEBIDownloader

        downloader = ChEBIDownloader(subset="3star")
        urls = downloader.get_download_urls("245")
        assert urls.keys() == {"secondary_ids", "names", "compounds"}
        assert "rel245" in urls["secondary_ids"]
        assert downloader.get_format("245") == "tsv"

    def test_era_for_boundary_versions(self) -> None:
        """244 is the last SDF release, 245 the first TSV release; no version means no era."""
        from pysec2pri.parsers.base import get_datasource_config

        config = get_datasource_config("chebi", config_package="pysec2pri.config")
        sdf_era = config.era_for("244")
        tsv_era = config.era_for("245")
        assert sdf_era is not None and sdf_era.id == "sdf"
        assert tsv_era is not None and tsv_era.id == "tsv"
        assert config.era_for(None) is None


class TestHMDBParsers:
    """Tests for the HMDB metabolite and protein parsers."""

    def test_metabolites_parse(self, hmdb_xml_path: Path) -> None:
        """parse() extracts sec->pri mappings; records with no secondary are skipped."""
        result = HMDBMetaboliteParser(show_progress=False).parse(hmdb_xml_path)
        assert isinstance(result, BaseMappingSet)
        mappings = result.mappings or []
        assert len(mappings) == 3
        subjects = {m.subject_id for m in mappings}
        objects = {m.object_id for m in mappings}
        assert {"HMDB:HMDB00001", "HMDB:HMDB0001001", "HMDB:HMDB00002"} <= subjects
        assert {"HMDB:HMDB0000001", "HMDB:HMDB0000002"} <= objects
        # HMDB0000003 has no secondary accessions: no mapping, but still a primary.
        assert "HMDB:HMDB0000003" not in subjects
        assert "HMDB:HMDB0000003" in result.to_pri_ids()

    def test_proteins_parse(self, hmdb_proteins_xml_path: Path) -> None:
        """parse() normalises bare numeric accessions and skips no-secondary records."""
        result = HMDBProteinParser(show_progress=False).parse(hmdb_proteins_xml_path)
        assert isinstance(result, BaseMappingSet)
        mappings = result.mappings or []
        assert len(mappings) == 3
        subjects = {m.subject_id for m in mappings}
        # bare "5229" normalises to "HMDBP:HMDBP05229"; full accessions pass through.
        assert {"HMDBP:HMDBP05229", "HMDBP:HMDBP05261"} <= subjects
        # HMDBP00003 has no secondary accessions.
        assert "HMDBP:HMDBP00003" not in subjects


class TestHGNCParser:
    """Tests for the HGNC parser."""

    def test_parse(self, hgnc_withdrawn_path: Path) -> None:
        """parse() returns a populated mapping set for withdrawn IDs."""
        result = HGNCParser(show_progress=False).parse(hgnc_withdrawn_path)
        assert isinstance(result, BaseMappingSet)
        assert len(result.mappings) > 0

    def test_label_routing_separates_synonyms_from_previous_labels(
        self, hgnc_complete_path: Path
    ) -> None:
        """Alias labels (hasExactSynonym) and previous labels (IAO:0100001) route differently.

        mock_hgnc_complete.tsv: BRCA1 alias=BRCC1/RNF53, prev=PSCP; TP53 alias=P53, prev=tumor_p53.
        ``to_name2synonym`` must only ever contain aliases; ``to_label_sec2pri`` contains both.
        """
        result = HGNCParser(show_progress=False).parse_labels(hgnc_complete_path)
        synonyms = set(result.to_name2synonym()["synonym"])
        all_secondary = set(result.to_label_sec2pri()["secondary_label"])
        assert {"BRCC1", "RNF53", "P53"} <= synonyms
        assert "PSCP" not in synonyms and "tumor_p53" not in synonyms
        assert {"BRCC1", "PSCP", "tumor_p53"} <= all_secondary


class TestVGNCParser:
    """Tests for the VGNC parser."""

    def test_parse(self, vgnc_withdrawn_path: Path) -> None:
        """parse() returns a populated mapping set for withdrawn IDs, never species-filtered.

        mock_vgnc_withdrawn.tsv: VGNC:9001 withdrawn with no replacement;
        VGNC:9002 merged into VGNC:9501/CHIMPGENE1.
        """
        result = VGNCParser(show_progress=False).parse(vgnc_withdrawn_path)
        assert isinstance(result, BaseMappingSet)
        by_subject = {m.subject_id: m for m in result.mappings}
        assert by_subject["VGNC:9001"].object_id == "sssom:NoTermFound"
        assert by_subject["VGNC:9001"].predicate_id == "oboInOwl:consider"
        assert by_subject["VGNC:9002"].object_id == "VGNC:9501"
        assert by_subject["VGNC:9002"].predicate_id == "IAO:0100001"

    def test_parse_with_gene_set_populates_primary_ids_across_all_species(
        self, vgnc_withdrawn_path: Path, vgnc_gene_set_path: Path
    ) -> None:
        """complete_set_path populates _primary_ids with every species' IDs, unfiltered."""
        result = VGNCParser(show_progress=False).parse(
            vgnc_withdrawn_path, complete_set_path=vgnc_gene_set_path
        )
        pri_ids = result.to_pri_ids()
        # mock_vgnc_gene_set.tsv has chimp (9598) and dog (9615) entries: both present.
        assert {"VGNC:9501", "VGNC:9502", "VGNC:9503", "VGNC:9504"} <= set(pri_ids)

    def test_species_scoping_avoids_cross_species_symbol_collision(
        self, vgnc_gene_set_path: Path
    ) -> None:
        """Symbols are scoped per species: the same symbol in two species isn't ambiguous.

        mock_vgnc_gene_set.tsv: "SHAREDSYM" is the approved symbol for both
        VGNC:9503 (chimp, 9598) and VGNC:9504 (dog, 9615) -- orthologous
        genes sharing a name. Filtered per species, each view sees only its
        own ID for that symbol.
        """
        chimp = VGNCParser(show_progress=False).parse_primary_labels(
            vgnc_gene_set_path, species="9598"
        )
        dog = VGNCParser(show_progress=False).parse_primary_labels(
            vgnc_gene_set_path, species="9615"
        )
        assert chimp._primary_labels["SHAREDSYM"] == {"VGNC:9503"}
        assert dog._primary_labels["SHAREDSYM"] == {"VGNC:9504"}

    def test_parse_labels_direction_and_species_filter(self, vgnc_gene_set_path: Path) -> None:
        """Label mappings are sec:pri (synonym is subject); the other species' rows are excluded."""
        result = VGNCParser(show_progress=False).parse_labels(vgnc_gene_set_path, species="9598")
        mappings = result.mappings or []
        by_label = {m.subject_label: m for m in mappings}
        assert {"ALTC1", "OLDCHIMP"} <= set(by_label)
        assert "OLDDOG" not in by_label and "ALTD1" not in by_label
        assert by_label["ALTC1"].predicate_id == "oboInOwl:hasExactSynonym"
        assert by_label["OLDCHIMP"].predicate_id == "IAO:0100001"
        assert by_label["OLDCHIMP"].mapping_date == "2015-06-01"
        assert by_label["OLDCHIMP"].object_label == "CHIMPGENE1"

    def test_species_is_folded_into_mapping_set_id_and_record_id(
        self, vgnc_gene_set_path: Path
    ) -> None:
        """Different species at the same release get distinct mapping_set_id/record_id."""
        chimp = VGNCParser(version="2026-06-21", show_progress=False).parse_labels(
            vgnc_gene_set_path, species="9598"
        )
        dog = VGNCParser(version="2026-06-21", show_progress=False).parse_labels(
            vgnc_gene_set_path, species="9615"
        )
        assert chimp.mapping_set_id != dog.mapping_set_id

    def test_all_species_labels_combines_and_flags_shared_symbol_ambiguous(
        self, vgnc_gene_set_path: Path
    ) -> None:
        """species="all" processes every species together; a shared symbol becomes ambiguous.

        Unlike per-species scoping (see
        test_species_scoping_avoids_cross_species_symbol_collision), combining
        every species means "SHAREDSYM" genuinely has two distinct primary IDs
        with no other context to disambiguate them.
        """
        result = VGNCParser(show_progress=False).parse_labels(vgnc_gene_set_path, species="all")
        mappings = result.mappings or []
        by_label = {m.subject_label: m for m in mappings}
        # Both chimp and dog rows are present together now.
        assert {"ALTC1", "OLDCHIMP", "ALTD1", "OLDDOG"} <= set(by_label)
        assert result._primary_labels["SHAREDSYM"] == {"VGNC:9503", "VGNC:9504"}

    def test_ids_subset_by_species_drops_unresolvable_withdrawn_entries(
        self, vgnc_withdrawn_path: Path, vgnc_gene_set_path: Path
    ) -> None:
        """species= subsets withdrawn mappings by resolving the replacement gene's taxon.

        mock_vgnc_withdrawn.tsv: VGNC:9001 has no replacement (unresolvable,
        dropped under subsetting); VGNC:9002 merges into chimp's VGNC:9501.
        """
        full = VGNCParser(show_progress=False).parse(
            vgnc_withdrawn_path, complete_set_path=vgnc_gene_set_path
        )
        chimp = VGNCParser(show_progress=False).parse(
            vgnc_withdrawn_path, complete_set_path=vgnc_gene_set_path, species="9598"
        )
        dog = VGNCParser(show_progress=False).parse(
            vgnc_withdrawn_path, complete_set_path=vgnc_gene_set_path, species="9615"
        )
        assert {m.subject_id for m in full.mappings} == {"VGNC:9001", "VGNC:9002"}
        assert {m.subject_id for m in chimp.mappings} == {"VGNC:9002"}
        assert dog.mappings == []

    def test_ids_species_requires_complete_set_path(self, vgnc_withdrawn_path: Path) -> None:
        """species= without complete_set_path raises -- there's no way to resolve taxon IDs."""
        with pytest.raises(ValueError, match="complete_set_path"):
            VGNCParser(show_progress=False).parse(vgnc_withdrawn_path, species="9598")


class TestNCBIParser:
    """Tests for the NCBI parser."""

    def test_parse(self, ncbi_history_path: Path) -> None:
        """parse() returns a populated mapping set for gene_history."""
        result = NCBIParser(show_progress=False).parse(ncbi_history_path, species="9606")
        assert isinstance(result, BaseMappingSet)
        assert len(result.mappings) > 0

    def test_parse_all_species_includes_every_organism(self, ncbi_history_path: Path) -> None:
        """species="all" skips the taxon filter, including non-human rows.

        mock_gene_history.tsv has 4 human (9606) rows and 1 mouse (10090) row.
        """
        human_only = NCBIParser(show_progress=False).parse(ncbi_history_path, species="9606")
        all_species = NCBIParser(show_progress=False).parse(ncbi_history_path, species="all")
        assert len(all_species.mappings) > len(human_only.mappings)
        subjects = {m.subject_id for m in all_species.mappings}
        assert "NCBIGene:999999" in subjects  # the mouse row's Discontinued_GeneID

    def test_parse_with_gene_info_populates_primary_ids(
        self, ncbi_history_path: Path, ncbi_info_path: Path
    ) -> None:
        """Passing gene_info_path populates _primary_ids with all current IDs."""
        result = NCBIParser(show_progress=False).parse(
            ncbi_history_path, species="9606", gene_info_path=ncbi_info_path
        )
        pri_ids = result.to_pri_ids()
        assert len(pri_ids) > 0
        assert all(id_.startswith("NCBIGene:") for id_ in pri_ids)

    def test_parse_labels_direction(self, ncbi_info_path: Path) -> None:
        """Label mappings are sec:pri: synonym is subject, current label is object."""
        result = NCBIParser(show_progress=False).parse_labels(ncbi_info_path, species="9606")
        mappings = result.mappings or []
        subject_labels = {m.subject_label for m in mappings}
        object_labels = {m.object_label for m in mappings}
        # mock_gene_info.tsv: GENE1 has synonyms ALT_SYM1/ALT_SYM2.
        assert {"ALT_SYM1", "ALT_SYM2"} <= subject_labels
        assert "GENE1" in object_labels
        assert "GENE1" not in subject_labels


class TestEnsemblParser:
    """Tests for the Ensembl parser.

    Fixtures are real-data slices from Ensembl release 115
    (homo_sapiens_core_115_38), not hand-authored, so the schema/semantics
    match the live FTP dumps exactly.
    """

    def test_parse_skips_identity_and_non_gene_rows(
        self,
        ensembl_stable_id_event_path: Path,
        ensembl_mapping_session_path: Path,
    ) -> None:
        """Identity rows (old==new) and non-gene types are excluded."""
        result = EnsemblParser(version="115", show_progress=False).parse(
            ensembl_stable_id_event_path,
            mapping_session_path=ensembl_mapping_session_path,
        )
        # mock file has 6 rows: 1 identity (skip), 1 retirement, 1 scored
        # rename, 1 new-gene with no old id (skip), 1 plain rename, 1
        # transcript-type row (skip) -> exactly 3 mappings survive.
        assert len(result.mappings) == 3
        subject_ids = {m.subject_id for m in result.mappings}
        assert "ENSEMBL:ENSG00000000003" not in subject_ids  # identity row

    def test_parse_retirement_uses_consider_predicate(
        self,
        ensembl_stable_id_event_path: Path,
        ensembl_mapping_session_path: Path,
    ) -> None:
        """A row with new_stable_id = N becomes a withdrawn-entry mapping."""
        result = EnsemblParser(version="115", show_progress=False).parse(
            ensembl_stable_id_event_path,
            mapping_session_path=ensembl_mapping_session_path,
        )
        withdrawn = [m for m in result.mappings if m.subject_id == "ENSEMBL:ENSG00000000893"]
        assert len(withdrawn) == 1
        assert withdrawn[0].predicate_id == "oboInOwl:consider"
        assert withdrawn[0].object_id == "sssom:NoTermFound"
        assert withdrawn[0].mapping_date == "2002-09-05"

    def test_parse_rename_uses_replaced_by_predicate_and_confidence(
        self,
        ensembl_stable_id_event_path: Path,
        ensembl_mapping_session_path: Path,
    ) -> None:
        """A row with old != new becomes a replaced-by mapping with score as confidence."""
        result = EnsemblParser(version="115", show_progress=False).parse(
            ensembl_stable_id_event_path,
            mapping_session_path=ensembl_mapping_session_path,
        )
        scored = [m for m in result.mappings if m.subject_id == "ENSEMBL:ENSG00000007565"]
        assert len(scored) == 1
        assert scored[0].predicate_id == "IAO:0100001"
        assert scored[0].predicate_label == "term replaced by"
        assert scored[0].object_id == "ENSEMBL:ENSG00000206171"
        assert scored[0].confidence == pytest.approx(0.973291)
        assert scored[0].mapping_date == "2006-03-10"

    def test_parse_notes_an_assembly_change_as_a_comment(
        self,
        ensembl_stable_id_event_path: Path,
        ensembl_mapping_session_path: Path,
    ) -> None:
        """A rename spanning an assembly change is noted in comment, not source_version."""
        result = EnsemblParser(version="115", show_progress=False).parse(
            ensembl_stable_id_event_path,
            mapping_session_path=ensembl_mapping_session_path,
        )
        by_subject = {m.subject_id: m for m in result.mappings}
        # session 361: NCBI35 (old) -> NCBI36 (new): a real assembly change.
        scored = by_subject["ENSEMBL:ENSG00000007565"]
        assert scored.comment == "Assembly changed from NCBI35 to NCBI36."
        assert scored.subject_source_version is None
        assert scored.object_source_version is None
        # session 388: GRCh37 -> GRCh37 (assembly-patch rename within one build): no comment.
        patched = by_subject["ENSEMBL:ASMPATCHG00000000170"]
        assert patched.comment is None
        # session 347: a withdrawal has no object side to diff against.
        withdrawn = by_subject["ENSEMBL:ENSG00000000893"]
        assert withdrawn.comment is None

    def test_set_level_source_version_is_always_the_release(
        self,
        ensembl_stable_id_event_path: Path,
        ensembl_mapping_session_path: Path,
    ) -> None:
        """The set-level source-version is the analyzed release, same as every datasource."""
        result = EnsemblParser(version="115", show_progress=False).parse(
            ensembl_stable_id_event_path,
            mapping_session_path=ensembl_mapping_session_path,
        )
        # default species is 9606 (human), which has a configured GRCh38 build,
        # but the set-level field is still the release, not the build.
        assert result.subject_source_version == "115"
        assert result.object_source_version == "115"

    def test_set_level_source_version_is_the_release_for_uncurated_species_too(
        self,
        ensembl_stable_id_event_path: Path,
        ensembl_mapping_session_path: Path,
    ) -> None:
        """A species with no configured build still gets the release as source-version."""
        result = EnsemblParser(version="115", show_progress=False, species=99999).parse(
            ensembl_stable_id_event_path,
            mapping_session_path=ensembl_mapping_session_path,
        )
        assert result.subject_source_version == "115"
        assert result.object_source_version == "115"

    def test_parse_without_mapping_session_has_no_per_row_date(
        self, ensembl_stable_id_event_path: Path
    ) -> None:
        """mapping_session is optional; per-row mapping_date is simply unset without it."""
        result = EnsemblParser(version="115", show_progress=False).parse(
            ensembl_stable_id_event_path
        )
        assert len(result.mappings) == 3
        assert all(m.mapping_date is None for m in result.mappings)

    def test_parse_with_gene_path_populates_primary_ids(
        self,
        ensembl_stable_id_event_path: Path,
        ensembl_mapping_session_path: Path,
        ensembl_gene_path: Path,
    ) -> None:
        """Passing gene_path populates _primary_ids with current gene IDs."""
        result = EnsemblParser(version="115", show_progress=False).parse(
            ensembl_stable_id_event_path,
            mapping_session_path=ensembl_mapping_session_path,
            gene_path=ensembl_gene_path,
        )
        pri_ids = result.to_pri_ids()
        assert len(pri_ids) == 3
        assert all(id_.startswith("ENSEMBL:") for id_ in pri_ids)

    def test_parse_labels_direction(
        self,
        ensembl_gene_path: Path,
        ensembl_xref_path: Path,
        ensembl_external_synonym_path: Path,
    ) -> None:
        """Label mappings are sec:pri: synonym is subject, current label is object."""
        result = EnsemblParser(version="115", show_progress=False).parse_labels(
            ensembl_gene_path, ensembl_xref_path, ensembl_external_synonym_path
        )
        mappings = result.mappings or []
        subject_labels = {m.subject_label for m in mappings}
        object_labels = {m.object_label for m in mappings}
        # mock_ensembl_external_synonym.txt: MT-RNR1 has synonyms 12S/MOTS-C/MTRNR1.
        assert {"12S", "MOTS-C", "MTRNR1"} <= subject_labels
        assert "MT-RNR1" in object_labels
        assert "MT-RNR1" not in subject_labels
        assert all(m.predicate_id == "oboInOwl:hasExactSynonym" for m in mappings)

    def test_parse_primary_ids(self, ensembl_gene_path: Path) -> None:
        """parse_primary_ids() populates only _primary_ids, no mappings."""
        result = EnsemblParser(version="115", show_progress=False).parse_primary_ids(
            ensembl_gene_path
        )
        assert result.mappings == []
        assert len(result.to_pri_ids()) == 3

    def test_parse_primary_labels(self, ensembl_gene_path: Path, ensembl_xref_path: Path) -> None:
        """parse_primary_labels() populates only _primary_labels, no mappings."""
        result = EnsemblParser(version="115", show_progress=False).parse_primary_labels(
            ensembl_gene_path, ensembl_xref_path
        )
        assert result.mappings == []
        assert "MT-RNR1" in result._primary_labels

    def test_species_is_folded_into_mapping_set_id_and_record_id(
        self,
        ensembl_stable_id_event_path: Path,
        ensembl_mapping_session_path: Path,
    ) -> None:
        """Different species at the same release get distinct mapping_set_id/record_id."""
        human = EnsemblParser(version="115", show_progress=False, species=9606).parse(
            ensembl_stable_id_event_path, mapping_session_path=ensembl_mapping_session_path
        )
        mouse = EnsemblParser(version="115", show_progress=False, species=10090).parse(
            ensembl_stable_id_event_path, mapping_session_path=ensembl_mapping_session_path
        )
        assert human.mapping_set_id != mouse.mapping_set_id
        assert "9606" in human.mapping_set_id
        assert "10090" in mouse.mapping_set_id
        human_rid = human.mappings[0].record_id
        mouse_rid = mouse.mappings[0].record_id
        assert human_rid != mouse_rid
        # Both still end in the same version-independent pair hash, since
        # the same (pri, sec) pair was parsed in both runs.
        assert human_rid[-16:] == mouse_rid[-16:]

    def test_parse_label_history_uses_replaced_by_predicate(self) -> None:
        """parse_label_history() builds IAO:0100001 mappings from precomputed transitions."""
        parser = EnsemblParser(version="115", show_progress=False, species=9606)
        result = parser.parse_label_history(
            [("ENSG00000211459", "OLD-SYM", "MT-RNR1", "2020-01-01")]
        )
        assert len(result.mappings) == 1
        m = result.mappings[0]
        assert m.predicate_id == "IAO:0100001"
        assert m.subject_label == "OLD-SYM"
        assert m.object_label == "MT-RNR1"
        assert m.object_id == "ENSEMBL:ENSG00000211459"
        assert m.mapping_date == "2020-01-01"
        assert "consolidate" in result.mapping_set_id


class TestUniProtParser:
    """Tests for the UniProt parser."""

    def test_parse(self, uniprot_sec_ac_path: Path) -> None:
        """parse() returns a populated mapping set."""
        result = UniProtParser(show_progress=False).parse(uniprot_sec_ac_path)
        assert isinstance(result, BaseMappingSet)
        assert len(result.mappings) >= 4


class TestParserIntegration:
    """Cross-parser invariants that every parser must satisfy."""

    def test_mappings_are_well_formed_across_all_parsers(
        self,
        chebi_sdf_path: Path,
        hmdb_xml_path: Path,
        hgnc_withdrawn_path: Path,
        ncbi_history_path: Path,
        uniprot_sec_ac_path: Path,
        ensembl_stable_id_event_path: Path,
        ensembl_mapping_session_path: Path,
    ) -> None:
        """Every parser produces Mapping instances with a subject/predicate, and a correct label."""
        results = [
            ChEBIParser(show_progress=False).parse(chebi_sdf_path),
            HMDBMetaboliteParser(show_progress=False).parse(hmdb_xml_path),
            HGNCParser(show_progress=False).parse(hgnc_withdrawn_path),
            NCBIParser(show_progress=False).parse(ncbi_history_path, species="9606"),
            UniProtParser(show_progress=False).parse(uniprot_sec_ac_path),
            EnsemblParser(show_progress=False).parse(
                ensembl_stable_id_event_path, mapping_session_path=ensembl_mapping_session_path
            ),
        ]
        for result in results:
            assert isinstance(result, BaseMappingSet)
            for m in result.mappings:
                assert isinstance(m, Mapping)
                assert m.subject_id is not None
                assert m.predicate_id is not None
                if m.predicate_id == "IAO:0100001":
                    assert m.predicate_label == "term replaced by"


class TestRecordIdVersioning:
    """record_id is release-scoped (an OWL Axiom IRI); _pair_hash is the join key.

    record_id must vary across releases of the same (pri, sec) pair, so that
    unioning several releases' SSSOM/RDF into one triplestore never asserts
    contradictory axioms under one IRI. Its trailing 16 hex characters are
    always the version-independent pair hash.
    """

    def test_pair_hash_is_version_independent(self) -> None:
        """The same (pri, sec) pair hashes identically regardless of parser version."""
        older = NCBIParser(version="2020-01-01", show_progress=False)
        newer = NCBIParser(version="2024-01-01", show_progress=False)
        assert older._pair_hash("A", "B") == newer._pair_hash("A", "B")

    def test_record_id_differs_across_versions_but_pair_hash_does_not(self) -> None:
        """record_id differs per version; its trailing pair hash stays stable."""
        older = ChEBIParser(version="200", show_progress=False)
        newer = ChEBIParser(version="245", show_progress=False)
        rid_older = older._record_id(older._record_namespace(), "CHEBI:1", "CHEBI:2")
        rid_newer = newer._record_id(newer._record_namespace(), "CHEBI:1", "CHEBI:2")
        assert rid_older != rid_newer
        assert rid_older[-16:] == rid_newer[-16:] == older._pair_hash("CHEBI:1", "CHEBI:2")

    def test_record_namespace_includes_version_and_product_slug(self) -> None:
        """_record_namespace mirrors mapping_set_id's {base}/{version}/{slug} ordering."""
        parser = NCBIParser(version="2026-01-01", show_progress=False)
        parser.species = "9606"
        assert parser._record_namespace() == "sec2pri:ncbigene/2026-01-01/9606/"

    def test_record_namespace_without_product_slug_omits_trailing_segment(self) -> None:
        """Parsers that don't override _product_slug only fold in the version."""
        parser = HGNCParser(version="2026-04-07", show_progress=False)
        assert parser._record_namespace() == "sec2pri:hgnc/2026-04-07/"


class TestMappingDate:
    """Tests for release-date driven SSSOM ``mapping_date`` resolution."""

    def test_resolve_mapping_date_priority_order(self) -> None:
        """release_date (datetime/date/str) wins; then an ISO version; then today."""
        from datetime import date, datetime

        parser = HGNCParser(show_progress=False)
        parser.release_date = datetime(2021, 4, 1, 13, 30, 0)
        assert parser._resolve_mapping_date() == "2021-04-01"

        parser.release_date = date(2022, 7, 15)
        assert parser._resolve_mapping_date() == "2022-07-15"

        parser.release_date = "2019-12-31"
        assert parser._resolve_mapping_date() == "2019-12-31"

        no_release_date_parser = HGNCParser(version="2020-10-01", show_progress=False)
        assert no_release_date_parser._resolve_mapping_date() == "2020-10-01"

        non_date_version_parser = ChEBIParser(version="245", show_progress=False)
        assert non_date_version_parser._resolve_mapping_date() == date.today().isoformat()

    def test_mapping_set_uses_release_date(self, hgnc_withdrawn_path: Path) -> None:
        """The generated mapping set's ``mapping_date`` is the parser's release date."""
        from datetime import datetime

        parser = HGNCParser(show_progress=False)
        parser.release_date = datetime(2021, 4, 1)
        result = parser.parse(hgnc_withdrawn_path)
        assert result.mapping_date == "2021-04-01"


class TestPerMappingDate:
    """Per-source ``mapping_date`` sourced from each datasource's own record dates."""

    def test_hgnc_dates_only_apply_to_previous_symbols(self, hgnc_complete_path: Path) -> None:
        """HGNC's date_symbol_changed applies to previous-symbol rows, not alias rows.

        mock_hgnc_complete.tsv: BRCA1 prev_symbol=PSCP, date_symbol_changed=1996-03-01;
        alias BRCC1 has no associated change date.
        """
        result = HGNCParser(show_progress=False).parse_labels(hgnc_complete_path)
        by_label = {m.subject_label: m for m in result.mappings}
        assert by_label["PSCP"].mapping_date == "1996-03-01"
        assert by_label["BRCC1"].mapping_date is None

    def test_ncbi_discontinued_uses_discontinue_date(self, ncbi_history_path: Path) -> None:
        """Discontinued-ID mappings carry NCBI's own Discontinue_Date, as ISO.

        mock_gene_history.tsv: GeneID=1001 -> Discontinued_GeneID=100001, Discontinue_Date=20230115.
        """
        result = NCBIParser(show_progress=False).parse(ncbi_history_path, species="9606")
        matches = [m for m in result.mappings if m.subject_id == "NCBIGene:100001"]
        assert len(matches) == 1
        assert matches[0].mapping_date == "2023-01-15"

    def test_chebi_id_mapping_uses_consolidated_date(
        self, chebi_sdf_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """ChEBI mappings carry the consolidated first-seen date, when known, else None.

        mock_chebi.sdf: CHEBI:10001 has secondary IDs CHEBI:99901, CHEBI:99902.
        Only the first is faked into the consolidated cache.
        """
        parser = ChEBIParser(show_progress=False)
        pair_key = parser._pair_hash("CHEBI:10001", "CHEBI:99901")
        monkeypatch.setattr(
            "pysec2pri.consolidate.load_mapping_dates",
            lambda *args, **kwargs: {pair_key: "2013-02-15"},
        )

        result = parser.parse(chebi_sdf_path)
        by_subject = {m.subject_id: m for m in result.mappings}
        assert by_subject["CHEBI:99901"].mapping_date == "2013-02-15"
        assert by_subject["CHEBI:99902"].mapping_date is None
