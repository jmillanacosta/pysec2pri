"""Tests for pysec2pri.diff module.

Tests cover:
- MappingDiff dataclass
- diff_mapping_sets function
- diff_sssom_files function
- summarize_diff function
"""

from pathlib import Path

import polars as pl
import pytest

from pysec2pri.diff import (
    MappingDiff,
    diff_mapping_sets,
    diff_sssom_files,
    mapping_set_to_dataframe,
    summarize_diff,
)
from pysec2pri.models import (
    BaseMapping,
    ChEBIMapping,
    MappingSet,
)

# Path to test data directory
TEST_DATA_DIR = Path(__file__).parent / "data"


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def sssom_v1_path() -> Path:
    """Path to fake SSSOM v1 file."""
    return TEST_DATA_DIR / Path("fake_sssom_v1.tsv")


@pytest.fixture
def sssom_v2_path() -> Path:
    """Path to fake SSSOM v2 file."""
    return TEST_DATA_DIR / Path("fake_sssom_v2.tsv")


@pytest.fixture
def empty_mapping_set() -> MappingSet:
    """Create an empty MappingSet."""
    return MappingSet(datasource_name="Test", version="1.0")


@pytest.fixture
def simple_mapping_set_v1() -> MappingSet:
    """Create a simple MappingSet (version 1)."""
    ms = MappingSet(datasource_name="Test", version="1.0")
    ms.add_mappings(
        [
            BaseMapping(
                subject_id="PRIMARY1",
                predicate_id="IAO:0100001",
                object_id="SECONDARY1",
            ),
            BaseMapping(
                subject_id="PRIMARY2",
                predicate_id="IAO:0100001",
                object_id="SECONDARY2",
            ),
            BaseMapping(
                subject_id="PRIMARY3",
                predicate_id="IAO:0100001",
                object_id="SECONDARY3",
            ),
        ]
    )
    return ms


@pytest.fixture
def simple_mapping_set_v2() -> MappingSet:
    """Create a simple MappingSet (version 2) with changes."""
    ms = MappingSet(datasource_name="Test", version="2.0")
    ms.add_mappings(
        [
            # Unchanged
            BaseMapping(
                subject_id="PRIMARY1",
                predicate_id="IAO:0100001",
                object_id="SECONDARY1",
            ),
            # Changed: SECONDARY2 now maps to PRIMARY4 instead of PRIMARY2
            BaseMapping(
                subject_id="PRIMARY4",
                predicate_id="IAO:0100001",
                object_id="SECONDARY2",
            ),
            # New mapping (added)
            BaseMapping(
                subject_id="PRIMARY5",
                predicate_id="IAO:0100001",
                object_id="SECONDARY5",
            ),
            # SECONDARY3 -> PRIMARY3 removed (not in v2)
        ]
    )
    return ms


# =============================================================================
# Tests for MappingDiff dataclass
# =============================================================================


class TestMappingDiff:
    """Tests for MappingDiff dataclass."""

    def test_create_empty_diff(self) -> None:
        """Test creating an empty diff."""
        diff = MappingDiff(
            old_version="1.0",
            new_version="2.0",
            datasource="Test",
        )
        assert diff.old_version == "1.0"
        assert diff.new_version == "2.0"
        assert diff.datasource == "Test"
        assert diff.added_count == 0
        assert diff.removed_count == 0
        assert diff.changed_count == 0

    def test_diff_properties(self) -> None:
        """Test diff property accessors."""
        diff = MappingDiff(
            old_version="1.0",
            new_version="2.0",
            datasource="Test",
            added=pl.DataFrame({"subject_id": ["A"], "object_id": ["1"]}),
            removed=pl.DataFrame({"subject_id": ["B", "C"], "object_id": ["2", "3"]}),
            changed=pl.DataFrame(
                {
                    "object_id": ["4"],
                    "old_subject_id": ["D"],
                    "new_subject_id": ["E"],
                }
            ),
        )
        assert diff.added_count == 1
        assert diff.removed_count == 2
        assert diff.changed_count == 1
        assert diff.total_changes == 4
        assert diff.has_changes is True

    def test_has_changes_false_when_empty(self) -> None:
        """Test has_changes is False for empty diff."""
        diff = MappingDiff(
            old_version="1.0",
            new_version="2.0",
            datasource="Test",
        )
        assert diff.has_changes is False


# =============================================================================
# Tests for mapping_set_to_dataframe
# =============================================================================


class TestMappingSetToDataframe:
    """Tests for mapping_set_to_dataframe function."""

    def test_converts_empty_set(self, empty_mapping_set: MappingSet) -> None:
        """Test converting empty MappingSet."""
        df = mapping_set_to_dataframe(empty_mapping_set)
        assert isinstance(df, pl.DataFrame)
        assert len(df) == 0
        assert "subject_id" in df.columns
        assert "object_id" in df.columns

    def test_converts_mappings(self, simple_mapping_set_v1: MappingSet) -> None:
        """Test converting MappingSet with mappings."""
        df = mapping_set_to_dataframe(simple_mapping_set_v1)
        assert len(df) == 3
        assert "subject_id" in df.columns
        assert "object_id" in df.columns

        # Check values
        subject_ids = df["subject_id"].to_list()
        assert "PRIMARY1" in subject_ids
        assert "PRIMARY2" in subject_ids


# =============================================================================
# Tests for diff_mapping_sets
# =============================================================================


class TestDiffMappingSets:
    """Tests for diff_mapping_sets function."""

    def test_diff_identical_sets(self, simple_mapping_set_v1: MappingSet) -> None:
        """Test diffing identical sets returns no changes."""
        diff = diff_mapping_sets(simple_mapping_set_v1, simple_mapping_set_v1)
        assert diff.has_changes is False
        assert diff.added_count == 0
        assert diff.removed_count == 0
        assert diff.changed_count == 0

    def test_diff_empty_to_populated(
        self, empty_mapping_set: MappingSet, simple_mapping_set_v1: MappingSet
    ) -> None:
        """Test diffing empty set to populated set."""
        diff = diff_mapping_sets(empty_mapping_set, simple_mapping_set_v1)
        assert diff.added_count == 3
        assert diff.removed_count == 0

    def test_diff_populated_to_empty(
        self, simple_mapping_set_v1: MappingSet, empty_mapping_set: MappingSet
    ) -> None:
        """Test diffing populated set to empty set."""
        diff = diff_mapping_sets(simple_mapping_set_v1, empty_mapping_set)
        assert diff.added_count == 0
        assert diff.removed_count == 3

    def test_diff_detects_added(
        self, simple_mapping_set_v1: MappingSet, simple_mapping_set_v2: MappingSet
    ) -> None:
        """Test that added mappings are detected."""
        diff = diff_mapping_sets(simple_mapping_set_v1, simple_mapping_set_v2)

        # SECONDARY5 -> PRIMARY5 is new in v2
        assert diff.added_count >= 1

        # Check added contains the new mapping
        added_objects = diff.added["object_id"].to_list()
        assert "SECONDARY5" in added_objects

    def test_diff_detects_removed(
        self, simple_mapping_set_v1: MappingSet, simple_mapping_set_v2: MappingSet
    ) -> None:
        """Test that removed mappings are detected."""
        diff = diff_mapping_sets(simple_mapping_set_v1, simple_mapping_set_v2)

        # SECONDARY3 -> PRIMARY3 is removed in v2
        assert diff.removed_count >= 1

        removed_objects = diff.removed["object_id"].to_list()
        assert "SECONDARY3" in removed_objects

    def test_diff_detects_changed(
        self, simple_mapping_set_v1: MappingSet, simple_mapping_set_v2: MappingSet
    ) -> None:
        """Test that changed mappings are detected."""
        diff = diff_mapping_sets(simple_mapping_set_v1, simple_mapping_set_v2)

        # SECONDARY2 changed from PRIMARY2 to PRIMARY4
        assert diff.changed_count >= 1

        changed_objects = diff.changed["object_id"].to_list()
        assert "SECONDARY2" in changed_objects

    def test_diff_captures_versions(
        self, simple_mapping_set_v1: MappingSet, simple_mapping_set_v2: MappingSet
    ) -> None:
        """Test that versions are captured in diff."""
        diff = diff_mapping_sets(simple_mapping_set_v1, simple_mapping_set_v2)
        assert diff.old_version == "1.0"
        assert diff.new_version == "2.0"

    def test_diff_captures_datasource(
        self, simple_mapping_set_v1: MappingSet, simple_mapping_set_v2: MappingSet
    ) -> None:
        """Test that datasource is captured in diff."""
        diff = diff_mapping_sets(simple_mapping_set_v1, simple_mapping_set_v2)
        assert diff.datasource == "Test"


# =============================================================================
# Tests for diff_sssom_files
# =============================================================================


class TestDiffSssomFiles:
    """Tests for diff_sssom_files function."""

    def test_diff_sssom_files(self, sssom_v1_path: str, sssom_v2_path: str) -> None:
        """Test diffing SSSOM files."""
        diff = diff_sssom_files(
            sssom_v1_path,
            sssom_v2_path,
            datasource="ChEBI",
        )

        assert isinstance(diff, MappingDiff)
        assert diff.datasource == "ChEBI"
        assert diff.has_changes is True

    def test_diff_detects_file_changes(self, sssom_v1_path: str, sssom_v2_path: str) -> None:
        """Test that file changes are detected."""
        diff = diff_sssom_files(
            sssom_v1_path,
            sssom_v2_path,
            datasource="ChEBI",
        )

        # v2 has:
        # - CHEBI:99905 added
        # - CHEBI:99904 removed
        # - CHEBI:99903 changed (10002 -> 10004)

        assert diff.added_count >= 1
        assert diff.removed_count >= 1
        assert diff.changed_count >= 1

    def test_diff_same_file(self, sssom_v1_path: str) -> None:
        """Test diffing a file with itself."""
        diff = diff_sssom_files(
            sssom_v1_path,
            sssom_v1_path,
            datasource="Test",
        )
        assert diff.has_changes is False

    def test_extracts_version_from_metadata(self, sssom_v1_path: str, sssom_v2_path: str) -> None:
        """Test that version is extracted from SSSOM metadata."""
        diff = diff_sssom_files(
            sssom_v1_path,
            sssom_v2_path,
            datasource="ChEBI",
        )
        # Versions should be extracted from #mapping_set_version
        assert diff.old_version == "v1.0"
        assert diff.new_version == "v2.0"


# =============================================================================
# Tests for summarize_diff
# =============================================================================


class TestSummarizeDiff:
    """Tests for summarize_diff function."""

    def test_summarize_empty_diff(self) -> None:
        """Test summarizing empty diff."""
        diff = MappingDiff(
            old_version="1.0",
            new_version="2.0",
            datasource="Test",
        )
        summary = summarize_diff(diff)

        assert "Test" in summary
        assert "1.0" in summary
        assert "2.0" in summary
        assert "0" in summary  # Counts are 0

    def test_summarize_with_changes(
        self, simple_mapping_set_v1: MappingSet, simple_mapping_set_v2: MappingSet
    ) -> None:
        """Test summarizing diff with changes."""
        diff = diff_mapping_sets(simple_mapping_set_v1, simple_mapping_set_v2)
        summary = summarize_diff(diff)

        assert "Test" in summary
        assert "Added" in summary
        assert "Removed" in summary
        assert "Changed" in summary

    def test_summarize_shows_details_for_small_diffs(self) -> None:
        """Test that small diffs show detailed changes."""
        diff = MappingDiff(
            old_version="1.0",
            new_version="2.0",
            datasource="Test",
            added=pl.DataFrame(
                {
                    "subject_id": ["NEW_PRI"],
                    "object_id": ["NEW_SEC"],
                }
            ),
        )
        summary = summarize_diff(diff)

        # Small diffs should show details
        assert "NEW_SEC" in summary or "Added:" in summary

    def test_summarize_format(self) -> None:
        """Test summary format is readable."""
        diff = MappingDiff(
            old_version="1.0",
            new_version="2.0",
            datasource="ChEBI",
        )
        summary = summarize_diff(diff)

        # Check it's a multi-line string
        lines = summary.split("\n")
        assert len(lines) > 1

        # Check it contains expected sections
        assert any("Diff Summary" in line for line in lines)
        assert any("Old version" in line for line in lines)
        assert any("New version" in line for line in lines)


# =============================================================================
# Integration Tests
# =============================================================================


class TestDiffIntegration:
    """Integration tests for diff functionality."""

    def test_diff_workflow(
        self, simple_mapping_set_v1: MappingSet, simple_mapping_set_v2: MappingSet
    ) -> None:
        """Test complete diff workflow."""
        # Create diff
        diff = diff_mapping_sets(simple_mapping_set_v1, simple_mapping_set_v2)

        # Generate summary
        summary = summarize_diff(diff)

        # Verify complete workflow
        assert diff.has_changes
        assert len(summary) > 0
        assert diff.total_changes > 0

    def test_diff_with_chebi_mappings(self) -> None:
        """Test diff with ChEBI-specific mappings."""
        ms1 = MappingSet(datasource_name="ChEBI", version="2024-01")
        ms1.add_mappings(
            [
                ChEBIMapping(
                    subject_id="CHEBI:10001",
                    predicate_id="IAO:0100001",
                    object_id="CHEBI:99901",
                ),
            ]
        )

        ms2 = MappingSet(datasource_name="ChEBI", version="2024-02")
        ms2.add_mappings(
            [
                ChEBIMapping(
                    subject_id="CHEBI:10001",
                    predicate_id="IAO:0100001",
                    object_id="CHEBI:99901",
                ),
                ChEBIMapping(
                    subject_id="CHEBI:10002",
                    predicate_id="IAO:0100001",
                    object_id="CHEBI:99902",
                ),
            ]
        )

        diff = diff_mapping_sets(ms1, ms2)
        assert diff.added_count == 1
        assert diff.removed_count == 0
