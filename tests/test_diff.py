"""Tests for the diff module.

Tests cover:
- MappingDiff dataclass properties
- mapping_set_to_dataframe conversion
- diff_mapping_sets comparison logic
- summarize_diff output formatting
"""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING

import polars as pl
import pytest

from pysec2pri.diff import (
    MappingDiff,
    diff_mapping_sets,
    diff_sssom_files,
    mapping_set_to_dataframe,
    summarize_diff,
)

if TYPE_CHECKING:
    pass

from sssom_schema import Mapping, MappingSet

# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def empty_mapping_set() -> MappingSet:
    """Create an empty MappingSet."""
    return MappingSet(mapping_set_id="test:empty", license="https://example.org/license")


@pytest.fixture
def simple_mapping_set_v1() -> MappingSet:
    """Create a simple MappingSet (version 1)."""
    ms = MappingSet(
        mapping_set_id="test:v1",
        license="https://example.org/license",
        mapping_set_version="1.0",
    )
    ms.mappings = [
        Mapping(
            subject_id="PRIMARY1",
            predicate_id="IAO:0100001",
            object_id="SECONDARY1",
            mapping_justification="semapv:BackgroundKnowledgeBasedMatching",
        ),
        Mapping(
            subject_id="PRIMARY2",
            predicate_id="IAO:0100001",
            object_id="SECONDARY2",
            mapping_justification="semapv:BackgroundKnowledgeBasedMatching",
        ),
        Mapping(
            subject_id="PRIMARY3",
            predicate_id="IAO:0100001",
            object_id="SECONDARY3",
            mapping_justification="semapv:BackgroundKnowledgeBasedMatching",
        ),
    ]
    return ms


@pytest.fixture
def simple_mapping_set_v2() -> MappingSet:
    """Create a simple MappingSet (version 2) with changes."""
    ms = MappingSet(
        mapping_set_id="test:v2",
        license="https://example.org/license",
        mapping_set_version="2.0",
    )
    ms.mappings = [
        # Unchanged
        Mapping(
            subject_id="PRIMARY1",
            predicate_id="IAO:0100001",
            object_id="SECONDARY1",
            mapping_justification="semapv:BackgroundKnowledgeBasedMatching",
        ),
        # Changed: SECONDARY2 now maps to PRIMARY4 instead of PRIMARY2
        Mapping(
            subject_id="PRIMARY4",
            predicate_id="IAO:0100001",
            object_id="SECONDARY2",
            mapping_justification="semapv:BackgroundKnowledgeBasedMatching",
        ),
        # New mapping (added)
        Mapping(
            subject_id="PRIMARY5",
            predicate_id="IAO:0100001",
            object_id="SECONDARY5",
            mapping_justification="semapv:BackgroundKnowledgeBasedMatching",
        ),
        # SECONDARY3 -> PRIMARY3 removed (not in v2)
    ]
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
        assert isinstance(diff.added, pl.DataFrame)
        assert isinstance(diff.removed, pl.DataFrame)
        assert isinstance(diff.changed, pl.DataFrame)

    def test_added_count(self) -> None:
        """Test added_count property."""
        diff = MappingDiff(
            old_version="1.0",
            new_version="2.0",
            datasource="Test",
            added=pl.DataFrame({"subject_id": ["A", "B"], "object_id": ["X", "Y"]}),
        )
        assert diff.added_count == 2

    def test_removed_count(self) -> None:
        """Test removed_count property."""
        diff = MappingDiff(
            old_version="1.0",
            new_version="2.0",
            datasource="Test",
            removed=pl.DataFrame({"subject_id": ["A"], "object_id": ["X"]}),
        )
        assert diff.removed_count == 1

    def test_changed_count(self) -> None:
        """Test changed_count property."""
        diff = MappingDiff(
            old_version="1.0",
            new_version="2.0",
            datasource="Test",
            changed=pl.DataFrame(
                {
                    "object_id": ["X"],
                    "old_subject_id": ["A"],
                    "new_subject_id": ["B"],
                }
            ),
        )
        assert diff.changed_count == 1

    def test_total_changes(self) -> None:
        """Test total_changes property."""
        diff = MappingDiff(
            old_version="1.0",
            new_version="2.0",
            datasource="Test",
            added=pl.DataFrame({"subject_id": ["A"], "object_id": ["X"]}),
            removed=pl.DataFrame({"subject_id": ["B"], "object_id": ["Y"]}),
            changed=pl.DataFrame(
                {
                    "object_id": ["Z"],
                    "old_subject_id": ["C"],
                    "new_subject_id": ["D"],
                }
            ),
        )
        assert diff.added_count == 1
        assert diff.removed_count == 1
        assert diff.changed_count == 1
        assert diff.total_changes == 3

    def test_has_changes_true(self) -> None:
        """Test has_changes is True when there are changes."""
        diff = MappingDiff(
            old_version="1.0",
            new_version="2.0",
            datasource="Test",
            added=pl.DataFrame({"subject_id": ["A"], "object_id": ["X"]}),
        )
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

    def test_diff_with_datasource(self, simple_mapping_set_v1: MappingSet) -> None:
        """Test diff includes datasource name."""
        diff = diff_mapping_sets(
            simple_mapping_set_v1,
            simple_mapping_set_v1,
            datasource="TestSource",
        )
        assert diff.datasource == "TestSource"

    def test_diff_versions(
        self, simple_mapping_set_v1: MappingSet, simple_mapping_set_v2: MappingSet
    ) -> None:
        """Test diff captures versions."""
        diff = diff_mapping_sets(simple_mapping_set_v1, simple_mapping_set_v2)
        assert diff.old_version == "1.0"
        assert diff.new_version == "2.0"


# =============================================================================
# Tests for diff_sssom_files
# =============================================================================


class TestDiffSssomFiles:
    """Tests for diff_sssom_files function."""

    def test_diff_sssom_files(self) -> None:
        """Test diffing SSSOM TSV files."""
        with TemporaryDirectory() as tmpdir:
            old_file = Path(tmpdir) / "old.tsv"
            new_file = Path(tmpdir) / "new.tsv"

            # Create old file
            old_content = """#mapping_set_version: "1.0"
subject_id\tobject_id
PRIMARY1\tSECONDARY1
PRIMARY2\tSECONDARY2
"""
            old_file.write_text(old_content)

            # Create new file with changes
            new_content = """#mapping_set_version: "2.0"
subject_id\tobject_id
PRIMARY1\tSECONDARY1
PRIMARY3\tSECONDARY3
"""
            new_file.write_text(new_content)

            diff = diff_sssom_files(old_file, new_file, datasource="test")

            assert diff.old_version == "1.0"
            assert diff.new_version == "2.0"
            # SECONDARY2 removed, SECONDARY3 added
            assert diff.removed_count == 1
            assert diff.added_count == 1

    def test_diff_sssom_files_no_version(self) -> None:
        """Test diffing SSSOM files without version metadata."""
        with TemporaryDirectory() as tmpdir:
            old_file = Path(tmpdir) / "old.tsv"
            new_file = Path(tmpdir) / "new.tsv"

            old_content = """subject_id\tobject_id
PRIMARY1\tSECONDARY1
"""
            old_file.write_text(old_content)

            new_content = """subject_id\tobject_id
PRIMARY1\tSECONDARY1
"""
            new_file.write_text(new_content)

            diff = diff_sssom_files(old_file, new_file)

            assert diff.old_version is None
            assert diff.new_version is None
            assert diff.has_changes is False


# =============================================================================
# Tests for summarize_diff
# =============================================================================


class TestSummarizeDiff:
    """Tests for summarize_diff function."""

    def test_summarize_empty_diff(self) -> None:
        """Test summarizing an empty diff."""
        diff = MappingDiff(
            old_version="1.0",
            new_version="2.0",
            datasource="TestDB",
        )
        summary = summarize_diff(diff)

        assert "TestDB" in summary
        assert "1.0" in summary
        assert "2.0" in summary
        assert "Added mappings:" in summary
        assert "0" in summary

    def test_summarize_with_changes(self) -> None:
        """Test summarizing a diff with changes."""
        diff = MappingDiff(
            old_version="1.0",
            new_version="2.0",
            datasource="TestDB",
            added=pl.DataFrame({"subject_id": ["A"], "object_id": ["X"]}),
            removed=pl.DataFrame({"subject_id": ["B"], "object_id": ["Y"]}),
        )
        summary = summarize_diff(diff)

        assert "Added mappings:" in summary
        assert "Removed mappings:" in summary
        assert "X" in summary  # Added object shown
        assert "Y" in summary  # Removed object shown

    def test_summarize_with_changed_mappings(self) -> None:
        """Test summarizing a diff with changed mappings."""
        diff = MappingDiff(
            old_version="1.0",
            new_version="2.0",
            datasource="TestDB",
            changed=pl.DataFrame(
                {
                    "object_id": ["SEC1"],
                    "old_subject_id": ["OLD_PRI"],
                    "new_subject_id": ["NEW_PRI"],
                }
            ),
        )
        summary = summarize_diff(diff)

        assert "Changed:" in summary
        assert "SEC1" in summary
        assert "OLD_PRI" in summary
        assert "NEW_PRI" in summary

    def test_summarize_truncates_large_lists(self) -> None:
        """Test that large change lists are not fully printed."""
        # Create more than 10 added mappings
        added = pl.DataFrame(
            {
                "subject_id": [f"PRI{i}" for i in range(20)],
                "object_id": [f"SEC{i}" for i in range(20)],
            }
        )
        diff = MappingDiff(
            old_version="1.0",
            new_version="2.0",
            datasource="TestDB",
            added=added,
        )
        summary = summarize_diff(diff)

        # Should show count but not list all 20
        assert "20" in summary
        # SEC19 should not be in the summary since it truncates after 10
        assert "SEC19" not in summary
