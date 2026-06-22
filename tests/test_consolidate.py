"""Tests for pysec2pri.consolidate (per-datasource mapping-date index)."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar

import pytest

from pysec2pri import consolidate as consolidate_module
from pysec2pri.consolidate import (
    SUPPORTED_DATASOURCES,
    consolidate_mapping_dates,
    default_cache_dir,
    load_mapping_dates,
)


class TestDefaultCacheDir:
    """Cache directory resolution."""

    def test_uses_env_var_when_set(self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
        """``$PYSEC2PRI_CACHE_DIR`` takes priority over the default."""
        monkeypatch.setenv("PYSEC2PRI_CACHE_DIR", str(tmp_path))
        assert default_cache_dir() == tmp_path

    def test_falls_back_to_home_cache(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Without the env var, defaults to ``~/.cache/pysec2pri``."""
        monkeypatch.delenv("PYSEC2PRI_CACHE_DIR", raising=False)
        assert default_cache_dir() == Path.home() / ".cache" / "pysec2pri"


class TestLoadMappingDates:
    """Reading the consolidated cache."""

    def test_missing_cache_returns_empty_dict(self, tmp_path: Path) -> None:
        """No behavior change until someone has run ``consolidate``."""
        assert load_mapping_dates("chebi", tmp_path) == {}

    def test_reads_first_seen_date_only(self, tmp_path: Path) -> None:
        """Only ``first_seen_date`` is exposed; last-seen columns are internal."""
        cache_path = consolidate_module._cache_path(tmp_path, "chebi", "3star", "ids")
        consolidate_module._write_cache(
            cache_path,
            {
                "sec2pri:chebi/aaa": {
                    "first_seen_version": "100",
                    "first_seen_date": "2013-02-15",
                    "last_seen_version": "245",
                    "last_seen_date": "2020-01-01",
                }
            },
        )
        result = load_mapping_dates("chebi", tmp_path, subset="3star", mapping_sets="ids")
        assert result == {"sec2pri:chebi/aaa": "2013-02-15"}

    def test_datasource_subset_and_mapping_sets_select_distinct_caches(
        self, tmp_path: Path
    ) -> None:
        """``ids``/``labels`` and datasource caches don't collide."""
        consolidate_module._write_cache(
            consolidate_module._cache_path(tmp_path, "chebi", "3star", "ids"),
            {"r1": {"first_seen_version": "1", "first_seen_date": "2001-01-01"}},
        )
        consolidate_module._write_cache(
            consolidate_module._cache_path(tmp_path, "chebi", "3star", "labels"),
            {"r2": {"first_seen_version": "1", "first_seen_date": "2002-02-02"}},
        )
        consolidate_module._write_cache(
            consolidate_module._cache_path(tmp_path, "uniprot", "3star", "ids"),
            {"r3": {"first_seen_version": "1", "first_seen_date": "2003-03-03"}},
        )
        assert load_mapping_dates("chebi", tmp_path, mapping_sets="ids") == {"r1": "2001-01-01"}
        assert load_mapping_dates("chebi", tmp_path, mapping_sets="labels") == {"r2": "2002-02-02"}
        assert load_mapping_dates("uniprot", tmp_path, mapping_sets="ids") == {"r3": "2003-03-03"}


class _FakeMappingSet:
    """Stand-in for a Sec2PriMappingSet exposing just what consolidate needs."""

    def __init__(
        self,
        record_ids: list[str],
        *,
        dates: dict[str, str] | None = None,
        mapping_set_version: str = "v1",
    ) -> None:
        """Store one SimpleNamespace mapping per record_id, with an optional date."""
        dates = dates or {}
        self.mappings = [
            SimpleNamespace(record_id=rid, mapping_date=dates.get(rid)) for rid in record_ids
        ]
        self.mapping_set_version = mapping_set_version


class _FakeDownloader:
    """Stand-in for a *Downloader class: fixed version list, no real network."""

    versions: ClassVar[list[str]] = []

    def __init__(self, *args: object, **kwargs: object) -> None:
        """Ignore all init args; this fake has no real config to load."""

    def list_versions(self) -> list[str]:
        """Return the fixed test version list."""
        return list(self.versions)


def _fake_release_date(datasource: str, version: str, subset: str = "3star") -> datetime:
    """Fake release date: numeric versions (ChEBI) or ISO dates (HGNC)."""
    try:
        return datetime(2000 + int(version), 1, 1)
    except ValueError:
        return datetime.strptime(version, "%Y-%m-%d")


class TestConsolidateMappingDatesValidation:
    """Upfront validation before any download/parse is attempted."""

    @pytest.mark.parametrize(
        ("kwargs", "match"),
        [
            ({"datasource": "not-a-real-source"}, "Unsupported datasource"),
            ({"datasource": "uniprot", "mapping_sets": "labels"}, "does not support mapping_sets"),
            ({"datasource": "chebi", "mode": "bogus"}, "mode must be"),
            ({"datasource": "ncbi", "mode": "release"}, "no versioned archive"),
        ],
    )
    def test_bad_arguments_raise(self, kwargs: dict[str, Any], match: str) -> None:
        """Unknown datasource/mapping_sets/mode, and release mode on an unarchived source."""
        with pytest.raises(ValueError, match=match):
            consolidate_mapping_dates(**kwargs)

    def test_all_supported_datasources_are_known_to_validation(self) -> None:
        """Every entry in SUPPORTED_DATASOURCES has a mapping_sets capability row."""
        for ds in SUPPORTED_DATASOURCES:
            assert ds in consolidate_module._SUPPORTED_MAPPING_SETS


class TestConsolidateByRelease:
    """Historical-walk 'release' mode, fully monkeypatched (no network)."""

    def test_first_and_last_seen_tracked_across_versions(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A mapping's first_seen is its earliest version; last_seen keeps bumping."""
        records_by_version = {"1": ["a", "b"], "2": ["b", "c"], "3": ["c"]}

        _FakeDownloader.versions = ["1", "2", "3"]
        monkeypatch.setattr("pysec2pri.parsers.chebi.ChEBIDownloader", _FakeDownloader)
        monkeypatch.setattr("pysec2pri.download.resolve_release_date", _fake_release_date)
        monkeypatch.setattr(
            consolidate_module,
            "_run_one_version",
            lambda datasource, version, subset, mapping_sets, tax_id: _FakeMappingSet(
                records_by_version[version]
            ),
        )

        cache_path = consolidate_mapping_dates("chebi", cache_dir=tmp_path, show_progress=False)
        records = consolidate_module._read_cache(cache_path)

        assert records["a"]["first_seen_version"] == "1"
        assert records["a"]["last_seen_version"] == "1"
        assert records["b"]["first_seen_version"] == "1"
        assert records["b"]["last_seen_version"] == "2"
        assert records["c"]["first_seen_version"] == "2"
        assert records["c"]["last_seen_version"] == "3"

        meta_path = consolidate_module._meta_path(tmp_path, "chebi", "3star", "ids")
        assert json.loads(meta_path.read_text())["last_version"] == "3"

    def test_resume_skips_seen_versions_unless_forced(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A second run re-scans nothing; ``force=True`` re-scans everything."""
        seen_versions: list[str] = []

        _FakeDownloader.versions = ["1", "2"]
        monkeypatch.setattr("pysec2pri.parsers.chebi.ChEBIDownloader", _FakeDownloader)
        monkeypatch.setattr("pysec2pri.download.resolve_release_date", _fake_release_date)

        def _run(
            datasource: str, version: str, subset: str, mapping_sets: str, tax_id: str
        ) -> _FakeMappingSet:
            seen_versions.append(version)
            return _FakeMappingSet(["a"])

        monkeypatch.setattr(consolidate_module, "_run_one_version", _run)

        consolidate_mapping_dates("chebi", cache_dir=tmp_path, show_progress=False)
        assert seen_versions == ["1", "2"]

        seen_versions.clear()
        consolidate_mapping_dates("chebi", cache_dir=tmp_path, show_progress=False)
        assert seen_versions == []

        seen_versions.clear()
        consolidate_mapping_dates("chebi", cache_dir=tmp_path, show_progress=False, force=True)
        assert seen_versions == ["1", "2"]

    def test_failed_version_is_skipped_not_fatal(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A per-version exception is logged and skipped; later versions still merge."""
        _FakeDownloader.versions = ["1", "2", "3"]
        monkeypatch.setattr("pysec2pri.parsers.chebi.ChEBIDownloader", _FakeDownloader)
        monkeypatch.setattr("pysec2pri.download.resolve_release_date", _fake_release_date)

        def _run(
            datasource: str, version: str, subset: str, mapping_sets: str, tax_id: str
        ) -> _FakeMappingSet:
            if version == "2":
                raise ValueError("simulated parse failure")
            return _FakeMappingSet(["a"] if version == "1" else ["a", "d"])

        monkeypatch.setattr(consolidate_module, "_run_one_version", _run)

        cache_path = consolidate_mapping_dates("chebi", cache_dir=tmp_path, show_progress=False)
        records = consolidate_module._read_cache(cache_path)

        assert records["a"]["first_seen_version"] == "1"
        assert records["a"]["last_seen_version"] == "3"
        assert records["d"]["first_seen_version"] == "3"

        meta_path = consolidate_module._meta_path(tmp_path, "chebi", "3star", "ids")
        assert json.loads(meta_path.read_text())["last_version"] == "3"

    def test_works_for_any_versioned_datasource_not_just_chebi(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """HGNC (a second versioned datasource) walks the same generic code path."""
        _FakeDownloader.versions = ["2024-01-01", "2024-04-01"]
        monkeypatch.setattr("pysec2pri.parsers.hgnc.HGNCDownloader", _FakeDownloader)
        monkeypatch.setattr("pysec2pri.download.resolve_release_date", _fake_release_date)
        monkeypatch.setattr(
            consolidate_module,
            "_run_one_version",
            lambda datasource, version, subset, mapping_sets, tax_id: _FakeMappingSet(["x"]),
        )

        cache_path = consolidate_mapping_dates("hgnc", cache_dir=tmp_path, show_progress=False)
        records = consolidate_module._read_cache(cache_path)
        assert records["x"]["first_seen_version"] == "2024-01-01"
        assert records["x"]["last_seen_version"] == "2024-04-01"


class TestConsolidateByDate:
    """Single-pass 'date' mode, fully monkeypatched (no network)."""

    def test_uses_real_per_row_date_directly(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """When real per-row dates exist, they're captured without any historical walk."""
        monkeypatch.setattr(
            consolidate_module,
            "_run_one_version",
            lambda datasource, version, subset, mapping_sets, tax_id: _FakeMappingSet(
                ["a", "b"], dates={"a": "1996-03-01", "b": "2010-05-12"}
            ),
        )

        cache_path = consolidate_mapping_dates(
            "hgnc", mode="date", cache_dir=tmp_path, show_progress=False
        )
        records = consolidate_module._read_cache(cache_path)
        assert records["a"]["first_seen_date"] == "1996-03-01"
        assert records["b"]["first_seen_date"] == "2010-05-12"

    def test_falls_back_to_release_with_warning_when_no_dates(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """ChEBI never produces real per-row dates: 'date' mode must warn and use 'release'."""
        _FakeDownloader.versions = ["1"]
        monkeypatch.setattr("pysec2pri.parsers.chebi.ChEBIDownloader", _FakeDownloader)
        monkeypatch.setattr("pysec2pri.download.resolve_release_date", _fake_release_date)
        monkeypatch.setattr(
            consolidate_module,
            "_run_one_version",
            lambda datasource, version, subset, mapping_sets, tax_id: _FakeMappingSet(["a"]),
        )

        with pytest.warns(UserWarning, match="falling back to mode='release'"):
            cache_path = consolidate_mapping_dates(
                "chebi", mode="date", cache_dir=tmp_path, show_progress=False
            )

        records = consolidate_module._read_cache(cache_path)
        assert records["a"]["first_seen_version"] == "1"
