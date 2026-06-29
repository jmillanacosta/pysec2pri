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
    _label_transitions,
    build_label_history,
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
        cache_path = consolidate_module._cache_path(tmp_path, "chebi", "ids", subset="3star")
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

    def test_excludes_records_with_no_resolved_date(self, tmp_path: Path) -> None:
        """A record walked but never assigned a real date is omitted, not '' or 'None'."""
        consolidate_module._write_cache(
            consolidate_module._cache_path(tmp_path, "chebi", "ids"),
            {
                "dated": {"first_seen_version": "245", "first_seen_date": "2020-01-01"},
                "undated": {"first_seen_version": "183", "first_seen_date": ""},
            },
        )
        assert load_mapping_dates("chebi", tmp_path, mapping_sets="ids") == {"dated": "2020-01-01"}

    def test_datasource_subset_and_mapping_sets_select_distinct_caches(
        self, tmp_path: Path
    ) -> None:
        """``ids``/``labels`` and datasource caches don't collide."""
        consolidate_module._write_cache(
            consolidate_module._cache_path(tmp_path, "chebi", "ids"),
            {"r1": {"first_seen_version": "1", "first_seen_date": "2001-01-01"}},
        )
        consolidate_module._write_cache(
            consolidate_module._cache_path(tmp_path, "chebi", "labels"),
            {"r2": {"first_seen_version": "1", "first_seen_date": "2002-02-02"}},
        )
        consolidate_module._write_cache(
            consolidate_module._cache_path(tmp_path, "uniprot", "ids"),
            {"r3": {"first_seen_version": "1", "first_seen_date": "2003-03-03"}},
        )
        assert load_mapping_dates("chebi", tmp_path, mapping_sets="ids") == {"r1": "2001-01-01"}
        assert load_mapping_dates("chebi", tmp_path, mapping_sets="labels") == {"r2": "2002-02-02"}
        assert load_mapping_dates("uniprot", tmp_path, mapping_sets="ids") == {"r3": "2003-03-03"}


class _FakeMappingSet:
    """Stand-in for a BaseMappingSet exposing just what consolidate needs."""

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
        monkeypatch.setattr("pysec2pri.downloads.ChEBIDownloader", _FakeDownloader)
        monkeypatch.setattr("pysec2pri.download.resolve_release_date", _fake_release_date)
        monkeypatch.setattr(
            consolidate_module,
            "_run_one_version",
            lambda datasource, version, mapping_sets, **kwargs: _FakeMappingSet(
                records_by_version[version]
            ),
        )

        cache_path, _ = consolidate_mapping_dates("chebi", cache_dir=tmp_path, show_progress=False)
        records = consolidate_module._read_cache(cache_path)

        assert records["a"]["first_seen_version"] == "1"
        assert records["a"]["last_seen_version"] == "1"
        assert records["b"]["first_seen_version"] == "1"
        assert records["b"]["last_seen_version"] == "2"
        assert records["c"]["first_seen_version"] == "2"
        assert records["c"]["last_seen_version"] == "3"

        meta_path = consolidate_module._meta_path(tmp_path, "chebi", "ids")
        assert json.loads(meta_path.read_text())["last_version"] == "3"

    def test_resume_skips_seen_versions_unless_forced(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A second run re-scans nothing; ``force=True`` re-scans everything."""
        seen_versions: list[str] = []

        _FakeDownloader.versions = ["1", "2"]
        monkeypatch.setattr("pysec2pri.downloads.ChEBIDownloader", _FakeDownloader)
        monkeypatch.setattr("pysec2pri.download.resolve_release_date", _fake_release_date)

        def _run(
            datasource: str, version: str, mapping_sets: str, **kwargs: object
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
        monkeypatch.setattr("pysec2pri.downloads.ChEBIDownloader", _FakeDownloader)
        monkeypatch.setattr("pysec2pri.download.resolve_release_date", _fake_release_date)

        def _run(
            datasource: str, version: str, mapping_sets: str, **kwargs: object
        ) -> _FakeMappingSet:
            if version == "2":
                raise ValueError("simulated parse failure")
            return _FakeMappingSet(["a"] if version == "1" else ["a", "d"])

        monkeypatch.setattr(consolidate_module, "_run_one_version", _run)

        cache_path, _ = consolidate_mapping_dates("chebi", cache_dir=tmp_path, show_progress=False)
        records = consolidate_module._read_cache(cache_path)

        assert records["a"]["first_seen_version"] == "1"
        assert records["a"]["last_seen_version"] == "3"
        assert records["d"]["first_seen_version"] == "3"

        meta_path = consolidate_module._meta_path(tmp_path, "chebi", "ids")
        assert json.loads(meta_path.read_text())["last_version"] == "3"

    def test_unresolvable_release_date_stores_empty_not_the_version_label(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """When no real release date is found, the version number is not stored as a date."""
        _FakeDownloader.versions = ["183"]
        monkeypatch.setattr("pysec2pri.downloads.ChEBIDownloader", _FakeDownloader)
        monkeypatch.setattr("pysec2pri.download.resolve_release_date", lambda *a, **kw: None)
        monkeypatch.setattr(
            consolidate_module,
            "_run_one_version",
            lambda datasource, version, mapping_sets, **kwargs: _FakeMappingSet(["a"]),
        )

        cache_path, _ = consolidate_mapping_dates("chebi", cache_dir=tmp_path, show_progress=False)
        records = consolidate_module._read_cache(cache_path)

        assert records["a"]["first_seen_version"] == "183"
        assert records["a"]["first_seen_date"] == ""
        assert records["a"]["last_seen_date"] == ""

    def test_works_for_any_versioned_datasource_not_just_chebi(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """HGNC (a second versioned datasource) walks the same generic code path."""
        _FakeDownloader.versions = ["2024-01-01", "2024-04-01"]
        monkeypatch.setattr("pysec2pri.downloads.HGNCDownloader", _FakeDownloader)
        monkeypatch.setattr("pysec2pri.download.resolve_release_date", _fake_release_date)
        monkeypatch.setattr(
            consolidate_module,
            "_run_one_version",
            lambda datasource, version, mapping_sets, **kwargs: _FakeMappingSet(["x"]),
        )

        cache_path, _ = consolidate_mapping_dates("hgnc", cache_dir=tmp_path, show_progress=False)
        records = consolidate_module._read_cache(cache_path)
        assert records["x"]["first_seen_version"] == "2024-01-01"
        assert records["x"]["last_seen_version"] == "2024-04-01"

    def test_output_writes_consolidated_sssom_to_custom_path(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """``output=`` also writes the consolidated mapping set to the given path."""
        _FakeDownloader.versions = ["1", "2"]
        monkeypatch.setattr("pysec2pri.downloads.ChEBIDownloader", _FakeDownloader)
        monkeypatch.setattr("pysec2pri.download.resolve_release_date", _fake_release_date)
        monkeypatch.setattr(
            consolidate_module,
            "_run_one_version",
            lambda datasource, version, mapping_sets, **kwargs: _FakeMappingSet(["a"]),
        )

        out = tmp_path / "subdir" / "custom_consolidated.sssom.tsv"
        consolidate_mapping_dates(
            "chebi", cache_dir=tmp_path, show_progress=False, output=out
        )
        assert out.exists()


class TestBuildConsolidatedMappingSet:
    """Materializing the cached field snapshots back into a real SSSOM mapping set."""

    def test_empty_first_seen_date_leaves_mapping_date_unset(self) -> None:
        """An empty (unresolvable) first_seen_date must not reach Mapping(mapping_date=...)."""
        fields_json = json.dumps(
            {
                "subject_id": "CHEBI:10001",
                "object_id": "CHEBI:99901",
                "predicate_id": "IAO:0100001",
                "mapping_justification": "semapv:BackgroundKnowledgeBasedMatching",
            }
        )
        records = {
            "a" * 16: {
                "first_seen_version": "183",
                "first_seen_date": "",
                "last_seen_version": "183",
                "last_seen_date": "",
                "fields_json": fields_json,
            }
        }

        mapping_set = consolidate_module._build_consolidated_mapping_set(
            "chebi", "ids", records, "183"
        )

        assert mapping_set.mappings[0].mapping_date is None


class TestConsolidateByDate:
    """Single-pass 'date' mode, fully monkeypatched (no network)."""

    def test_uses_real_per_row_date_directly(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """When real per-row dates exist, they're captured without any historical walk."""
        monkeypatch.setattr(
            consolidate_module,
            "_run_one_version",
            lambda datasource, version, mapping_sets, **kwargs: _FakeMappingSet(
                ["a", "b"], dates={"a": "1996-03-01", "b": "2010-05-12"}
            ),
        )

        cache_path, _ = consolidate_mapping_dates(
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
        monkeypatch.setattr("pysec2pri.downloads.ChEBIDownloader", _FakeDownloader)
        monkeypatch.setattr("pysec2pri.download.resolve_release_date", _fake_release_date)
        monkeypatch.setattr(
            consolidate_module,
            "_run_one_version",
            lambda datasource, version, mapping_sets, **kwargs: _FakeMappingSet(["a"]),
        )

        with pytest.warns(UserWarning, match="falling back to mode='release'"):
            cache_path, _ = consolidate_mapping_dates(
                "chebi", mode="date", cache_dir=tmp_path, show_progress=False
            )

        records = consolidate_module._read_cache(cache_path)
        assert records["a"]["first_seen_version"] == "1"


class TestLabelTransitions:
    """Pure-function tests for the Ensembl label-history diff logic."""

    def test_changed_label_is_reported(self) -> None:
        """A gene present in both snapshots with a different label is reported."""
        prev = {"G1": "OLD1"}
        curr = {"G1": "NEW1"}
        assert _label_transitions(prev, curr) == [("G1", "OLD1", "NEW1")]

    def test_unchanged_label_is_not_reported(self) -> None:
        """A gene whose label didn't change produces no transition."""
        prev = {"G1": "SAME"}
        curr = {"G1": "SAME"}
        assert _label_transitions(prev, curr) == []

    def test_gene_only_in_one_snapshot_is_not_reported(self) -> None:
        """A brand-new gene, or one absent from the later snapshot, has nothing to diff."""
        prev = {"G1": "OLD1", "G2": "GONE"}
        curr = {"G1": "OLD1", "G3": "NEW"}
        assert _label_transitions(prev, curr) == []

    def test_multiple_changes_all_reported(self) -> None:
        """Every changed gene gets its own transition tuple."""
        prev = {"G1": "OLD1", "G2": "OLD2", "G3": "SAME"}
        curr = {"G1": "NEW1", "G2": "NEW2", "G3": "SAME"}
        assert set(_label_transitions(prev, curr)) == {
            ("G1", "OLD1", "NEW1"),
            ("G2", "OLD2", "NEW2"),
        }


class TestBuildLabelHistory:
    """Cross-release Ensembl label-history walk, fully monkeypatched (no network)."""

    def test_rejects_all_species(self, tmp_path: Path) -> None:
        """species='all' (the config default) is rejected with a clear message."""
        with pytest.raises(ValueError, match="requires an explicit single species"):
            build_label_history(species="all", cache_dir=tmp_path, show_progress=False)

    def _patch_walk(
        self,
        monkeypatch: pytest.MonkeyPatch,
        versions: list[str],
        label_maps: dict[str, dict[str, str]],
        seen_versions: list[str] | None = None,
    ) -> None:
        from pysec2pri.parsers.ensembl import EnsemblParser

        monkeypatch.setitem(
            consolidate_module._LIST_VERSIONS_FNS, "ensembl", lambda **kwargs: versions
        )
        monkeypatch.setattr("pysec2pri.download.resolve_release_date", lambda *a, **kw: None)

        def _fake_download(
            datasource: str,
            output_dir: Path,
            version: str | None = None,
            species: object = None,
            keys: list[str] | None = None,
            **kwargs: object,
        ) -> tuple[dict[str, Path], None]:
            return {"gene": Path(f"gene_{version}.txt"), "xref": Path(f"xref_{version}.txt")}, None

        monkeypatch.setattr("pysec2pri.download.download_datasource_with_release", _fake_download)

        def _fake_snapshot(self: EnsemblParser, gene_path: Path, xref_path: Path) -> dict[str, str]:
            version = Path(gene_path).stem.split("_")[-1]
            if seen_versions is not None:
                seen_versions.append(version)
            return label_maps[version]

        monkeypatch.setattr(EnsemblParser, "current_label_snapshot", _fake_snapshot)

    def test_detects_transitions_across_releases(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Diffing consecutive release snapshots emits one transition per changed gene."""
        label_maps = {
            "113": {"G1": "OLD1", "G2": "SAME"},
            "114": {"G1": "NEW1", "G2": "SAME"},
            "115": {"G1": "NEW1", "G2": "SAME", "G3": "BRANDNEW"},
        }
        self._patch_walk(monkeypatch, ["113", "114", "115"], label_maps)

        ms = build_label_history(species=9606, cache_dir=tmp_path, show_progress=False)
        seen = {(m.object_id, m.subject_label, m.object_label) for m in ms.mappings}
        assert seen == {("ENSEMBL:G1", "OLD1", "NEW1")}
        assert all(m.predicate_id == "IAO:0100001" for m in ms.mappings)

    def test_from_to_version_bounds_the_walk(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """from_version/to_version restrict which releases are walked."""
        label_maps = {
            "113": {"G1": "A"},
            "114": {"G1": "B"},
            "115": {"G1": "C"},
        }
        seen: list[str] = []
        self._patch_walk(monkeypatch, ["113", "114", "115"], label_maps, seen_versions=seen)

        build_label_history(
            species=9606,
            from_version="114",
            to_version="115",
            cache_dir=tmp_path,
            show_progress=False,
        )
        assert seen == ["114", "115"]

    def test_resume_skips_already_walked_releases(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A second run only walks releases past the resumed last_version."""
        label_maps = {
            "113": {"G1": "A"},
            "114": {"G1": "B"},
            "115": {"G1": "C"},
        }
        seen: list[str] = []
        self._patch_walk(monkeypatch, ["113", "114", "115"], label_maps, seen_versions=seen)

        build_label_history(species=9606, cache_dir=tmp_path, show_progress=False)
        assert seen == ["113", "114", "115"]

        seen.clear()
        ms = build_label_history(species=9606, cache_dir=tmp_path, show_progress=False)
        assert seen == []
        # Resuming still returns the full set of transitions found so far.
        assert {(m.subject_label, m.object_label) for m in ms.mappings} == {("A", "B"), ("B", "C")}

    def test_force_rewalks_every_release(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """force=True ignores resume state and re-walks from scratch."""
        label_maps = {"113": {"G1": "A"}, "114": {"G1": "B"}}
        seen: list[str] = []
        self._patch_walk(monkeypatch, ["113", "114"], label_maps, seen_versions=seen)

        build_label_history(species=9606, cache_dir=tmp_path, show_progress=False)
        seen.clear()
        build_label_history(species=9606, cache_dir=tmp_path, show_progress=False, force=True)
        assert seen == ["113", "114"]
