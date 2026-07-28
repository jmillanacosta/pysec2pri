"""Tests for the dynamic generate-command factory in pysec2pri.cli."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import click
import pytest
from click.testing import CliRunner

from pysec2pri.cli import _make_generate_cmd


class _FakeMappingSet:
    """Minimal stand-in for a BaseMappingSet, just enough for save()/_emit()."""

    def __init__(self, version: str | None) -> None:
        """Store the version and an empty mappings list."""
        self.mapping_set_version = version
        self.mappings: list[object] = []

    def save(self, fmt: str, output_path: Path, **kwargs: object) -> Path:
        """Write a placeholder file and return its path, mirroring the real save method."""
        Path(output_path).write_text("fake")
        return Path(output_path)


def _fake_generate(*_: Any, **__: Any) -> _FakeMappingSet:
    return _FakeMappingSet(version="115")


class TestSpeciesInOutputFilename:
    """A species-aware command's default output filename must include the species code.

    Otherwise two runs for different species (e.g. human vs. dog) would
    overwrite the same default file.
    """

    def test_species_in_default_filename(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """--species 9615 produces a filename containing "9615"."""
        import pysec2pri.api as api

        monkeypatch.setattr(api, "_generate", _fake_generate)
        cmd_fn = _make_generate_cmd("ensembl", "ids", [click.option("--species", default="9606")])
        cmd = click.command(name="ids")(cmd_fn)

        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(cmd, ["--species", "9615"])
            assert result.exit_code == 0, result.output
            assert (api.DEFAULT_OUTPUT_DIR / "ensembl_ids_9615_115.sssom.tsv").exists()

    def test_no_species_option_omits_species_segment(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """A command with no --species option (e.g. hgnc) keeps the old filename shape."""
        import pysec2pri.api as api

        monkeypatch.setattr(api, "_generate", _fake_generate)
        cmd_fn = _make_generate_cmd("hgnc", "ids", [])
        cmd = click.command(name="ids")(cmd_fn)

        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(cmd, [])
            assert result.exit_code == 0, result.output
            assert (api.DEFAULT_OUTPUT_DIR / "hgnc_ids_115.sssom.tsv").exists()


class TestConsolidateInOutputFilename:
    """A consolidated run is a distinct product and must not overwrite the plain run."""

    def test_consolidate_is_folded_into_default_filename(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """--consolidate produces its own file, leaving the plain default untouched."""
        import pysec2pri.api as api

        monkeypatch.setattr(api, "_generate", _fake_generate)
        cmd_fn = _make_generate_cmd("hgnc", "ids", [])
        cmd = click.command(name="ids")(cmd_fn)

        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(cmd, ["--consolidate"])
            assert result.exit_code == 0, result.output
            assert (api.DEFAULT_OUTPUT_DIR / "hgnc_ids_115_consolidate.sssom.tsv").exists()
            assert not (api.DEFAULT_OUTPUT_DIR / "hgnc_ids_115.sssom.tsv").exists()
