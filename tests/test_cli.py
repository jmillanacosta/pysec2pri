"""Tests for the dynamic generate-command factory in pysec2pri.cli."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import click
from click.testing import CliRunner

from pysec2pri.cli import _make_generate_cmd


class _FakeMappingSet:
    """Minimal stand-in for a Sec2PriMappingSet, just enough for save()/_emit()."""

    def __init__(self, version: str | None) -> None:
        """Store the version and an empty mappings list."""
        self.mapping_set_version = version
        self.mappings: list[object] = []

    def save(self, fmt: str, output_path: Path, **kwargs: object) -> Path:
        """Write a placeholder file and return its path, mirroring the real save method."""
        Path(output_path).write_text("fake")
        return Path(output_path)


def _fake_generate(**_: Any) -> _FakeMappingSet:
    return _FakeMappingSet(version="115")


class TestSpeciesInOutputFilename:
    """A species-aware command's default output filename must include the species code.

    Otherwise two runs for different species (e.g. human vs. dog) would
    silently overwrite the same default file.
    """

    def test_species_is_folded_into_default_filename(self) -> None:
        """--species 9615 produces a filename containing "9615", not just the version."""
        cmd_fn = _make_generate_cmd(
            "ensembl", "ids", _fake_generate, [click.option("--species", default="9606")]
        )
        cmd = click.command(name="ids")(cmd_fn)

        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(cmd, ["--species", "9615"])
            assert result.exit_code == 0, result.output
            assert Path("ensembl_ids_9615_115_sssom.tsv").exists()

    def test_no_species_option_omits_species_segment(self) -> None:
        """A command with no --species option (e.g. hgnc) keeps the old filename shape."""
        cmd_fn = _make_generate_cmd("hgnc", "ids", _fake_generate, [])
        cmd = click.command(name="ids")(cmd_fn)

        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(cmd, [])
            assert result.exit_code == 0, result.output
            assert Path("hgnc_ids_115_sssom.tsv").exists()
