"""Integration test: every shipped datasource config must pass schema validation."""

from __future__ import annotations

from pathlib import Path

import pytest
from mapkgsutils.config.schema import validate_config_file

from pysec2pri.parsers.base import CONFIG_DIR


@pytest.mark.parametrize("path", sorted(CONFIG_DIR.glob("*.yaml")), ids=lambda p: p.name)
def test_shipped_config_validates(path: Path) -> None:
    """Every config/*.yaml file shipped with the package is schema-valid."""
    validate_config_file(path)
