"""Integration test: every shipped datasource config must pass schema validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from pysec2pri.config.schema import (
    ConfigValidationError,
    validate_config_dict,
    validate_config_file,
)
from pysec2pri.parsers.base import CONFIG_DIR


@pytest.mark.parametrize("path", sorted(CONFIG_DIR.glob("*.yaml")), ids=lambda p: p.name)
def test_shipped_config_validates(path: Path) -> None:
    """Every config/*.yaml file shipped with the package is schema-valid."""
    validate_config_file(path)


def test_broken_config_fails_validation() -> None:
    """A config dict missing a required field is rejected, not silently accepted."""
    with pytest.raises(ConfigValidationError):
        broken = {"prefix": "X", "curie_base_url": "http://example.org/"}
        validate_config_dict(broken, "broken.yaml")
