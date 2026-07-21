"""Test the copy-paste example in the docs."""

from __future__ import annotations

import importlib.util
import inspect
import shutil
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml
from mapkgsutils.config.schema import validate_config_file

EXAMPLE_DIR = Path(__file__).parent.parent / "docs" / "source" / "example"
EXAMPLE_CONFIG = EXAMPLE_DIR / "example.yaml"


def _load(name: str, path: Path) -> Any:
    """Import a module from an explicit path (the docs are not a package)."""
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def example_parser_cls(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Any:
    """Return ExampleParser bound to the example config.

    Copying a source's config into ``pysec2pri/config/`` is step 1 of adding
    one; here it goes into an importable throwaway package instead, so the
    example never registers itself as a real datasource.
    """
    pkg = tmp_path / "example_config"
    pkg.mkdir()
    (pkg / "__init__.py").touch()
    shutil.copy(EXAMPLE_CONFIG, pkg / "example.yaml")
    monkeypatch.syspath_prepend(str(tmp_path))

    module = _load("docs_example_parser", EXAMPLE_DIR / "example.py")
    return type("BoundExampleParser", (module.ExampleParser,), {"config_package": "example_config"})


@pytest.fixture
def example_files(tmp_path: Path) -> tuple[Path, Path]:
    """Write the two files the example's config and parser describe."""
    withdrawn = tmp_path / "withdrawn.tsv"
    withdrawn.write_text(
        "example_id\treplacement_id\tstatus\n"
        "EX:101\tEX:900\tmerged\n"
        "EX:102\tEX:900\tmerged\n"
        "EX:103\t\twithdrawn\n"
    )
    complete = tmp_path / "complete.tsv"
    complete.write_text(
        "example_id\tsymbol\tprev_symbols\talias_symbols\n"
        "EX:900\tABC1\tABC-old\tABCX|ABCY\n"
        "EX:901\tDEF2\t\t\n"
    )
    return withdrawn, complete


def test_example_config_is_valid() -> None:
    """The example config passes the same schema check as the real ones."""
    validate_config_file(EXAMPLE_CONFIG)


def test_example_config_matches_its_parser(example_parser_cls: Any) -> None:
    """Every method and input the example config names exists on its parser."""
    raw = yaml.safe_load(EXAMPLE_CONFIG.read_text(encoding="utf-8"))
    assert raw["parser_class"] == example_parser_cls.__mro__[1].__name__
    for kind, spec in raw["mapping_sets"].items():
        method = getattr(example_parser_cls, spec["method"], None)
        assert method is not None, f"{kind}.method={spec['method']!r} is not on the parser"
        params = set(inspect.signature(method).parameters)
        for key, param in (spec.get("inputs") or {}).items():
            assert param in params, f"{kind}.inputs[{key!r}] -> {param!r} is not a parameter"
        for key in spec.get("inputs") or {}:
            assert key in raw["download_urls"], f"{kind}.inputs[{key!r}] is not a download key"


def test_example_parses_ids(example_parser_cls: Any, example_files: tuple[Path, Path]) -> None:
    """Retired IDs map to their replacement, or to NoTermFound."""
    withdrawn, complete = example_files
    ms = example_parser_cls(show_progress=False).parse(withdrawn, complete)
    got = {(m.subject_id, m.predicate_id, m.object_id) for m in ms.mappings}
    assert got == {
        ("EX:101", "IAO:0100001", "EX:900"),
        ("EX:102", "IAO:0100001", "EX:900"),
        ("EX:103", "oboInOwl:consider", "sssom:NoTermFound"),
    }
    # The complete set gives the full current list.
    assert object.__getattribute__(ms, "_primary_ids") == {"EX:900", "EX:901"}


def test_example_parses_labels(example_parser_cls: Any, example_files: tuple[Path, Path]) -> None:
    """Previous symbols and aliases land in one set, under their own predicates."""
    _, complete = example_files
    ms = example_parser_cls(show_progress=False).parse_labels(complete)
    got = {(m.subject_label, m.predicate_id, m.object_label) for m in ms.mappings}
    assert got == {
        ("ABC-old", "IAO:0100001", "ABC1"),
        ("ABCX", "oboInOwl:hasExactSynonym", "ABC1"),
        ("ABCY", "oboInOwl:hasExactSynonym", "ABC1"),
    }


def test_example_parses_primary_ids(
    example_parser_cls: Any, example_files: tuple[Path, Path]
) -> None:
    """The primary_ids set carries the current IDs and no mappings."""
    _, complete = example_files
    ms = example_parser_cls(show_progress=False).parse_primary_ids(complete)
    assert not ms.mappings
    assert object.__getattribute__(ms, "_primary_ids") == {"EX:900", "EX:901"}


def test_example_downloader_builds_urls() -> None:
    """The example downloader's URLs match the config's download_urls keys."""
    module = _load("docs_example_downloads", EXAMPLE_DIR / "example_downloads.py")
    urls = module._urls_for_version("2026-01-01")
    raw = yaml.safe_load(EXAMPLE_CONFIG.read_text(encoding="utf-8"))
    assert set(urls) == set(raw["download_urls"])
    assert all("2026-01-01" in u for u in urls.values())
