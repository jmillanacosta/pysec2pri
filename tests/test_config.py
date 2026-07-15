"""Integration tests: each datasource config must match the code it drives."""

from __future__ import annotations

import inspect
import re
from pathlib import Path
from typing import Any

import pytest
from mapkgsutils.config.schema import validate_config_file
from mapkgsutils.parsers.config import get_datasource_config, product_dimensions

from pysec2pri.api import _resolve_parser_class
from pysec2pri.parsers.base import CONFIG_DIR

_CONFIG_PATHS = sorted(CONFIG_DIR.glob("*.yaml"))


def _mapping_set_specs() -> list[tuple[str, str, dict[str, Any]]]:
    """Return ``(config_id, kind, spec)`` for every mapping set of every config."""
    import yaml

    specs: list[tuple[str, str, dict[str, Any]]] = []
    for path in _CONFIG_PATHS:
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
        for kind, spec in (raw.get("mapping_sets") or {}).items():
            specs.append((path.stem, kind, spec))
    return specs


_SPECS = _mapping_set_specs()
_SPEC_IDS = [f"{cfg}-{kind}" for cfg, kind, _ in _SPECS]


@pytest.mark.parametrize("path", _CONFIG_PATHS, ids=lambda p: p.name)
def test_shipped_config_validates(path: Path) -> None:
    """Every config/*.yaml file in the package is schema-valid."""
    validate_config_file(path)


@pytest.mark.parametrize("path", _CONFIG_PATHS, ids=lambda p: p.name)
def test_parser_class_resolves(path: Path) -> None:
    """Every config's ``parser_class`` names an importable parser."""
    _resolve_parser_class(path.stem)


@pytest.mark.parametrize(("config_id", "kind", "spec"), _SPECS, ids=_SPEC_IDS)
def test_mapping_set_method_exists(config_id: str, kind: str, spec: dict[str, Any]) -> None:
    """Each mapping set's ``method`` is a real method of its parser."""
    parser_cls = _resolve_parser_class(config_id)
    method = spec["method"]
    assert hasattr(parser_cls, method), (
        f"{config_id}.yaml declares {kind}.method={method!r}, absent from {parser_cls.__name__}"
    )


@pytest.mark.parametrize(("config_id", "kind", "spec"), _SPECS, ids=_SPEC_IDS)
def test_mapping_set_inputs_bind_to_parameters(
    config_id: str, kind: str, spec: dict[str, Any]
) -> None:
    """Each ``inputs`` value names a parameter of the mapping set's method."""
    parser_cls = _resolve_parser_class(config_id)
    method = getattr(parser_cls, spec["method"])
    params = set(inspect.signature(method).parameters)
    for key, param in (spec.get("inputs") or {}).items():
        assert param in params, (
            f"{config_id}.yaml binds {kind}.inputs[{key!r}] to {param!r}, "
            f"not a parameter of {parser_cls.__name__}.{spec['method']}"
        )


def test_wikidata_query_keys_match_registries() -> None:
    """Wikidata's declared entity types line up with the registries they name.

    ``wikidata.yaml``'s ``queries`` keys *are* the entity types, so a key with
    no query text or no column mapping (or vice versa) is a drift.
    """
    from pysec2pri.parsers.wikidata import WikidataParser
    from pysec2pri.queries import QUERY_COLUMNS, WIKIDATA_QUERIES, WIKIDATA_TEST_QUERIES

    declared = set(WikidataParser.entity_types())
    assert declared, "wikidata.yaml declares no queries"
    assert declared == set(WIKIDATA_QUERIES), "config queries != WIKIDATA_QUERIES"
    assert declared == set(WIKIDATA_TEST_QUERIES), "config queries != WIKIDATA_TEST_QUERIES"
    assert declared == set(QUERY_COLUMNS), "config queries != QUERY_COLUMNS"


@pytest.mark.parametrize(("config_id", "kind", "spec"), _SPECS, ids=_SPEC_IDS)
def test_mapping_set_inputs_are_download_keys(
    config_id: str, kind: str, spec: dict[str, Any]
) -> None:
    """Each ``inputs`` key is downloadable: it appears in ``download_urls``."""
    import yaml

    raw = yaml.safe_load((CONFIG_DIR / f"{config_id}.yaml").read_text(encoding="utf-8"))
    keys = set(raw.get("download_urls") or {})
    for era in raw.get("distribution_eras") or []:
        keys |= set(era.get("download_urls") or {})
    for key in spec.get("inputs") or {}:
        assert key in keys, (
            f"{config_id}.yaml lists {kind}.inputs[{key!r}] but no download_urls entry provides it"
        )


def test_docs_list_every_source() -> None:
    """The docs' supported-databases table matches what the configs declare.

    A source someone added but never documented is invisible; one documented
    but removed is a lie. Either way the table has to track the configs.
    """
    from pysec2pri.api import sources

    index = (CONFIG_DIR.parent.parent.parent / "docs" / "source" / "index.rst").read_text(
        encoding="utf-8"
    )
    table = index.split("Supported databases")[1].split("Quick Start")[0]
    for source in sources():
        assert f"``{source}``" in table, f"{source!r} is not in the docs table"
    documented = set(re.findall(r"``([a-z_]+)``\n     - ``ids``", table))
    assert documented <= set(sources()), f"documented but gone: {documented - set(sources())}"


def _parser_for(config_id: str) -> Any:
    """Return an instance of *config_id*'s parser, at an arbitrary version."""
    return _resolve_parser_class(config_id)(version="1", show_progress=False)


def _dimensions(config_id: str) -> list[str]:
    """Return the names *config_id*'s config lists under ``products``."""
    return product_dimensions(get_datasource_config(config_id, config_package="pysec2pri.config"))


@pytest.mark.parametrize("path", _CONFIG_PATHS, ids=lambda p: p.name)
def test_every_product_dimension_reaches_the_iris(path: Path) -> None:
    """Each name in ``products`` must change the record namespace when it changes.

    A dimension the IRIs ignore means two different datasets of the same
    release collide on ``mapping_set_id`` and ``record_id``.
    """
    config_id = path.stem
    parser = _parser_for(config_id)
    dimensions = _dimensions(config_id)

    for name in dimensions:
        namespaces = set()
        for value in ("aaa", "bbb"):
            setattr(parser, name, value)
            namespaces.add(parser._record_namespace())
        assert len(namespaces) == 2, f"{config_id}: {name!r} does not reach the record namespace"


@pytest.mark.parametrize("path", _CONFIG_PATHS, ids=lambda p: p.name)
def test_product_dimensions_are_parser_arguments(path: Path) -> None:
    """Each name in ``products`` must be an argument one of the parser's methods takes."""
    config_id = path.stem
    cls = _resolve_parser_class(config_id)
    dimensions = _dimensions(config_id)

    accepted = {
        param
        for method in vars(cls).values()
        if callable(method)
        for param in inspect.signature(method).parameters
    }
    for name in dimensions:
        assert name in accepted, f"{config_id}: no parser method takes {name!r}"


@pytest.mark.parametrize("path", _CONFIG_PATHS, ids=lambda p: p.name)
def test_source_without_products_has_no_slug(path: Path) -> None:
    """A source declaring no ``products`` puts nothing but the version in its IRIs."""
    config_id = path.stem
    if _dimensions(config_id):
        pytest.skip("declares product dimensions")
    assert _parser_for(config_id)._product_slug() is None
