"""Secondary to primary identifier mapping.

This package provides tools for converting secondary (retired/withdrawn)
identifiers to primary (current) identifiers across various biological
databases including ChEBI, HMDB, HGNC, NCBI Gene, and UniProt.

Public names are imported lazily (PEP 562): ``import pysec2pri`` stays cheap and
only pulls the heavy ``sssom``/``linkml``/``polars`` stack when a symbol that
needs it is first accessed.
"""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

# name -> submodule providing it (resolved on first attribute access)
_LAZY_EXPORTS: dict[str, str] = {
    # Disambiguation context (symbol/id/xref) — from mapkgsutils
    "ContextSpec": "mapkgsutils.context",
    "DecisionRecord": "mapkgsutils.context",
    "XrefMapping": "mapkgsutils.context",
    "XrefRecord": "mapkgsutils.context",
    "load_xref_mapping": "mapkgsutils.context",
    # Mapping-set classes and sentinels
    "WITHDRAWN_ENTRY": "pysec2pri.parsers.base",
    "WITHDRAWN_ENTRY_LABEL": "pysec2pri.parsers.base",
    "AmbiguousMappingSet": "pysec2pri.parsers.base",
    "BaseMappingSet": "pysec2pri.parsers.base",
    "IdMappingSet": "pysec2pri.parsers.base",
    "LabelMappingSet": "pysec2pri.parsers.base",
    # ID resolution
    "build_alias_index": "pysec2pri.update",
    "build_label_lookup": "pysec2pri.update",
    "build_lookup": "pysec2pri.update",
    "build_primary_token_to_id": "pysec2pri.update",
    "resolve_ambiguous_with_hints": "pysec2pri.update",
    "update_ids": "pysec2pri.update",
    "update_labels": "pysec2pri.update",
    # High-level API (generate_*/load_*/write_*/resolve_*)
    "find_ambiguous": "pysec2pri.api",
    "list_versions": "pysec2pri.api",
    "generate_ids": "pysec2pri.api",
    "generate_labels": "pysec2pri.api",
    "generate_primary_ids": "pysec2pri.api",
    "generate_primary_labels": "pysec2pri.api",
    "sources": "pysec2pri.api",
    "supports_consolidate": "pysec2pri.api",
    "load_label_mapping": "pysec2pri.api",
    "load_mapping": "pysec2pri.api",
    "resolve_ids": "pysec2pri.api",
    "resolve_labels": "pysec2pri.api",
    "write_json": "pysec2pri.api",
    "write_label_sec2pri": "pysec2pri.api",
    "write_name2synonym": "pysec2pri.api",
    "write_output": "pysec2pri.api",
    "write_owl": "pysec2pri.api",
    "write_rdf": "pysec2pri.api",
    "write_sec2pri": "pysec2pri.api",
    "write_sssom": "pysec2pri.api",
}


def __getattr__(name: str) -> object:
    """Import and return a public symbol on first access (PEP 562)."""
    module_name = _LAZY_EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(importlib.import_module(module_name), name)
    globals()[name] = value  # cache so subsequent access skips __getattr__
    return value


def __dir__() -> list[str]:
    """Expose lazily-available names to ``dir()`` and tab completion."""
    return sorted(__all__)


if TYPE_CHECKING:
    # Imports for type checkers
    from mapkgsutils.context import (
        ContextSpec,
        DecisionRecord,
        XrefMapping,
        XrefRecord,
        load_xref_mapping,
    )

    from pysec2pri.api import (
        find_ambiguous,
        generate_ids,
        generate_labels,
        generate_primary_ids,
        generate_primary_labels,
        list_versions,
        load_label_mapping,
        load_mapping,
        resolve_ids,
        resolve_labels,
        sources,
        supports_consolidate,
        write_json,
        write_label_sec2pri,
        write_name2synonym,
        write_output,
        write_owl,
        write_rdf,
        write_sec2pri,
        write_sssom,
    )
    from pysec2pri.parsers.base import (
        WITHDRAWN_ENTRY,
        WITHDRAWN_ENTRY_LABEL,
        AmbiguousMappingSet,
        BaseMappingSet,
        IdMappingSet,
        LabelMappingSet,
    )
    from pysec2pri.update import (
        build_alias_index,
        build_label_lookup,
        build_lookup,
        build_primary_token_to_id,
        resolve_ambiguous_with_hints,
        update_ids,
        update_labels,
    )

__all__ = [
    # Sentinel values
    "WITHDRAWN_ENTRY",
    "WITHDRAWN_ENTRY_LABEL",
    # Mapping Set classes
    "AmbiguousMappingSet",
    "BaseMappingSet",
    # Disambiguation context (symbol/id/xref)
    "ContextSpec",
    "DecisionRecord",
    "IdMappingSet",
    "LabelMappingSet",
    "XrefMapping",
    "XrefRecord",
    # ID resolution
    "build_alias_index",
    "build_label_lookup",
    "build_lookup",
    "build_primary_token_to_id",
    "find_ambiguous",
    "generate_ids",
    "generate_labels",
    "generate_primary_ids",
    "generate_primary_labels",
    # generate_* (download + parse in one call)
    "list_versions",
    "load_label_mapping",
    "load_mapping",
    "load_xref_mapping",
    "resolve_ambiguous_with_hints",
    "resolve_ids",
    "resolve_labels",
    "sources",
    "supports_consolidate",
    "update_ids",
    "update_labels",
    # Export functions
    "write_json",
    "write_label_sec2pri",
    "write_name2synonym",
    "write_output",
    "write_owl",
    "write_rdf",
    "write_sec2pri",
    "write_sssom",
]
