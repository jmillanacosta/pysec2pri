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
    "build_alias_index": "pysec2pri.update_ids",
    "build_label_lookup": "pysec2pri.update_ids",
    "build_lookup": "pysec2pri.update_ids",
    "build_primary_token_to_id": "pysec2pri.update_ids",
    "resolve_ambiguous_with_hints": "pysec2pri.update_ids",
    "update_ids": "pysec2pri.update_ids",
    "update_labels": "pysec2pri.update_ids",
    # High-level API (generate_*/load_*/write_*/resolve_*/crosswalk)
    "crosswalk": "pysec2pri.api",
    "find_ambiguous": "pysec2pri.api",
    "generate_chebi": "pysec2pri.api",
    "generate_chebi_synonyms": "pysec2pri.api",
    "generate_ensembl": "pysec2pri.api",
    "generate_ensembl_labels": "pysec2pri.api",
    "generate_hgnc": "pysec2pri.api",
    "generate_hgnc_labels": "pysec2pri.api",
    "generate_hgnc_primary_ids": "pysec2pri.api",
    "generate_hmdb": "pysec2pri.api",
    "generate_hmdb_proteins": "pysec2pri.api",
    "generate_ncbi": "pysec2pri.api",
    "generate_ncbi_labels": "pysec2pri.api",
    "generate_uniprot": "pysec2pri.api",
    "generate_vgnc": "pysec2pri.api",
    "generate_vgnc_labels": "pysec2pri.api",
    "generate_wikidata": "pysec2pri.api",
    "generate_wikidata_labels": "pysec2pri.api",
    "list_versions": "pysec2pri.api",
    "load_chebi": "pysec2pri.api",
    "load_ensembl": "pysec2pri.api",
    "load_hgnc": "pysec2pri.api",
    "load_hmdb": "pysec2pri.api",
    "load_hmdb_proteins": "pysec2pri.api",
    "load_label_mapping": "pysec2pri.api",
    "load_mapping": "pysec2pri.api",
    "load_ncbi": "pysec2pri.api",
    "load_uniprot": "pysec2pri.api",
    "load_vgnc": "pysec2pri.api",
    "load_wikidata": "pysec2pri.api",
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
        crosswalk,
        find_ambiguous,
        generate_chebi,
        generate_chebi_synonyms,
        generate_ensembl,
        generate_ensembl_labels,
        generate_hgnc,
        generate_hgnc_labels,
        generate_hgnc_primary_ids,
        generate_hmdb,
        generate_hmdb_proteins,
        generate_ncbi,
        generate_ncbi_labels,
        generate_uniprot,
        generate_vgnc,
        generate_vgnc_labels,
        generate_wikidata,
        generate_wikidata_labels,
        list_versions,
        load_chebi,
        load_ensembl,
        load_hgnc,
        load_hmdb,
        load_hmdb_proteins,
        load_label_mapping,
        load_mapping,
        load_ncbi,
        load_uniprot,
        load_vgnc,
        load_wikidata,
        resolve_ids,
        resolve_labels,
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
    from pysec2pri.update_ids import (
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
    "crosswalk",
    "find_ambiguous",
    # generate_* (download + parse in one call)
    "generate_chebi",
    "generate_chebi_synonyms",
    "generate_ensembl",
    "generate_ensembl_labels",
    "generate_hgnc",
    "generate_hgnc_labels",
    "generate_hgnc_primary_ids",
    "generate_hmdb",
    "generate_hmdb_proteins",
    "generate_ncbi",
    "generate_ncbi_labels",
    "generate_uniprot",
    "generate_vgnc",
    "generate_vgnc_labels",
    "generate_wikidata",
    "generate_wikidata_labels",
    "list_versions",
    "load_chebi",
    "load_ensembl",
    "load_hgnc",
    "load_hmdb",
    "load_hmdb_proteins",
    "load_label_mapping",
    "load_mapping",
    "load_ncbi",
    "load_uniprot",
    "load_vgnc",
    "load_wikidata",
    "load_xref_mapping",
    "resolve_ambiguous_with_hints",
    "resolve_ids",
    "resolve_labels",
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
