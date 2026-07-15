"""Resolve secondary identifiers and previous labels to their primary form.

Applies mod:`mapkgsutils.resolve`.

Takes a single string delimited by ``|``, ``,``, ``;``, or whitespace::

    from pysec2pri import generate_ids, update_ids

    ms = generate_ids("hgnc")
    update_ids("HGNC:1234|HGNC:5678", ms)

DataFrame, annotating one or more columns::

    update_ids(df, ms, at="gene_id")
    update_ids(df, ms, at=["gene_id", "alt_id"])
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Union, overload

from mapkgsutils import resolve as _resolve
from mapkgsutils.context import ContextSpec, XrefMapping, XrefRecord, load_xref_mapping
from mapkgsutils.resolve import build_primary_token_to_id, resolve_ambiguous_with_hints

if TYPE_CHECKING:
    from pathlib import Path

    import pandas as pd

    from pysec2pri.parsers.base import BaseMappingSet

IdsInput = Union[str, list[str], "pd.DataFrame"]

# Predicate carrying the deprecation relation resolved here (``term replaced by``).
_DEPRECATION_PREDICATE = "IAO:0100001"


def build_alias_index(mapping_set: BaseMappingSet) -> dict[str, list[str]]:
    """Return ``{primary_id: [alias label]}``, excluding deprecation edges."""
    return _resolve.build_alias_index(mapping_set, exclude_predicates={_DEPRECATION_PREDICATE})


def build_lookup(mapping_set: BaseMappingSet) -> dict[str, str]:
    """Return ``{secondary_id: primary_id}``."""
    return _resolve.build_lookup(mapping_set, on="id")


def build_label_lookup(mapping_set: BaseMappingSet) -> dict[str, str]:
    """Return ``{previous_label: current_label}``."""
    return _resolve.build_lookup(mapping_set, on="label")


def build_ambiguous_set(mapping_set: BaseMappingSet) -> set[str]:
    """Return secondary IDs that are also current primary IDs."""
    return _resolve.build_ambiguous(mapping_set, on="id")


def build_ambiguous_labels_set(mapping_set: BaseMappingSet) -> set[str]:
    """Return previous labels that are also current primary labels."""
    return _resolve.build_ambiguous(mapping_set, on="label")


# docstr-coverage:excused `overload`
@overload
def update_ids(
    ids: str,
    mapping_set: BaseMappingSet,
    *,
    at: None = ...,
    suffix: str = ...,
    lookup: dict[str, str] | None = ...,
    ambiguous: set[str] | None = ...,
) -> dict[str, str]: ...


# docstr-coverage:excused `overload`
@overload
def update_ids(
    ids: list[str],
    mapping_set: BaseMappingSet,
    *,
    at: None = ...,
    suffix: str = ...,
    lookup: dict[str, str] | None = ...,
    ambiguous: set[str] | None = ...,
) -> dict[str, str]: ...


# docstr-coverage:excused `overload`
@overload
def update_ids(
    ids: pd.DataFrame,
    mapping_set: BaseMappingSet,
    *,
    at: str | list[str],
    suffix: str = ...,
    lookup: dict[str, str] | None = ...,
    ambiguous: set[str] | None = ...,
    synonyms: str | list[str] | None = ...,
    label_mapping_set: BaseMappingSet | None = ...,
    xref: str | None = ...,
    xref_mapping: XrefMapping | None = ...,
    xref_predicates: set[str] | None = ...,
    report_path: Path | str | None = ...,
    context: ContextSpec | list[ContextSpec] | None = ...,
) -> pd.DataFrame: ...


def update_ids(
    ids: IdsInput,
    mapping_set: BaseMappingSet,
    *,
    at: str | list[str] | None = None,
    suffix: str = "_primary",
    lookup: dict[str, str] | None = None,
    ambiguous: set[str] | None = None,
    synonyms: str | list[str] | None = None,
    label_mapping_set: BaseMappingSet | None = None,
    xref: str | None = None,
    xref_mapping: XrefMapping | None = None,
    xref_predicates: set[str] | None = None,
    report_path: Path | str | None = None,
    context: ContextSpec | list[ContextSpec] | None = None,
) -> dict[str, str] | pd.DataFrame:
    """Resolve secondary identifiers in *ids* to primary identifiers via *mapping_set*.

    *label_mapping_set* supplies the alias index for ``synonyms`` hint
    resolution. See :func:`mapkgsutils.resolve.resolve`.
    """
    return _resolve.resolve(  # type: ignore[call-overload, misc]
        ids,
        mapping_set,
        on="id",
        at=at,
        suffix=suffix,
        lookup=lookup,
        ambiguous=ambiguous,
        synonyms=synonyms,
        alias_mapping_set=label_mapping_set,
        xref=xref,
        xref_mapping=xref_mapping,
        xref_predicates=xref_predicates,
        relation_predicates={_DEPRECATION_PREDICATE},
        report_path=report_path,
        context=context,
    )


# docstr-coverage:excused `overload`
@overload
def update_labels(
    labels: str,
    mapping_set: BaseMappingSet,
    *,
    at: None = ...,
    suffix: str = ...,
    lookup: dict[str, str] | None = ...,
    ambiguous: set[str] | None = ...,
    synonyms: str | list[str] | None = ...,
) -> dict[str, str]: ...


# docstr-coverage:excused `overload`
@overload
def update_labels(
    labels: list[str],
    mapping_set: BaseMappingSet,
    *,
    at: None = ...,
    suffix: str = ...,
    lookup: dict[str, str] | None = ...,
    ambiguous: set[str] | None = ...,
    synonyms: str | list[str] | None = ...,
) -> dict[str, str]: ...


# docstr-coverage:excused `overload`
@overload
def update_labels(
    labels: pd.DataFrame,
    mapping_set: BaseMappingSet,
    *,
    at: str | list[str],
    suffix: str = ...,
    lookup: dict[str, str] | None = ...,
    ambiguous: set[str] | None = ...,
    synonyms: str | list[str] | None = ...,
    xref: str | None = ...,
    xref_mapping: XrefMapping | None = ...,
    xref_predicates: set[str] | None = ...,
    report_path: Path | str | None = ...,
    context: ContextSpec | list[ContextSpec] | None = ...,
) -> pd.DataFrame: ...


def update_labels(
    labels: IdsInput,
    mapping_set: BaseMappingSet,
    *,
    at: str | list[str] | None = None,
    suffix: str = "_current",
    lookup: dict[str, str] | None = None,
    ambiguous: set[str] | None = None,
    synonyms: str | list[str] | None = None,
    xref: str | None = None,
    xref_mapping: XrefMapping | None = None,
    xref_predicates: set[str] | None = None,
    report_path: Path | str | None = None,
    context: ContextSpec | list[ContextSpec] | None = None,
) -> dict[str, str] | pd.DataFrame:
    """Resolve previous/alias labels in *labels* to current labels via *mapping_set*.

    The alias index for ``synonyms`` hint resolution comes from *mapping_set*.
    See :func:`mapkgsutils.resolve.resolve`.
    """
    return _resolve.resolve(  # type: ignore[call-overload, misc]
        labels,
        mapping_set,
        on="label",
        at=at,
        suffix=suffix,
        lookup=lookup,
        ambiguous=ambiguous,
        synonyms=synonyms,
        xref=xref,
        xref_mapping=xref_mapping,
        xref_predicates=xref_predicates,
        relation_predicates={_DEPRECATION_PREDICATE},
        report_path=report_path,
        context=context,
    )


__all__ = [
    "ContextSpec",
    "XrefMapping",
    "XrefRecord",
    "build_alias_index",
    "build_ambiguous_labels_set",
    "build_ambiguous_set",
    "build_label_lookup",
    "build_lookup",
    "build_primary_token_to_id",
    "load_xref_mapping",
    "resolve_ambiguous_with_hints",
    "update_ids",
    "update_labels",
]
