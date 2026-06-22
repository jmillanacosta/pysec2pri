"""Resolve secondary identifiers to primary identifiers using a MappingSet.

Typical usage
-------------
Single string (possibly separated by commas/semicolons/pipes/whitespace)::

    from pysec2pri import generate_hgnc
    from pysec2pri.update_ids import update_ids

    ms = generate_hgnc()
    update_ids("HGNC:1234|HGNC:5678", ms)
    # {'HGNC:1234': 'HGNC:9999', 'HGNC:5678': 'HGNC:5678'}

List of strings::

    update_ids(["HGNC:1234", "HGNC:5678"], ms)

Pandas DataFrame, annotate one or more columns::

    import pandas as pd

    df = pd.DataFrame({"gene_id": ["HGNC:1234", "HGNC:5678"]})
    update_ids(df, ms, at="gene_id")
    # returns df with an extra column "gene_id_primary"

    # Multiple columns at once:
    update_ids(df, ms, at=["gene_id", "alt_id"])
    # returns df with "gene_id_primary" and "alt_id_primary" columns added

Notes
-----
* Identifiers that are not found in the mapping set are returned/kept
  as-is.
* Identifiers separated by common delimiters (``|``, ``,``, ``;``,
  whitespace) inside a single string are each looked up individually.
* The mapping look-up is done once against the full set of unique IDs
  to avoid repeated scans of large mapping sets.
* **Ambiguous identifiers**, those that appear both as a secondary ID
  in the mapping set and as a current primary ID, are left blank in
  the resolved output.  A warning is emitted listing every ambiguous
  token so the user can resolve them manually.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING, Union, overload

if TYPE_CHECKING:
    from pathlib import Path

    import pandas as pd

from pysec2pri.context import (
    ContextSpec,
    DecisionRecord,
    XrefMapping,
    XrefRecord,
    load_xref_mapping,
    resolve_ambiguous_with_xref,
    write_decision_log,
)
from pysec2pri.logging import logger
from pysec2pri.parsers.base import Sec2PriMappingSet

_IAO_DEPRECATION = "IAO:0100001"

# Separator pattern: pipe, comma, semicolon, or whitespace
_SEP = re.compile(r"[|,;\s]+")

# Type alias for the flexible input
IdsInput = Union[str, list[str], "pd.DataFrame"]

# Sentinel for "this ID is ambiguous, do not replace"
_AMBIGUOUS = object()


# Helpers


def _build_lookup(mapping_set: Sec2PriMappingSet) -> dict[str, str]:
    """Return a ``{subject_id: object_id}`` dict from *mapping_set*.

    Only mappings with a non-empty ``subject_id`` are included.  If a
    secondary ID maps to the withdrawn sentinel it is kept as-is so
    callers can decide what to do with it.
    """
    lookup: dict[str, str] = {}
    for m in mapping_set.mappings or []:
        sec = str(getattr(m, "subject_id", None) or "")
        pri = str(getattr(m, "object_id", None) or "")
        if sec:
            lookup[sec] = pri
    return lookup


# Synonym-hint helpers (public)


def build_alias_index(mapping_set: Sec2PriMappingSet) -> dict[str, list[str]]:
    """Return ``{object_id: [subject_labels linked via non-IAO predicates]}``.

    Builds an index of all *non-deprecation* alias mappings in a
    :class:`~pysec2pri.parsers.base.LabelMappingSet`.  Only entries whose
    ``predicate_id`` is **not** ``IAO:0100001`` are included; deprecation
    (``IAO:0100001`` / "term replaced by") mappings are deliberately excluded
    because they express history, not active aliasing.

    This index is used by :func:`resolve_ambiguous_with_hints` to confirm
    whether a user-supplied alias belongs to the secondary mapping's target
    (confirming secondary usage) or to the entity's own primary entry
    (confirming primary usage).

    Args:
        mapping_set: A :class:`~pysec2pri.parsers.base.LabelMappingSet`
            (e.g. the result of ``generate_hgnc_labels()``).

    Returns:
        Dict mapping each ``object_id`` to the list of ``subject_label``
        values that point to it via a non-IAO predicate.
    """
    index: dict[str, list[str]] = {}
    for m in mapping_set.mappings or []:
        pred = str(getattr(m, "predicate_id", None) or "")
        if pred == _IAO_DEPRECATION:
            continue
        obj_id = str(getattr(m, "object_id", None) or "")
        subj_lbl = str(getattr(m, "subject_label", None) or "")
        if obj_id and subj_lbl:
            index.setdefault(obj_id, []).append(subj_lbl)
    return index


def build_primary_token_to_id(mapping_set: Sec2PriMappingSet) -> dict[str, str]:
    """Return ``{primary_label: primary_id}`` from a label mapping set.

    Collects every ``(object_label, object_id)`` pair seen in the mappings
    and, where available, the ``_primary_labels`` store.  Useful for
    translating a primary label string into its CURIE so that
    :func:`build_alias_index` (keyed by object_id) can be looked up.

    Args:
        mapping_set: A :class:`~pysec2pri.parsers.base.LabelMappingSet`.

    Returns:
        Dict ``{primary_label: primary_id}``.
    """
    result: dict[str, str] = {}
    for m in mapping_set.mappings or []:
        lbl = str(getattr(m, "object_label", None) or "")
        oid = str(getattr(m, "object_id", None) or "")
        if lbl and oid:
            result[lbl] = oid
    stored: dict[str, set[str]] | None = getattr(mapping_set, "_primary_labels", None)
    if isinstance(stored, dict):
        for sym, ids in stored.items():
            if sym and ids and sym not in result:
                result[sym] = next(iter(ids))
    return result


def resolve_ambiguous_with_hints(
    ambiguous_token: str,
    user_aliases: list[str],
    lkp: dict[str, str],
    alias_index: dict[str, list[str]],
    token_to_id: dict[str, str] | None = None,
) -> tuple[str, str | None]:
    """Attempt to resolve an ambiguous label or ID using user-provided alias hints.

    An ambiguous token appears both as a current primary entry **and** as a
    secondary (subject) in a mapping that points to a different primary.

    Two resolution cases are checked:

    1. **Secondary usage**: at least one of the user-supplied aliases matches
       the target token itself (its primary label or primary ID), or appears
       among the non-IAO aliases of the mapping's *target*
       (``lkp[ambiguous_token]``).  This confirms the token is being used as
       a secondary alias of the target, returns ``(target_token, target_id)``.

    2. **Primary usage**: at least one of the user-supplied aliases appears
       among the non-IAO aliases of the token's *own* primary entry.  This
       confirms the token is being used as a standalone primary,
       returns ``(ambiguous_token, own_id)``.

    If neither case applies the ambiguity cannot be resolved and
    ``("", None)`` is returned.

    Args:
        ambiguous_token:
            The label or ID that is both primary and secondary.
        user_aliases:
            Alias strings provided by the user (e.g. from a ``known_alias``
            column) to help determine which entity is actually meant.
        lkp:
            ``{secondary_token: resolved_token}`` lookup (from
            :func:`build_label_lookup` or :func:`build_lookup`).
        alias_index:
            ``{primary_id: [non-IAO alias strings]}`` built via
            :func:`build_alias_index`.
        token_to_id:
            ``{primary_token: primary_id}``.  When ``None`` the token is
            treated as its own ID, which is appropriate for ID mapping sets
            where the token already is a CURIE (e.g. ``"HGNC:53564"``).

    Returns:
        A ``(resolved_token, resolved_id)`` tuple:

        * Secondary case -> ``(target_token, target_id)``
        * Primary case   -> ``(ambiguous_token, own_id)``
        * Unresolvable   -> ``("", None)``
    """
    if not user_aliases:
        return ("", None)

    target_token = lkp.get(ambiguous_token, "")

    def _id(tok: str) -> str:
        return (token_to_id or {}).get(tok, tok)

    target_id = _id(target_token) if target_token else ""

    user_set = set(user_aliases)

    # Case 1: secondary usage.
    # The user alias matches the target's primary label, its primary ID,
    # *or* one of its non-IAO alias strings.  Matching the primary label
    # directly covers the common case where the user writes the current
    # label (e.g. "RNF141") rather than a secondary alias.
    if target_token and (
        target_token in user_set
        or (target_id and target_id in user_set)
        or (target_id and user_set & set(alias_index.get(target_id, [])))
    ):
        return (target_token, target_id if target_id else None)

    # Case 2: primary usage. User aliases overlap own primary's alias list
    own_id = _id(ambiguous_token)
    if own_id and user_set & set(alias_index.get(own_id, [])):
        return (ambiguous_token, own_id if own_id else None)

    return ("", None)


def build_ambiguous_set(mapping_set: Sec2PriMappingSet) -> set[str]:
    """Return the set of *ambiguous* subject IDs in *mapping_set*.

    An identifier is ambiguous when it appears both as a ``subject_id``
    (i.e. a secondary/previous term) and as a current primary
    identifier, either in the explicitly stored ``_primary_ids`` set or
    among the ``object_id`` values of the mappings.

    When such overlap exists a naïve replacement could silently corrupt
    references that already use the current entity, so the resolver
    intentionally leaves those cells blank.

    Args:
        mapping_set: Any :class:`~pysec2pri.parsers.base.Sec2PriMappingSet`.

    Returns:
        Set of ID strings that are both secondary and primary.  Empty set
        when no ambiguity is detected.
    """
    stored: set[str] = (
        object.__getattribute__(mapping_set, "_primary_ids")
        if hasattr(mapping_set, "_primary_ids")
        else set()
    )
    primary_ids: set[str] = stored or {
        str(getattr(m, "object_id", None) or "") for m in (mapping_set.mappings or [])
    } - {""}
    subject_ids: set[str] = {
        str(getattr(m, "subject_id", None) or "") for m in (mapping_set.mappings or [])
    } - {""}
    return subject_ids & primary_ids


def build_ambiguous_labels_set(mapping_set: Sec2PriMappingSet) -> set[str]:
    """Return the set of *ambiguous* subject labels in *mapping_set*.

    Analogous to :func:`build_ambiguous_set` but operates on
    ``subject_label`` / ``object_label`` (label) mappings.

    Args:
        mapping_set: Any :class:`~pysec2pri.parsers.base.LabelMappingSet`.

    Returns:
        Set of label strings that are both secondary and primary.  Empty
        set when no ambiguity is detected.
    """
    stored: dict[str, set[str]] | set[str] = (
        object.__getattribute__(mapping_set, "_primary_labels")
        if hasattr(mapping_set, "_primary_labels")
        else set()
    )
    # _primary_labels is stored as dict[label, set[ids]]; extract the key set.
    primary_labels: set[str] = (set(stored.keys()) if isinstance(stored, dict) else stored) or {
        str(getattr(m, "object_label", None) or "") for m in (mapping_set.mappings or [])
    } - {""}
    subject_labels: set[str] = {
        str(getattr(m, "subject_label", None) or "") for m in (mapping_set.mappings or [])
    } - {""}
    return subject_labels & primary_labels


def _warn_ambiguous(ambiguous_found: set[str], kind: str = "ID") -> None:
    """Emit a structured warning for every ambiguous token encountered.

    Called after a resolver pass to report all tokens that were left blank
    because their resolution was ambiguous.

    Args:
        ambiguous_found: Set of token strings that were ambiguous.
        kind: Human-readable label for the token type (``"ID"`` or
            ``"label"``), used in the warning message.
    """
    if not ambiguous_found:
        return
    count = len(ambiguous_found)
    listed = ", ".join(sorted(ambiguous_found))
    logger.warning(
        "%d ambiguous %s(s) were left blank in the resolved output because "
        "the same %s appears both as a current primary %s in this datasource "
        "AND as a secondary/previous %s that maps to a different entry. "
        "Automatic replacement would be unsafe, please resolve these manually. "
        "Ambiguous %s(s): %s",
        count,
        kind,
        kind.lower(),
        kind.lower(),
        kind.lower(),
        kind,
        listed,
    )


def _split_ids(raw: str) -> list[str]:
    """Split *raw* on common delimiters; strips empty tokens."""
    return [tok for tok in _SEP.split(raw.strip()) if tok]


def _resolve_tokens(
    tokens: list[str],
    lookup: dict[str, str],
    ambiguous: set[str],
    ambiguous_found: set[str],
) -> list[str]:
    """Map each token to its primary ID, blank for ambiguous, self for unknown."""
    result: list[str] = []
    for tok in tokens:
        if tok in ambiguous:
            ambiguous_found.add(tok)
            result.append("")
        else:
            result.append(lookup.get(tok, tok))
    return result


def _resolve_string(
    raw: str,
    lookup: dict[str, str],
    ambiguous: set[str],
    ambiguous_found: set[str],
) -> str:
    """Resolve all IDs inside *raw*, rejoining with the original separator."""
    sep_match = re.search(r"[|,;\s]", raw)
    sep = sep_match.group(0) if sep_match else ""
    tokens = _split_ids(raw)
    resolved = _resolve_tokens(tokens, lookup, ambiguous, ambiguous_found)
    return sep.join(resolved)


# Public API


def build_lookup(mapping_set: Sec2PriMappingSet) -> dict[str, str]:
    """Return a ``{secondary_id: primary_id}`` dictionary.

    Useful when you want to apply the look-up yourself or cache it for
    repeated calls.

    Args:
        mapping_set: A :class:`~pysec2pri.parsers.base.Sec2PriMappingSet`
            (e.g. the object returned by ``generate_hgnc()``).

    Returns:
        Dictionary mapping every secondary ID to its current primary ID.
    """
    return _build_lookup(mapping_set)


# --- dispatch helpers (keep public functions below C901 threshold) ---


def _update_str(
    ids: str,
    lkp: dict[str, str],
    amb: set[str],
    kind: str,
) -> dict[str, str]:
    tokens = _split_ids(ids)
    unique = dict.fromkeys(tokens)
    ambiguous_found: set[str] = set()
    result: dict[str, str] = {}
    for tok in unique:
        if tok in amb:
            ambiguous_found.add(tok)
            result[tok] = ""
        else:
            result[tok] = lkp.get(tok, tok)
    _warn_ambiguous(ambiguous_found, kind=kind)
    return result


def _update_list(
    ids: list[str],
    lkp: dict[str, str],
    amb: set[str],
    kind: str,
) -> dict[str, str]:
    unique: dict[str, None] = {}
    for item in ids:
        for tok in _split_ids(item):
            unique[tok] = None
    ambiguous_found: set[str] = set()
    result: dict[str, str] = {}
    for tok in unique:
        if tok in amb:
            ambiguous_found.add(tok)
            result[tok] = ""
        else:
            result[tok] = lkp.get(tok, tok)
    _warn_ambiguous(ambiguous_found, kind=kind)
    return result


def _update_dataframe(
    df: pd.DataFrame,
    at: str | list[str] | None,
    suffix: str,
    lkp: dict[str, str],
    amb: set[str],
    kind: str,
    col_label: str,
) -> pd.DataFrame:
    import pandas as pd

    if not isinstance(df, pd.DataFrame):
        raise TypeError(
            f"'{col_label}' must be a str, list[str], or pandas.DataFrame, "
            f"got {type(df).__name__!r}."
        )
    if at is None:
        raise ValueError(
            f"When '{col_label}' is a DataFrame you must specify 'at' "
            "(column name or list of names)."
        )
    columns: list[str] = [at] if isinstance(at, str) else list(at)
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise KeyError(f"Column(s) not found in DataFrame: {missing}")
    ambiguous_found: set[str] = set()
    result = df.copy()
    for col in columns:
        result[col + suffix] = (
            result[col]
            .astype(str)
            .map(
                lambda cell, _lkp=lkp, _amb=amb, _af=ambiguous_found: _resolve_string(
                    cell, _lkp, _amb, _af
                )
            )
        )
    _warn_ambiguous(ambiguous_found, kind=kind)
    return result


def _raw_context_value(df: pd.DataFrame, idx: int, column: str) -> str:
    """Return the row's evidence value for *column*, or ``""`` when absent/NaN."""
    if column not in df.columns:
        return ""
    raw = str(df.iloc[idx][column])
    return raw if raw.strip() not in {"", "nan"} else ""


def _hint_decision(
    ambiguous_token: str,
    user_aliases: list[str],
    lkp: dict[str, str],
    alias_index: dict[str, list[str]],
    token_to_id: dict[str, str] | None,
    stage: str,
) -> tuple[str, str | None, DecisionRecord]:
    """Resolve via :func:`resolve_ambiguous_with_hints`, with a decision record.

    Thin wrapper so the ``label``/``id`` context kinds produce the same
    :class:`~pysec2pri.context.DecisionRecord` shape as
    :func:`~pysec2pri.context.resolve_ambiguous_with_xref`, for a unified
    decision log.
    """
    token = "|".join(user_aliases)
    resolved_tok, resolved_id = resolve_ambiguous_with_hints(
        ambiguous_token, user_aliases, lkp, alias_index, token_to_id
    )
    if not resolved_tok:
        reason = (
            "no alias hints given" if not user_aliases else "no alias hint matched either identity"
        )
        return "", None, DecisionRecord(stage, token, None, None, False, reason)
    reason = (
        "alias hint matched secondary target"
        if resolved_tok != ambiguous_token
        else "alias hint matched own primary identity"
    )
    return resolved_tok, resolved_id, DecisionRecord(stage, token, None, resolved_tok, True, reason)


def _resolve_cell_with_context(
    cell: str,
    lkp: dict[str, str],
    amb: set[str],
    specs: list[ContextSpec],
    raw_values: list[str],
    xref_indices: list[dict[str, list[XrefRecord]] | None],
    alias_index: dict[str, list[str]],
    token_to_id: dict[str, str] | None,
    ambiguous_found: set[str],
    decisions: list[DecisionRecord],
) -> tuple[str, str | None]:
    """Resolve one cell, trying each of *specs* in order on ambiguous tokens.

    The first :class:`~pysec2pri.context.ContextSpec` that resolves a given
    ambiguous token wins; every attempt (successful or not) is appended to
    *decisions* for an auditable log.

    Returns a ``(resolved_value, resolved_id)`` tuple where *resolved_id* is
    the primary ID of the resolved token, and ``None`` when unresolved.
    """
    sep_match = re.search(r"[|,;\s]", cell)
    sep = sep_match.group(0) if sep_match else ""
    tokens = _split_ids(cell)
    resolved: list[str] = []
    resolved_id: str | None = None
    for tok in tokens:
        if tok in amb:
            hint_tok = ""
            hint_id: str | None = None
            for spec, raw_val, xref_index in zip(specs, raw_values, xref_indices, strict=True):
                if spec.kind == "xref":
                    if xref_index is None:
                        continue
                    hint_tok, hint_id, decision = resolve_ambiguous_with_xref(
                        tok,
                        raw_val,
                        lkp,
                        xref_index,
                        token_to_id,
                        accepted_predicates=spec.predicates,
                    )
                else:
                    user_aliases = _split_ids(raw_val) if raw_val else []
                    hint_tok, hint_id, decision = _hint_decision(
                        tok, user_aliases, lkp, alias_index, token_to_id, stage=spec.kind
                    )
                decisions.append(decision)
                if hint_tok:
                    break
            if hint_tok:
                resolved.append(hint_tok)
                resolved_id = hint_id
            else:
                ambiguous_found.add(tok)
                resolved.append("")
        else:
            resolved_tok = lkp.get(tok, tok)
            resolved.append(resolved_tok)
            if resolved_id is None:
                # In ID mode (token_to_id=None) the resolved token IS the ID;
                # in label mode look up the corresponding primary ID.
                resolved_id = (token_to_id or {}).get(resolved_tok) or resolved_tok
    return sep.join(resolved), resolved_id


def _update_dataframe_with_context(
    df: pd.DataFrame,
    at: str | list[str] | None,
    suffix: str,
    lkp: dict[str, str],
    amb: set[str],
    kind: str,
    col_label: str,
    specs: list[ContextSpec],
    alias_index: dict[str, list[str]],
    token_to_id: dict[str, str] | None,
    report_path: Path | str | None,
) -> pd.DataFrame:
    """Like :func:`_update_dataframe` but resolves ambiguous cells using *specs*.

    Each :class:`~pysec2pri.context.ContextSpec` names a column carrying
    per-row evidence (an alias string for ``"label"``/``"id"``, or a
    cross-reference token for ``"xref"``). When *report_path* is given, every
    disambiguation attempt is written there as a TSV decision log.
    """
    import pandas as pd

    if not isinstance(df, pd.DataFrame):
        raise TypeError(
            f"'{col_label}' must be a str, list[str], or pandas.DataFrame, "
            f"got {type(df).__name__!r}."
        )
    if at is None:
        raise ValueError(
            f"When '{col_label}' is a DataFrame you must specify 'at' "
            "(column name or list of names)."
        )
    columns: list[str] = [at] if isinstance(at, str) else list(at)
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise KeyError(f"Column(s) not found in DataFrame: {missing}")

    xref_indices: list[dict[str, list[XrefRecord]] | None] = [
        spec.xref_mapping.by_subject() if spec.kind == "xref" and spec.xref_mapping else None
        for spec in specs
    ]

    result = df.copy()
    decisions: list[DecisionRecord] = []
    for col in columns:
        ambiguous_found: set[str] = set()
        new_values: list[str] = []
        new_ids: list[str | None] = []
        for idx in range(len(result)):
            cell = str(result.iloc[idx][col])
            raw_values = [_raw_context_value(result, idx, spec.column) for spec in specs]
            val, rid = _resolve_cell_with_context(
                cell,
                lkp,
                amb,
                specs,
                raw_values,
                xref_indices,
                alias_index,
                token_to_id,
                ambiguous_found,
                decisions,
            )
            new_values.append(val)
            new_ids.append(rid)
        result[col + suffix] = new_values
        result[col + suffix + "_id"] = new_ids
        _warn_ambiguous(ambiguous_found, kind=kind)

    if report_path is not None:
        write_decision_log(decisions, report_path)

    return result


def _update_dataframe_with_synonyms(
    df: pd.DataFrame,
    at: str | list[str] | None,
    suffix: str,
    lkp: dict[str, str],
    amb: set[str],
    kind: str,
    col_label: str,
    synonyms_col: str,
    alias_index: dict[str, list[str]],
    token_to_id: dict[str, str] | None,
) -> pd.DataFrame:
    """Like :func:`_update_dataframe_with_context` for a single ``label`` hint column."""
    spec = ContextSpec(kind="label", column=synonyms_col)
    return _update_dataframe_with_context(
        df,
        at,
        suffix,
        lkp,
        amb,
        kind,
        col_label,
        [spec],
        alias_index,
        token_to_id,
        report_path=None,
    )


def _build_context_specs(
    synonyms: str | list[str] | None,
    xref: str | None,
    xref_mapping: XrefMapping | None,
    xref_predicates: set[str] | None,
    context: ContextSpec | list[ContextSpec] | None,
) -> list[ContextSpec]:
    """Assemble the ``synonyms=``/``xref=``/``context=`` convenience kwargs into specs.

    *specs* are tried in order per ambiguous cell; the first one that
    resolves a given token wins (see :func:`_resolve_cell_with_context`).
    """
    specs: list[ContextSpec] = []
    if isinstance(synonyms, str):
        specs.append(ContextSpec(kind="label", column=synonyms))
    if xref is not None:
        if xref_mapping is None:
            raise ValueError("'xref_mapping' is required when 'xref' is given.")
        specs.append(
            ContextSpec(
                kind="xref", column=xref, xref_mapping=xref_mapping, predicates=xref_predicates
            )
        )
    if context is not None:
        specs.extend([context] if isinstance(context, ContextSpec) else list(context))
    return specs


@overload
def update_ids(
    ids: str,
    mapping_set: Sec2PriMappingSet,
    *,
    at: None = ...,
    suffix: str = ...,
    lookup: dict[str, str] | None = ...,
    ambiguous: set[str] | None = ...,
) -> dict[str, str]:
    """Update IDs from a string."""
    ...


@overload
def update_ids(
    ids: list[str],
    mapping_set: Sec2PriMappingSet,
    *,
    at: None = ...,
    suffix: str = ...,
    lookup: dict[str, str] | None = ...,
    ambiguous: set[str] | None = ...,
) -> dict[str, str]:
    """Update IDs from list[str]."""
    ...


@overload
def update_ids(
    ids: pd.DataFrame,
    mapping_set: Sec2PriMappingSet,
    *,
    at: str | list[str],
    suffix: str = ...,
    lookup: dict[str, str] | None = ...,
    ambiguous: set[str] | None = ...,
    synonyms: str | list[str] | None = ...,
    label_mapping_set: Sec2PriMappingSet | None = ...,
    xref: str | None = ...,
    xref_mapping: XrefMapping | None = ...,
    xref_predicates: set[str] | None = ...,
    report_path: Path | str | None = ...,
    context: ContextSpec | list[ContextSpec] | None = ...,
) -> pd.DataFrame:
    """Update IDs from pd.DataFrame."""
    ...


def update_ids(
    ids: IdsInput,
    mapping_set: Sec2PriMappingSet,
    *,
    at: str | list[str] | None = None,
    suffix: str = "_primary",
    lookup: dict[str, str] | None = None,
    ambiguous: set[str] | None = None,
    synonyms: str | list[str] | None = None,
    label_mapping_set: Sec2PriMappingSet | None = None,
    xref: str | None = None,
    xref_mapping: XrefMapping | None = None,
    xref_predicates: set[str] | None = None,
    report_path: Path | str | None = None,
    context: ContextSpec | list[ContextSpec] | None = None,
) -> dict[str, str] | pd.DataFrame:
    """Resolve secondary identifiers to primary identifiers.

    Parameters
    ----------
    ids:
        One of:

        * **str**: a single identifier, or multiple identifiers joined by
          ``|``, ``,``, ``;``, or whitespace.
        * **list[str]**: a list of identifier strings (each may itself
          contain multiple IDs separated by the delimiters above).
        * **pandas.DataFrame**: a DataFrame; you must also supply *at*.

    mapping_set:
        The :class:`~pysec2pri.parsers.base.Sec2PriMappingSet` to look up
        against (e.g. the result of ``generate_hgnc()``).

    at:
        *DataFrame mode only.* Column name or list of column names that
        contain identifiers.  For each column ``col`` a new column named
        ``col + suffix`` is added to the returned DataFrame.

    suffix:
        Suffix appended to column names in DataFrame mode (default
        ``"_primary"``).

    lookup:
        Pre-built ``{secondary_id: primary_id}`` dictionary.  Pass the
        result of :func:`build_lookup` to avoid rebuilding on repeated
        calls.

    ambiguous:
        Pre-built set of ambiguous IDs (see :func:`build_ambiguous_set`).
        When ``None``, it is computed automatically from *mapping_set*.
        Pass an explicit set (including an empty one) to skip the
        computation.

    synonyms:
        *DataFrame mode only.* Name of a column in the DataFrame that
        contains user-supplied alias strings (delimited by ``|``, ``,``,
        ``;``, or whitespace) to help resolve ambiguous identifiers.
        When provided, :func:`resolve_ambiguous_with_hints` is called for
        every ambiguous cell using that row's alias list.

    label_mapping_set:
        A :class:`~pysec2pri.parsers.base.LabelMappingSet` used to build
        the alias index when *synonyms* is provided.  When ``None`` and
        *synonyms* is set, hint-based resolution falls back to direct
        token/ID matches only and a warning is emitted.

    xref:
        *DataFrame mode only.* Name of a column in the DataFrame that
        contains a per-row cross-reference token (e.g. an Ensembl ID).
        When provided (with *xref_mapping*), :func:`pysec2pri.context.
        resolve_ambiguous_with_xref` is called for every ambiguous cell,
        looking the token up in *xref_mapping* to confirm whether the row
        is using the secondary or the primary identity.

    xref_mapping:
        The :class:`~pysec2pri.context.XrefMapping` crosswalk table to
        resolve *xref* tokens against.  Required when *xref* is given.

    xref_predicates:
        Equivalence predicates accepted from *xref_mapping* records.
        ``None`` (default) accepts any predicate (and unannotated records,
        see :data:`pysec2pri.context.DEFAULT_TRUST_UNANNOTATED`).

    report_path:
        *DataFrame mode only.* When given, every disambiguation attempt
        (from *synonyms*, *xref*, and/or *context*) is written to this path
        as a TSV decision log (see :func:`pysec2pri.context.write_decision_log`).

    context:
        *DataFrame mode only.* One or more explicit
        :class:`~pysec2pri.context.ContextSpec` entries, for cases not
        covered by the *synonyms*/*xref* convenience kwargs (e.g. multiple
        cross-reference columns). Specs are tried in order per ambiguous
        cell; combined with any spec implied by *synonyms*/*xref*.

    Returns
    -------
    dict[str, str]
        When *ids* is a ``str`` or ``list[str]``: a dictionary mapping
        each **unique** input identifier to its resolved primary ID.
        Identifiers not found in the mapping set are returned unchanged.
        **Ambiguous identifiers are mapped to an empty string** and a
        warning is emitted.

    pandas.DataFrame
        When *ids* is a ``DataFrame``: a copy of the DataFrame with one
        new ``<col><suffix>`` column per entry in *at*.  Ambiguous cells
        are set to ``""``; a warning is emitted after all columns are
        processed.  When *synonyms*, *xref*, or *context* is given, a
        companion ``<col><suffix>_id`` column is also added.
    """
    lkp = lookup if lookup is not None else _build_lookup(mapping_set)
    amb = ambiguous if ambiguous is not None else build_ambiguous_set(mapping_set)

    if isinstance(ids, str):
        return _update_str(ids, lkp, amb, kind="ID")
    if isinstance(ids, list):
        return _update_list(ids, lkp, amb, kind="ID")

    # DataFrame mode
    specs = _build_context_specs(synonyms, xref, xref_mapping, xref_predicates, context)
    if specs:
        alias_index: dict[str, list[str]] = {}
        if any(spec.kind in ("label", "id") for spec in specs):
            if label_mapping_set is None:
                logger.warning(
                    "update_ids: a 'label'/'id' context spec was given but no "
                    "'label_mapping_set' was provided; hint-based ambiguity "
                    "resolution will be limited to direct token/ID matches.",
                )
            else:
                alias_index = build_alias_index(label_mapping_set)
        return _update_dataframe_with_context(
            ids,
            at,
            suffix,
            lkp,
            amb,
            kind="ID",
            col_label="ids",
            specs=specs,
            alias_index=alias_index,
            token_to_id=None,  # IDs are their own keys
            report_path=report_path,
        )
    return _update_dataframe(ids, at, suffix, lkp, amb, kind="ID", col_label="ids")


def build_label_lookup(mapping_set: Sec2PriMappingSet) -> dict[str, str]:
    """Return a ``{secondary_label: primary_label}`` dictionary.

    Useful when you want to apply the look-up yourself or cache it for
    repeated calls.

    Args:
        mapping_set: A :class:`~pysec2pri.parsers.base.LabelMappingSet`
            (e.g. the result of ``generate_hgnc_labels()``).

    Returns:
        Dictionary mapping every previous/alias label to its current label.
    """
    lookup: dict[str, str] = {}
    for m in mapping_set.mappings or []:
        sec = str(getattr(m, "subject_label", None) or "")
        pri = str(getattr(m, "object_label", None) or "")
        if sec:
            lookup[sec] = pri
    return lookup


@overload
def update_labels(
    labels: str,
    mapping_set: Sec2PriMappingSet,
    *,
    at: None = ...,
    suffix: str = ...,
    lookup: dict[str, str] | None = ...,
    ambiguous: set[str] | None = ...,
    synonyms: str | list[str] | None = ...,
) -> dict[str, str]:
    """Input as str."""
    ...


@overload
def update_labels(
    labels: list[str],
    mapping_set: Sec2PriMappingSet,
    *,
    at: None = ...,
    suffix: str = ...,
    lookup: dict[str, str] | None = ...,
    ambiguous: set[str] | None = ...,
    synonyms: str | list[str] | None = ...,
) -> dict[str, str]:
    """Input as list[str]."""
    ...


@overload
def update_labels(
    labels: pd.DataFrame,
    mapping_set: Sec2PriMappingSet,
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
) -> pd.DataFrame:
    """Input as pd.DataFrame."""
    ...


def update_labels(
    labels: IdsInput,
    mapping_set: Sec2PriMappingSet,
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
    """Resolve previous/alias gene labels to current labels.

    Same as :func:`update_ids` but resolves via the
    ``subject_label`` to ``object_label`` mapping rather than IDs.

    Parameters
    ----------
    labels:
        One of:

        * **str**: a single label, or multiple labels joined by
          ``|``, ``,``, ``;``, or whitespace.
        * **list[str]**: a list of label strings.
        * **pandas.DataFrame**: a DataFrame; you must also supply *at*.

    mapping_set:
        A :class:`~pysec2pri.parsers.base.LabelMappingSet`
        (e.g. the result of ``generate_hgnc_labels()``).

    at:
        *DataFrame mode only.* Column name or list of column names that
        contain labels.  For each column ``col`` a new column named
        ``col + suffix`` is added to the returned DataFrame.

    suffix:
        Suffix appended to column names in DataFrame mode (default
        ``"_current"``).

    lookup:
        Pre-built ``{previous_label: current_label}`` dictionary.
        Pass the result of :func:`build_label_lookup` to avoid rebuilding
        on repeated calls.

    ambiguous:
        Pre-built set of ambiguous labels (see
        :func:`build_ambiguous_labels_set`).  When ``None``, it is
        computed automatically from *mapping_set*.

    synonyms:
        *DataFrame mode only.* Name of a column in the DataFrame that
        contains user-supplied alias strings (delimited by ``|``, ``,``,
        ``;``, or whitespace) to help resolve ambiguous labels.
        When provided, :func:`resolve_ambiguous_with_hints` is called for
        every ambiguous cell using that row's alias list.  The alias index
        is built from *mapping_set* itself (non-IAO entries).

    xref:
        *DataFrame mode only.* Name of a column in the DataFrame that
        contains a per-row cross-reference token (e.g. an Ensembl ID).
        When provided (with *xref_mapping*), :func:`pysec2pri.context.
        resolve_ambiguous_with_xref` is called for every ambiguous cell,
        looking the token up in *xref_mapping* to confirm whether the row
        is using the secondary or the primary identity.

    xref_mapping:
        The :class:`~pysec2pri.context.XrefMapping` crosswalk table to
        resolve *xref* tokens against.  Required when *xref* is given.

    xref_predicates:
        Equivalence predicates accepted from *xref_mapping* records.
        ``None`` (default) accepts any predicate (and unannotated records,
        see :data:`pysec2pri.context.DEFAULT_TRUST_UNANNOTATED`).

    report_path:
        *DataFrame mode only.* When given, every disambiguation attempt
        (from *synonyms*, *xref*, and/or *context*) is written to this path
        as a TSV decision log (see :func:`pysec2pri.context.write_decision_log`).

    context:
        *DataFrame mode only.* One or more explicit
        :class:`~pysec2pri.context.ContextSpec` entries, for cases not
        covered by the *synonyms*/*xref* convenience kwargs (e.g. multiple
        cross-reference columns). Specs are tried in order per ambiguous
        cell; combined with any spec implied by *synonyms*/*xref*.

    Returns
    -------
    dict[str, str]
        When *labels* is a ``str`` or ``list[str]``: a dictionary mapping
        each unique input label to its resolved current label.  Symbols
        not found in the mapping set are returned unchanged.  **Ambiguous
        labels are mapped to an empty string** and a warning is emitted.

    pandas.DataFrame
        When *labels* is a ``DataFrame``: a copy of the DataFrame with one
        new ``<col><suffix>`` column per entry in *at*.  Ambiguous cells
        are set to ``""``; a warning is emitted after all columns are
        processed.  When *synonyms*, *xref*, or *context* is given, a
        companion ``<col><suffix>_id`` column is also added.
    """
    lkp = lookup if lookup is not None else build_label_lookup(mapping_set)
    amb = ambiguous if ambiguous is not None else build_ambiguous_labels_set(mapping_set)

    if isinstance(labels, str):
        return _update_str(labels, lkp, amb, kind="label")
    if isinstance(labels, list):
        return _update_list(labels, lkp, amb, kind="label")

    # DataFrame mode
    specs = _build_context_specs(synonyms, xref, xref_mapping, xref_predicates, context)
    if specs:
        token_to_id = build_primary_token_to_id(mapping_set)
        alias_index: dict[str, list[str]] = (
            build_alias_index(mapping_set)
            if any(spec.kind in ("label", "id") for spec in specs)
            else {}
        )
        return _update_dataframe_with_context(
            labels,
            at,
            suffix,
            lkp,
            amb,
            kind="label",
            col_label="labels",
            specs=specs,
            alias_index=alias_index,
            token_to_id=token_to_id,
            report_path=report_path,
        )
    return _update_dataframe(labels, at, suffix, lkp, amb, kind="label", col_label="labels")


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
