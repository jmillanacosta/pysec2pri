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
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING, Union, overload

if TYPE_CHECKING:
    import pandas as pd

from pysec2pri.parsers.base import Sec2PriMappingSet

# Separator pattern: pipe, comma, semicolon, or whitespace
_SEP = re.compile(r"[|,;\s]+")

# Type alias for the flexible input
IdsInput = Union[str, list[str], "pd.DataFrame"]


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


def _split_ids(raw: str) -> list[str]:
    """Split *raw* on common delimiters; strips empty tokens."""
    return [tok for tok in _SEP.split(raw.strip()) if tok]


def _resolve_tokens(tokens: list[str], lookup: dict[str, str]) -> list[str]:
    """Map each token to its primary ID, falling back to the token itself."""
    return [lookup.get(tok, tok) for tok in tokens]


def _resolve_string(raw: str, lookup: dict[str, str]) -> str:
    """Resolve all IDs inside *raw* and rejoin with the original separator.

    The original separator character (first one found) is preserved.  If
    multiple different separators are used inside *raw* the first one wins.
    """
    sep_match = re.search(r"[|,;\s]", raw)
    sep = sep_match.group(0) if sep_match else ""
    tokens = _split_ids(raw)
    resolved = _resolve_tokens(tokens, lookup)
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


@overload
def update_ids(
    ids: str,
    mapping_set: Sec2PriMappingSet,
    *,
    at: None = ...,
    suffix: str = ...,
    lookup: dict[str, str] | None = ...,
) -> dict[str, str]:
    """Overload: string input to ``{id: primary_id}`` dict."""
    ...


@overload
def update_ids(
    ids: list[str],
    mapping_set: Sec2PriMappingSet,
    *,
    at: None = ...,
    suffix: str = ...,
    lookup: dict[str, str] | None = ...,
) -> dict[str, str]:
    """Overload: list of strings input to ``{id: primary_id}`` dict."""
    ...


@overload
def update_ids(
    ids: pd.DataFrame,
    mapping_set: Sec2PriMappingSet,
    *,
    at: str | list[str],
    suffix: str = ...,
    lookup: dict[str, str] | None = ...,
) -> pd.DataFrame:
    """Overload: DataFrame input to DataFrame with added primary-ID columns."""
    ...


def update_ids(
    ids: IdsInput,
    mapping_set: Sec2PriMappingSet,
    *,
    at: str | list[str] | None = None,
    suffix: str = "_primary",
    lookup: dict[str, str] | None = None,
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
        * **pandas.DataFrame**: a DataFrame; you must also supply at.

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

    Returns
    -------
    dict[str, str]
        When *ids* is a ``str`` or ``list[str]``: a dictionary mapping
        each **unique** input identifier to its resolved primary ID.
        Identifiers not found in the mapping set are returned unchanged.

    pandas.DataFrame
        When *ids* is a ``DataFrame``: a copy of the DataFrame with one
        new ``<col><suffix>`` column per entry in *at*.

    Examples
    --------
    Setup::

        ms = generate_hgnc()

    Single string::

        update_ids("HGNC:1234", ms)
        # {'HGNC:1234': 'HGNC:9999'}

    Pipe-separated string::

        update_ids("HGNC:1234|HGNC:5678", ms)
        # {'HGNC:1234': 'HGNC:9999', 'HGNC:5678': 'HGNC:5678'}

    List::

        update_ids(["HGNC:1234", "HGNC:5678", "HGNC:1234"], ms)
        # {'HGNC:1234': 'HGNC:9999', 'HGNC:5678': 'HGNC:5678'}

    DataFrame::

        import pandas as pd

        df = pd.DataFrame({"gene": ["HGNC:1234", "HGNC:5678"]})
        update_ids(df, ms, at="gene")
        #        gene  gene_primary
        # 0  HGNC:1234  HGNC:9999
        # 1  HGNC:5678  HGNC:5678
    """
    lkp = lookup if lookup is not None else _build_lookup(mapping_set)

    # a str
    if isinstance(ids, str):
        tokens = _split_ids(ids)
        unique = dict.fromkeys(tokens)  # preserves insertion order, deduplicates
        return {tok: lkp.get(tok, tok) for tok in unique}

    # a list[str]
    if isinstance(ids, list):
        # Flatten: each element may itself be a multi-ID string
        unique_l: dict[str, None] = {}
        for item in ids:
            for tok in _split_ids(item):
                unique_l[tok] = None
        return {tok: lkp.get(tok, tok) for tok in unique_l}

    # a pd.DataFrame
    import pandas as pd

    if not isinstance(ids, pd.DataFrame):
        raise TypeError(
            f"'ids' must be a str, list[str], or pandas.DataFrame, got {type(ids).__name__!r}."
        )

    if at is None:
        raise ValueError(
            "When 'ids' is a DataFrame you must specify 'at' (column name or list of names)."
        )

    columns: list[str] = [at] if isinstance(at, str) else list(at)
    missing = [c for c in columns if c not in ids.columns]
    if missing:
        raise KeyError(f"Column(s) not found in DataFrame: {missing}")

    result = ids.copy()
    for col in columns:
        new_col = col + suffix
        result[new_col] = (
            result[col].astype(str).map(lambda cell, _lkp=lkp: _resolve_string(cell, _lkp))
        )

    return result


__all__ = [
    "build_lookup",
    "update_ids",
]
