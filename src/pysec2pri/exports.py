"""Export functions for writing mapping sets to various file formats."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING

from mapkgsutils.exports import write_json, write_owl, write_rdf, write_sssom

if TYPE_CHECKING:
    import pandas as pd

    from pysec2pri.parsers.base import BaseMappingSet, IdMappingSet, LabelMappingSet

__all__ = [
    "WRITERS",
    "write_json",
    "write_label2prev",
    "write_label_sec2pri",
    "write_name2synonym",
    "write_output",
    "write_owl",
    "write_pri_ids",
    "write_pri_labels",
    "write_rdf",
    "write_sec2pri",
    "write_secondary",
    "write_sssom",
]


def write_pri_ids(
    mapping_set: IdMappingSet,
    output_path: Path | str,
) -> Path:
    """Write unique primary IDs to a text file, one per line.

    Args:
        mapping_set: The mapping set to read primary IDs from.
        output_path: Destination file path.

    Returns:
        Path to the written file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ids = mapping_set.to_pri_ids()
    output_path.write_text("\n".join(ids) + "\n", encoding="utf-8")
    return output_path


def write_pri_labels(
    mapping_set: LabelMappingSet,
    output_path: Path | str,
) -> Path:
    """Write unique ``primary_id``/``label`` pairs to a two-column TSV.

    Args:
        mapping_set: The mapping set to read primary labels from.
        output_path: Destination file path.

    Returns:
        Path to the written file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pairs = mapping_set.to_pri_labels()
    body = "\n".join(f"{pri_id}\t{label}" for pri_id, label in pairs)
    output_path.write_text("id\tlabel\n" + body + "\n", encoding="utf-8")
    return output_path


def _write_tsv(df: pd.DataFrame, output_path: Path | str, *, header: bool = True) -> Path:
    """Write *df* to *output_path* as a TSV, creating parent directories as needed.

    Shared by every ``write_*``/``to_*`` pair below, so a mapping set's tabular
    formats are built once (as a DataFrame) and written the same way whether
    the caller wants the in-memory object or the file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False, header=header)
    return output_path


def _sec2pri_frame(mapping_set: BaseMappingSet) -> pd.DataFrame:
    """Build the secondary-to-primary ID table for ``write_sec2pri``/``to_sec2pri``.

    Columns: ``primary_id`` (object_id), ``secondary_id`` (subject_id),
    ``predicate_id``, ``mapping_cardinality``.
    """
    import pandas as pd

    rows = [
        {
            "primary_id": str(getattr(m, "object_id", "") or ""),
            "secondary_id": str(getattr(m, "subject_id", "") or ""),
            "predicate_id": str(getattr(m, "predicate_id", "") or ""),
            "mapping_cardinality": str(getattr(m, "mapping_cardinality", "") or ""),
        }
        for m in (mapping_set.mappings or [])
    ]
    return pd.DataFrame(
        rows, columns=["primary_id", "secondary_id", "predicate_id", "mapping_cardinality"]
    )


def _label_sec2pri_frame(mapping_set: BaseMappingSet) -> pd.DataFrame:
    """Build the full label sec2pri table for ``write_label_sec2pri``/``to_label_sec2pri``.

    Every mapping row is included (both deprecation and synonym predicates).
    Columns: ``secondary_id``, ``secondary_label``, ``primary_id``,
    ``primary_label``, ``predicate_id``, ``mapping_cardinality``.
    """
    import pandas as pd

    rows = [
        {
            "secondary_id": str(getattr(m, "subject_id", "") or ""),
            "secondary_label": str(getattr(m, "subject_label", "") or ""),
            "primary_id": str(getattr(m, "object_id", "") or ""),
            "primary_label": str(getattr(m, "object_label", "") or ""),
            "predicate_id": str(getattr(m, "predicate_id", "") or ""),
            "mapping_cardinality": str(getattr(m, "mapping_cardinality", "") or ""),
        }
        for m in (mapping_set.mappings or [])
    ]
    return pd.DataFrame(
        rows,
        columns=[
            "secondary_id",
            "secondary_label",
            "primary_id",
            "primary_label",
            "predicate_id",
            "mapping_cardinality",
        ],
    )


def _name2synonym_frame(mapping_set: BaseMappingSet) -> pd.DataFrame:
    """Build the name-to-synonym table for ``write_name2synonym``/``to_name2synonym``.

    Only ``oboInOwl:hasExactSynonym`` rows are included; deprecation rows
    (``IAO:0100001``) belong in ``label_sec2pri``/``label2prev``, not here.
    Columns: ``primary_id``, ``name``, ``synonym``.
    """
    import pandas as pd

    rows = [
        {
            "primary_id": str(getattr(m, "object_id", "") or ""),
            "name": str(getattr(m, "object_label", "") or ""),
            "synonym": str(getattr(m, "subject_label", "") or ""),
        }
        for m in (mapping_set.mappings or [])
        if getattr(m, "predicate_id", None) == "oboInOwl:hasExactSynonym"
        and (getattr(m, "subject_label", None) or getattr(m, "object_label", None))
    ]
    return pd.DataFrame(rows, columns=["primary_id", "name", "synonym"])


def _label2prev_frame(mapping_set: BaseMappingSet) -> pd.DataFrame:
    """Build the label-to-previous-label table for ``write_label2prev``/``to_label2prev``.

    Only ``IAO:0100001`` (``"term replaced by"``) rows are included; synonym
    rows (``oboInOwl:hasExactSynonym``) belong in ``name2synonym``, not here.
    Columns: ``primary_id``, ``primary_label``, ``previous_label``,
    ``mapping_cardinality``.
    """
    import pandas as pd

    rows = [
        {
            "primary_id": str(getattr(m, "object_id", "") or ""),
            "primary_label": str(getattr(m, "object_label", "") or ""),
            "previous_label": str(getattr(m, "subject_label", "") or ""),
            "mapping_cardinality": str(getattr(m, "mapping_cardinality", "") or ""),
        }
        for m in (mapping_set.mappings or [])
        if getattr(m, "predicate_id", None) != "oboInOwl:hasExactSynonym"
        and (getattr(m, "subject_label", None) or getattr(m, "object_label", None))
    ]
    return pd.DataFrame(
        rows, columns=["primary_id", "primary_label", "previous_label", "mapping_cardinality"]
    )


def _secondary_frame(mapping_set: BaseMappingSet) -> pd.DataFrame:
    """Build the unique-secondary-ID table for ``write_secondary``/``to_secondary``."""
    import pandas as pd

    ids = sorted(
        {str(getattr(m, "subject_id", None) or "") for m in (mapping_set.mappings or [])} - {""}
    )
    return pd.DataFrame({"secondary_id": ids})


def write_sec2pri(
    mapping_set: BaseMappingSet,
    output_path: Path | str,
) -> Path:
    """Write secondary to primary ID mappings to a TSV file.

    Columns: ``primary_id`` (object_id), ``secondary_id`` (subject_id),
    ``predicate_id``, ``mapping_cardinality``.

    Args:
        mapping_set: The mapping set to write.
        output_path: Destination file path (e.g. ``sec2pri.tsv``).

    Returns:
        Path to the written file.
    """
    return _write_tsv(_sec2pri_frame(mapping_set), output_path)


def write_label_sec2pri(
    mapping_set: BaseMappingSet,
    output_path: Path | str,
) -> Path:
    """Write the full previous/alias-label to current-label table to a TSV file.

    Every mapping row is written (both deprecation and synonym predicates).
    Columns: ``secondary_id``, ``secondary_label``, ``primary_id``,
    ``primary_label``, ``predicate_id``, ``mapping_cardinality``.

    Args:
        mapping_set: The mapping set to write.
        output_path: Destination file path (e.g. ``label_sec2pri.tsv``).

    Returns:
        Path to the written file.
    """
    return _write_tsv(_label_sec2pri_frame(mapping_set), output_path)


def write_name2synonym(
    mapping_set: BaseMappingSet,
    output_path: Path | str,
) -> Path:
    """Write name to synonym mappings to a TSV file.

    Only ``oboInOwl:hasExactSynonym`` rows are written; deprecation rows
    (``IAO:0100001``) are excluded because they belong in the
    ``label2prev`` output.  Columns: ``primary_id``, ``name``, ``synonym``.

    Args:
        mapping_set: The mapping set to write.
        output_path: Destination file path (e.g. ``name2synonym.tsv``).

    Returns:
        Path to the written file.
    """
    return _write_tsv(_name2synonym_frame(mapping_set), output_path)


def write_label2prev(
    mapping_set: BaseMappingSet,
    output_path: Path | str,
) -> Path:
    """Write label to previous (deprecated) label mappings to a TSV file.

    Only ``IAO:0100001`` (``"term replaced by"``) rows are written; synonym
    rows (``oboInOwl:hasExactSynonym``) are excluded because they belong in
    the ``name2synonym`` output.  Columns: ``primary_id``, ``primary_label``,
    ``previous_label``, ``mapping_cardinality``.

    Args:
        mapping_set: The mapping set to write.
        output_path: Destination file path (e.g. ``label2prev.tsv``).

    Returns:
        Path to the written file.
    """
    return _write_tsv(_label2prev_frame(mapping_set), output_path)


def write_secondary(
    mapping_set: BaseMappingSet,
    output_path: Path | str,
) -> Path:
    """Write unique secondary IDs (subject_id) to a text file, one per line.

    Args:
        mapping_set: The mapping set to read secondary IDs from.
        output_path: Destination file path.

    Returns:
        Path to the written file.
    """
    return _write_tsv(_secondary_frame(mapping_set), output_path, header=False)


# Registry mapping config output names to writer functions.
WRITERS: dict[str, Callable[..., Path]] = {
    "sssom": write_sssom,
    "sec2pri": write_sec2pri,
    "secID2priID": write_sec2pri,
    "pri_ids": write_pri_ids,
    "priIDs": write_pri_ids,
    "pri_labels": write_pri_labels,
    "secIDs": write_secondary,
    "name2synonym": write_name2synonym,
    "label_sec2pri": write_label_sec2pri,
    "label2prev": write_label2prev,
    "rdf": write_rdf,
    "json": write_json,
    "owl": write_owl,
}


def write_output(
    mapping_set: BaseMappingSet,
    output_format: str,
    output_path: Path | str,
) -> Path:
    """Write a mapping set in any registered output format.

    Args:
        mapping_set: The mapping set to write.
        output_format: Format name (must be a key in WRITERS).
        output_path: Path to write to.

    Returns:
        Path to the written file.

    Raises:
        ValueError: If output_format is not recognized.
    """
    writer = WRITERS.get(output_format)
    if writer is None:
        msg = f"Unknown output format: {output_format!r}. Available: {sorted(WRITERS)}"
        raise ValueError(msg)
    return writer(mapping_set, output_path)
