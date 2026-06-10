"""Export functions for writing mapping sets to various file formats."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING

from sssom import MappingSetDataFrame

if TYPE_CHECKING:
    from pysec2pri.parsers.base import Sec2PriMappingSet

__all__ = [
    "WRITERS",
    "write_json",
    "write_label2prev",
    "write_name2synonym",
    "write_output",
    "write_owl",
    "write_pri_ids",
    "write_rdf",
    "write_sec2pri",
    "write_secondary",
    "write_sssom",
]


def write_sssom(
    mapping_set: Sec2PriMappingSet,
    output_path: Path | str,
) -> Path:
    """Write a mapping set to an SSSOM TSV file.

    Args:
        mapping_set: The mapping set to write.
        output_path: Destination ``.sssom.tsv`` file path.

    Returns:
        Path to the written file.
    """
    import codecs
    import re
    from typing import cast

    import curies
    from sssom.parsers import to_mapping_set_dataframe  # type: ignore[attr-defined]
    from sssom.sssom_document import MappingSetDocument
    from sssom.writers import write_table

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build a curies.Converter from the curie_map stored on the mapping set.
    # Using Converter(records=...) preserves prefix casing exactly as declared
    # and handles pydantic-wrapped Prefix objects (which expose .prefix_url).
    raw_curie_map: object = mapping_set.curie_map or {}
    records: list[curies.Record] = []
    if isinstance(raw_curie_map, dict):
        for k, v in raw_curie_map.items():
            if isinstance(v, str):
                uri_prefix: str = v
            elif hasattr(v, "prefix_url"):
                uri_prefix = cast(str, v.prefix_url)
            else:
                continue
            records.append(curies.Record(prefix=k, uri_prefix=uri_prefix))
    converter = curies.Converter(records=records)
    doc = MappingSetDocument(mapping_set=mapping_set, converter=converter)
    msdf = to_mapping_set_dataframe(doc)

    with output_path.open("w", encoding="utf-8") as f:
        write_table(msdf, f)

    # Fix escaped unicode in YAML header (sssom issue)
    content = output_path.read_text(encoding="utf-8")
    content = re.sub(
        r"\\x([0-9a-fA-F]{2})",
        lambda m: codecs.decode(bytes([int(m.group(1), 16)]), "latin-1"),
        content,
    )
    output_path.write_text(content, encoding="utf-8")

    return output_path


def _to_msdf_via_sssom_parser(mapping_set: Sec2PriMappingSet) -> MappingSetDataFrame | None:
    """Write to a temporary SSSOM TSV then parse back with sssom's own parser.

    Args:
        mapping_set: The mapping set to convert.

    Returns:
        A fully-validated ``MappingSetDataFrame`` ready for RDF/JSON/OWL serialisation.
    """
    import tempfile

    from sssom.parsers import parse_sssom_table

    with tempfile.NamedTemporaryFile(
        suffix=".sssom.tsv", mode="w", encoding="utf-8", delete=False
    ) as tmp:
        tmp_path = Path(tmp.name)

    try:
        write_sssom(mapping_set, tmp_path)
        return parse_sssom_table(str(tmp_path))
    finally:
        tmp_path.unlink(missing_ok=True)


def write_rdf(
    mapping_set: Sec2PriMappingSet,
    output_path: Path | str,
    serialisation: str = "turtle",
) -> Path:
    """Write a mapping set to an RDF file.

    Args:
        mapping_set: The mapping set to write.
        output_path: Destination file path (e.g. ``mappings.ttl``).
        serialisation: RDFLib serialisation format.

    Returns:
        Path to the written file.
    """
    from sssom.writers import write_rdf as _sssom_write_rdf

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    msdf = _to_msdf_via_sssom_parser(mapping_set)
    if msdf is None:
        raise ValueError("Failed to parse mapping set for RDF serialisation.")
    with output_path.open("w", encoding="utf-8") as f:
        _sssom_write_rdf(msdf, f, serialisation=serialisation)
    return output_path


def write_json(
    mapping_set: Sec2PriMappingSet,
    output_path: Path | str,
) -> Path:
    """Write a mapping set to an SSSOM JSON file.

    Args:
        mapping_set: The mapping set to write.
        output_path: Destination file path (e.g. ``mappings.json``).

    Returns:
        Path to the written file.
    """
    from sssom.writers import write_json as _sssom_write_json

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    msdf = _to_msdf_via_sssom_parser(mapping_set)
    if msdf is None:
        raise ValueError("Failed to parse mapping set for JSON serialisation.")
    with output_path.open("w", encoding="utf-8") as f:
        _sssom_write_json(msdf, f)
    return output_path


def write_owl(
    mapping_set: Sec2PriMappingSet,
    output_path: Path | str,
    serialisation: str = "turtle",
) -> Path:
    """Write a mapping set to an OWL/RDF file (default: Turtle).

    Args:
        mapping_set: The mapping set to write.
        output_path: Destination file path (e.g. ``mappings_owl.ttl``).
        serialisation: RDFLib serialisation format.

    Returns:
        Path to the written file.
    """
    from sssom.writers import write_owl as _sssom_write_owl

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    msdf = _to_msdf_via_sssom_parser(mapping_set)
    if msdf is None:
        raise ValueError("Failed to parse mapping set for OWL serialisation.")
    with output_path.open("w", encoding="utf-8") as f:
        _sssom_write_owl(msdf, f, serialisation=serialisation)
    return output_path


def write_pri_ids(
    mapping_set: Sec2PriMappingSet,
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

    # Use the authoritative set when available (never appears in other outputs)
    all_ids: set[str] = getattr(mapping_set, "_primary_ids", set()) or set()

    if not all_ids:
        # Fall back to extracting from mappings
        for m in mapping_set.mappings or []:  # type: ignore[has-type]
            obj_id = getattr(m, "object_id", None)
            if obj_id:
                all_ids.add(str(obj_id))

    with output_path.open("w", encoding="utf-8") as f:
        for pri_id in sorted(all_ids):
            f.write(f"{pri_id}\n")

    return output_path


def write_pri_labels(
    mapping_set: Sec2PriMappingSet,
    output_path: Path | str,
) -> Path:
    """Write unique primary labels to a text file, one per line.

    Args:
        mapping_set: The mapping set to read primary IDs from.
        output_path: Destination file path.

    Returns:
        Path to the written file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Use the authoritative set when available (never appears in other outputs)
    all_sym: set[str] = getattr(mapping_set, "_primary_labels", set()) or set()

    if not all_sym:
        # Fall back to extracting from mappings
        for m in mapping_set.mappings or []:  # type: ignore[has-type]
            obj_label = getattr(m, "object_label", None)
            if obj_label:
                all_sym.add(str(obj_label))

    with output_path.open("w", encoding="utf-8") as f:
        for pri_sym in sorted(all_sym):
            f.write(f"{pri_sym}\n")

    return output_path


def write_sec2pri(
    mapping_set: Sec2PriMappingSet,
    output_path: Path | str,
) -> Path:
    """Write secondary to primary ID mappings to a TSV file.

    Columns: ``subject_id``, ``object_id``, ``predicate_id``, ``mapping_cardinality``.

    Args:
        mapping_set: The mapping set to write.
        output_path: Destination file path (e.g. ``sec2pri.tsv``).

    Returns:
        Path to the written file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Use explicit primary/secondary column names for clarity
    # Columns: primary_id (object_id), secondary_id (subject_id), predicate_id, mapping_cardinality
    columns = ["primary_id", "secondary_id", "predicate_id", "mapping_cardinality"]

    with output_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for m in mapping_set.mappings or []:  # type: ignore[has-type]
            values = [
                str(getattr(m, "object_id", "") or ""),
                str(getattr(m, "subject_id", "") or ""),
                str(getattr(m, "predicate_id", "") or ""),
                str(getattr(m, "mapping_cardinality", "") or ""),
            ]
            f.write("\t".join(values) + "\n")

    return output_path


def write_name2synonym(
    mapping_set: Sec2PriMappingSet,
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
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Columns: primary_id (object), name (object_label/primary), synonym (subject_label)
    columns = ["primary_id", "name", "synonym"]

    with output_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for m in mapping_set.mappings or []:  # type: ignore[has-type]
            # Only emit synonym rows; skip deprecation (IAO:0100001) rows
            if getattr(m, "predicate_id", None) != "oboInOwl:hasExactSynonym":
                continue
            subject_label = getattr(m, "subject_label", None)
            object_label = getattr(m, "object_label", None)
            if subject_label or object_label:
                # object_label is the primary (name), subject_label is the synonym
                values = [
                    str(getattr(m, "object_id", "") or ""),
                    str(object_label or ""),
                    str(subject_label or ""),
                ]
                f.write("\t".join(values) + "\n")

    return output_path


def write_label2prev(
    mapping_set: Sec2PriMappingSet,
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
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Columns: primary_id (object), primary_label (object_label), previous_label (subject_label)
    columns = ["primary_id", "primary_label", "previous_label", "mapping_cardinality"]

    with output_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for m in mapping_set.mappings or []:  # type: ignore[has-type]
            # Only emit deprecation rows; skip synonym (hasExactSynonym) rows
            if getattr(m, "predicate_id", None) == "oboInOwl:hasExactSynonym":
                continue
            subject_label = getattr(m, "subject_label", None)
            object_label = getattr(m, "object_label", None)
            if subject_label or object_label:
                # object_label is the primary label, subject_label is the previous label
                values = [
                    str(getattr(m, "object_id", "") or ""),
                    str(object_label or ""),
                    str(subject_label or ""),
                    str(getattr(m, "mapping_cardinality", "") or ""),
                ]
                f.write("\t".join(values) + "\n")

    return output_path


def write_secondary(
    mapping_set: Sec2PriMappingSet,
    output_path: Path | str,
) -> Path:
    """Write unique secondary IDs (subject_id) to a text file, one per line.

    Args:
        mapping_set: The mapping set to read secondary IDs from.
        output_path: Destination file path.

    Returns:
        Path to the written file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sec_ids: set[str] = set()
    for m in mapping_set.mappings or []:  # type: ignore[has-type]
        subj_id = getattr(m, "subject_id", None)
        if subj_id:
            sec_ids.add(str(subj_id))

    with output_path.open("w", encoding="utf-8") as f:
        for sec_id in sorted(sec_ids):
            f.write(f"{sec_id}\n")

    return output_path


# Registry mapping config output names to writer functions.
WRITERS: dict[str, Callable[..., Path]] = {
    "sssom": write_sssom,
    "sec2pri": write_sec2pri,
    "secID2priID": write_sec2pri,
    "pri_ids": write_pri_ids,
    "priIDs": write_pri_ids,
    "secIDs": write_secondary,
    "name2synonym": write_name2synonym,
    "label2prev": write_label2prev,
    "rdf": write_rdf,
    "json": write_json,
    "owl": write_owl,
}


def write_output(
    mapping_set: Sec2PriMappingSet,
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
