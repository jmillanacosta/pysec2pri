"""Export functions for writing mapping sets to various file formats."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pysec2pri.parsers.base import Sec2PriMappingSet

__all__ = [
    "write_name2synonym",
    "write_pri_ids",
    "write_sec2pri",
    "write_sssom",
    "write_symbol2prev",
]


def write_sssom(
    mapping_set: Sec2PriMappingSet,
    output_path: Path | str,
) -> Path:
    """Write a MappingSet to an SSSOM TSV file."""
    import codecs
    import re
    from typing import Any, cast

    import curies
    from sssom.parsers import (  # type: ignore[attr-defined]
        to_mapping_set_dataframe,
    )
    from sssom.sssom_document import MappingSetDocument
    from sssom.writers import write_table

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build converter from curie_map
    raw_curie_map: Any = mapping_set.curie_map or {}
    # Convert curie_map to format expected by curies
    # Values may be strings or Prefix objects (sssom_schema converts them)
    records: list[curies.Record] = []
    if isinstance(raw_curie_map, dict):
        for k, v in raw_curie_map.items():
            if isinstance(v, str):
                uri_prefix = v
            elif hasattr(v, "prefix_url"):
                uri_prefix = cast(str, v.prefix_url)
            else:
                continue
            records.append(curies.Record(prefix=k, uri_prefix=uri_prefix))
    converter = curies.Converter(records=records)

    # Wrap in MappingSetDocument for sssom library compatibility
    doc = MappingSetDocument(mapping_set=mapping_set, converter=converter)
    msdf = to_mapping_set_dataframe(doc)

    with output_path.open("w", encoding="utf-8") as f:
        write_table(msdf, f)

    # Fix escaped unicode in YAML header (sssom issue: )
    content = output_path.read_text(encoding="utf-8")
    # Replace \xNN escapes with actual unicode characters
    content = re.sub(
        r"\\x([0-9a-fA-F]{2})",
        lambda m: codecs.decode(bytes([int(m.group(1), 16)]), "latin-1"),
        content,
    )
    output_path.write_text(content, encoding="utf-8")

    return output_path


def write_pri_ids(
    mapping_set: Sec2PriMappingSet,
    output_path: Path | str,
) -> Path:
    """Write unique primary IDs (object_id) to a text file, one per line."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Collect unique primary IDs
    pri_ids: set[str] = set()
    for m in mapping_set.mappings or []:
        obj_id = getattr(m, "object_id", None)
        if obj_id:
            pri_ids.add(str(obj_id))

    with output_path.open("w", encoding="utf-8") as f:
        for pri_id in sorted(pri_ids):
            f.write(f"{pri_id}\n")

    return output_path


def write_sec2pri(
    mapping_set: Sec2PriMappingSet,
    output_path: Path | str,
) -> Path:
    """Write secondary to primary ID mappings to a TSV file."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    columns = ["subject_id", "object_id", "predicate_id", "mapping_cardinality"]

    with output_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for m in mapping_set.mappings or []:
            values = [
                str(getattr(m, "subject_id", "") or ""),
                str(getattr(m, "object_id", "") or ""),
                str(getattr(m, "predicate_id", "") or ""),
                str(getattr(m, "mapping_cardinality", "") or ""),
            ]
            f.write("\t".join(values) + "\n")

    return output_path


def write_name2synonym(
    mapping_set: Sec2PriMappingSet,
    output_path: Path | str,
) -> Path:
    """Write name to synonym mappings to a TSV file."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    columns = ["subject_id", "subject_label", "object_label"]

    with output_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for m in mapping_set.mappings or []:
            subject_label = getattr(m, "subject_label", None)
            object_label = getattr(m, "object_label", None)
            if subject_label or object_label:
                values = [
                    str(getattr(m, "subject_id", "") or ""),
                    str(subject_label or ""),
                    str(object_label or ""),
                ]
                f.write("\t".join(values) + "\n")

    return output_path


def write_symbol2prev(
    mapping_set: Sec2PriMappingSet,
    output_path: Path | str,
) -> Path:
    """Write symbol to previous symbol mappings to a TSV file."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    columns = ["subject_id", "subject_label", "object_label", "mapping_cardinality"]

    with output_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for m in mapping_set.mappings or []:
            subject_label = getattr(m, "subject_label", None)
            object_label = getattr(m, "object_label", None)
            if subject_label or object_label:
                values = [
                    str(getattr(m, "subject_id", "") or ""),
                    str(subject_label or ""),
                    str(object_label or ""),
                    str(getattr(m, "mapping_cardinality", "") or ""),
                ]
                f.write("\t".join(values) + "\n")

    return output_path
