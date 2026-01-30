"""Export functions for writing mapping sets to various file formats.

This module provides functions for writing MappingSet objects to:
- SSSOM TSV format (standard semantic mapping format)
- Legacy TSV formats matching the original R script outputs
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from pysec2pri.constants import (
    DatasourceConfig,
    OutputType,
    get_datasource_config,
    get_output_filename,
)
from pysec2pri.models import BaseMapping

if TYPE_CHECKING:
    from pysec2pri.models import MappingSet

__all__ = [
    "write_name2synonym",
    # Multi-output helper
    "write_outputs",
    "write_sec2pri",
    # SSSOM output
    "write_sssom",
    # Legacy output formats
    "write_subject_ids",
    "write_symbol2prev",
]


# =============================================================================
# SSSOM Output
# =============================================================================


def write_sssom(
    mapping_set: MappingSet,
    output_path: Path | str,
    mapping_date: str | None = None,
) -> None:
    """Write a MappingSet to an SSSOM TSV file.

    Uses the official sssom-py library for correct SSSOM format output.

    Args:
        mapping_set: The pysec2pri MappingSet to write.
        output_path: Path to the output file.
        mapping_date: Date string (YYYY-MM-DD). Defaults to today.
    """
    # Avoid circular imports
    from sssom import write_tsv

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    msdf = mapping_set.to_sssom_dataframe(mapping_date)

    with output_path.open("w", encoding="utf-8") as f:
        write_tsv(msdf, f)


# =============================================================================
# Legacy Output Format Functions
# =============================================================================


def write_subject_ids(
    mapping_set: MappingSet,
    output_path: Path | str,
) -> None:
    """Write primary IDs to a TSV file (legacy format).

    Outputs a file with columns: primaryID, primarySymbol (if available).
    Excludes withdrawn entries (sec2pri:WithdrawnEntry).

    Args:
        mapping_set: The pysec2pri MappingSet to write.
        output_path: Path to the output file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Collect unique primary IDs with their labels/symbols
    seen: set[str] = set()
    rows: list[dict[str, str]] = []
    for m in mapping_set.mappings:
        if m.subject_id not in seen and m.subject_id != BaseMapping.withdrawn_entry:
            seen.add(m.subject_id)
            legacy = m.to_legacy_dict()
            primary_id = legacy["primaryID"] if legacy["primaryID"] is not None else ""
            symbol = legacy.get("primarySymbol") or legacy.get("primaryLabel") or ""
            row: dict[str, str] = {"primaryID": primary_id, "primarySymbol": symbol}
            rows.append(row)

    # Determine columns - check if any row has non-empty symbol
    has_symbol = any(r.get("primarySymbol") for r in rows)
    columns = ["primaryID", "primarySymbol"] if has_symbol else ["primaryID"]

    with output_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            values = [row.get(col, "") or "" for col in columns]
            f.write("\t".join(values) + "\n")


def write_sec2pri(
    mapping_set: MappingSet,
    output_path: Path | str,
) -> None:
    """Write secondary to primary ID mappings to a TSV file (legacy format).

    Outputs the raw mapping data similar to the legacy R scripts.

    Args:
        mapping_set: The pysec2pri MappingSet to write.
        output_path: Path to the output file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    has_symbol = False
    for m in mapping_set.mappings:
        legacy = m.to_legacy_dict()
        # Ensure all values are str, never None
        row = {k: (v if v is not None else "") for k, v in legacy.items()}
        # Check if any symbol data is present (non-empty)
        if row.get("primarySymbol") or row.get("secondarySymbol"):
            has_symbol = True
        rows.append(row)

    # Determine columns based on data
    if has_symbol:
        columns = [
            "primaryID",
            "primarySymbol",
            "secondaryID",
            "secondarySymbol",
            "predicateID",
            "mapping_cardinality_sec2pri",
            "comment",
            "source",
        ]
    else:
        columns = [
            "primaryID",
            "secondaryID",
            "predicateID",
            "mapping_cardinality_sec2pri",
            "comment",
            "source",
        ]

    with output_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            values = [row.get(col, "") or "" for col in columns]
            f.write("\t".join(values) + "\n")


def write_name2synonym(
    mapping_set: MappingSet,
    output_path: Path | str,
) -> None:
    """Write name to synonym mappings to a TSV file (legacy format).

    For use with ChEBI and HMDB which have name/synonym data.

    Args:
        mapping_set: The pysec2pri MappingSet to write.
        output_path: Path to the output file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    columns = ["primaryID", "name", "synonym"]
    rows: list[dict[str, str]] = []

    for m in mapping_set.mappings:
        legacy = m.to_legacy_dict()
        # Get name and synonym from label fields
        name = legacy.get("primaryLabel") or ""
        synonym = legacy.get("secondaryLabel") or ""
        if name or synonym:
            rows.append(
                {
                    "primaryID": legacy["primaryID"] if legacy["primaryID"] is not None else "",
                    "name": name,
                    "synonym": synonym,
                }
            )

    with output_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            values = [row.get(col, "") for col in columns]
            f.write("\t".join(values) + "\n")


def write_symbol2prev(
    mapping_set: MappingSet,
    output_path: Path | str,
) -> None:
    """Write symbol to alias/previous symbol mappings (legacy format).

    For use with HGNC and NCBI which have symbol data.

    Args:
        mapping_set: The pysec2pri MappingSet to write.
        output_path: Path to the output file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    columns = [
        "primaryID",
        "primarySymbol",
        "secondarySymbol",
        "predicateID",
        "mapping_cardinality_sec2pri",
        "comment",
        "source",
    ]
    rows: list[dict[str, str]] = []

    for m in mapping_set.mappings:
        legacy = m.to_legacy_dict()
        # Only include rows that have symbol data (non-empty)
        if legacy.get("primarySymbol") or legacy.get("secondarySymbol"):
            rows.append(
                {
                    "primaryID": legacy["primaryID"] if legacy["primaryID"] is not None else "",
                    "primarySymbol": legacy.get("primarySymbol") or "",
                    "secondarySymbol": legacy.get("secondarySymbol") or "",
                    "predicateID": (
                        legacy["predicateID"] if legacy["predicateID"] is not None else ""
                    ),
                    "mapping_cardinality_sec2pri": legacy.get("mapping_cardinality_sec2pri") or "",
                    "comment": legacy.get("comment") or "",
                    "source": legacy.get("source") or "",
                }
            )

    with output_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            values = [row.get(col, "") for col in columns]
            f.write("\t".join(values) + "\n")


# =============================================================================
# Multi-output helper
# =============================================================================


def resolve_output_types(
    *,
    output_types: list[str] | None,
    config: DatasourceConfig | None,
) -> list[OutputType]:
    """Resolve requested output types into concrete OutputType values.

    If output_types is None, defaults to all outputs available for the
    datasource config. If no config is available, defaults to SSSOM only.

    String values are resolved first by enum value, then by case-insensitive
    enum name match. Unknown values are silently dropped.

    Args:
        output_types: User-requested output type identifiers, or None.
        config: Datasource configuration object, or None.

    Returns:
        List of resolved OutputType values, filtered to those supported
        by the datasource if a config is provided.
    """
    if output_types is None:
        if config:
            return list(config.available_outputs)
        return [OutputType.SSSOM]

    resolved: list[OutputType] = []

    for ot_name in output_types:
        try:
            resolved.append(OutputType(ot_name))
            continue
        except ValueError:
            pass

        for ot in OutputType:
            if ot.name.lower() == ot_name.lower():
                resolved.append(ot)
                break

    if config:
        resolved = [t for t in resolved if t in config.available_outputs]

    return resolved


def write_output(
    *,
    mapping_set: MappingSet,
    output_type: str,
    output_path: Path,
    mapping_date: str | None,
) -> None:
    """Write a single output file for a specific output type.

    Dispatches to the correct writer implementation based on output_type.
    Raises if an unsupported output type is provided.

    Args:
        mapping_set: MappingSet to serialize.
        output_type: OutputType specifying which output to write.
        output_path: Destination file path.
        mapping_date: Optional mapping date for formats that require it.

    Returns:
        None

    Raises:
        ValueError: If output_type is not handled explicitly.
    """
    if output_type == OutputType.SSSOM:
        write_sssom(mapping_set, output_path, mapping_date=mapping_date)
    elif output_type == OutputType.PRIMARY_IDS:
        write_subject_ids(mapping_set, output_path)
    elif output_type == OutputType.SEC2PRI:
        write_sec2pri(mapping_set, output_path)
    elif output_type == OutputType.NAME2SYNONYM:
        write_name2synonym(mapping_set, output_path)
    elif output_type == OutputType.SYMBOL2PREV:
        write_symbol2prev(mapping_set, output_path)
    else:
        raise ValueError(f"Unhandled output type: {output_type}")


def resolve_output_path(
    *,
    output_dir: Path,
    datasource_name: str,
    output_type: OutputType,
    config: DatasourceConfig,
) -> Path:
    """Determine the filesystem path for a given output type.

    Uses datasource configuration to generate filenames when available.
    Falls back to a default naming scheme otherwise.

    Args:
        output_dir: Directory where outputs are written.
        datasource_name: Datasource name used for fallback filenames.
        output_type: OutputType being written.
        config: Datasource configuration object, or None.

    Returns:
        Full Path to the output file.
    """
    if config:
        filename = get_output_filename(config, output_type)
    else:
        filename = f"{datasource_name}_{output_type.value}.tsv"

    return output_dir / filename


def write_outputs(
    mapping_set: MappingSet,
    output_dir: Path | str,
    output_types: list[str] | None = None,
    datasource_name: str | None = None,
    mapping_date: str | None = None,
) -> dict[str, Path]:
    """Write one or more output files for a MappingSet.

    Args:
        mapping_set: The MappingSet to write.
        output_dir: Directory to write output files to.
        output_types: Optional list of output type identifiers.
        datasource_name: Optional datasource name override.
        mapping_date: Optional mapping date string for outputs that require it.

    Returns:
        Dictionary mapping output type values to their written file paths.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ds_name = datasource_name or mapping_set.datasource_name
    config = get_datasource_config(ds_name)
    if config is None:
        raise ValueError(f"No datasource config found for '{ds_name}'")

    types_to_write = resolve_output_types(
        output_types=output_types,
        config=config,
    )

    written_files: dict[str, Path] = {}

    for output_type in types_to_write:
        output_path = resolve_output_path(
            output_dir=output_dir,
            datasource_name=ds_name,
            output_type=output_type,
            config=config,
        )

        write_output(
            mapping_set=mapping_set,
            output_type=output_type,
            output_path=output_path,
            mapping_date=mapping_date,
        )

        written_files[output_type.value] = output_path

    return written_files
