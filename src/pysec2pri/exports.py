"""Export functions for writing mapping sets to various file formats.

This module provides functions for writing MappingSet objects to:
- SSSOM TSV format (standard semantic mapping format)
- Legacy TSV formats matching the original R script outputs
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from pysec2pri.models import BaseMapping

if TYPE_CHECKING:
    from pysec2pri.models import MappingSet

__all__ = [
    # SSSOM output
    "write_sssom",
    # Legacy output formats
    "write_primary_ids",
    "write_sec2pri",
    "write_name2synonym",
    "write_symbol2prev",
    # Multi-output helper
    "write_outputs",
]


# =============================================================================
# SSSOM Output
# =============================================================================


def write_sssom(
    mapping_set: "MappingSet",
    output_path: Path | str,
    mapping_date: str | None = None,
) -> None:
    """Write a MappingSet to an SSSOM TSV file.

    Uses the official sssom-py library for correct SSSOM format output.

    Args:
        mapping_set: The pysec2pri MappingSet to write.
        output_path: Path to the output file.
        mapping_date: Date string (YYYY-MM-DD). Defaults to today.

    Example:
        >>> from pysec2pri import parse_chebi
        >>> from pysec2pri.exports import write_sssom
        >>> mapping_set = parse_chebi("ChEBI_complete_3star.sdf")
        >>> write_sssom(mapping_set, "chebi_mappings.sssom.tsv")
    """
    # Import here to avoid circular imports
    from pysec2pri.api import to_sssom_dataframe
    from sssom import write_tsv

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    msdf = to_sssom_dataframe(mapping_set, mapping_date)

    with output_path.open("w", encoding="utf-8") as f:
        write_tsv(msdf, f)


# =============================================================================
# Legacy Output Format Functions
# =============================================================================


def write_primary_ids(
    mapping_set: "MappingSet",
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
        if (m.primary_id not in seen
                and m.primary_id != BaseMapping.WITHDRAWN_SENTINEL):
            seen.add(m.primary_id)
            legacy = m.to_legacy_dict()
            row = {"primaryID": legacy["primaryID"]}
            # Get symbol/label from legacy dict
            if "primarySymbol" in legacy and legacy["primarySymbol"]:
                row["primarySymbol"] = legacy["primarySymbol"]
            elif "primaryLabel" in legacy and legacy["primaryLabel"]:
                row["primarySymbol"] = legacy["primaryLabel"]
            rows.append(row)

    # Determine columns
    has_symbol = any("primarySymbol" in r for r in rows)
    columns = ["primaryID", "primarySymbol"] if has_symbol else ["primaryID"]

    with output_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            values = [row.get(col, "") or "" for col in columns]
            f.write("\t".join(values) + "\n")


def write_sec2pri(
    mapping_set: "MappingSet",
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
        row = m.to_legacy_dict()
        if "primarySymbol" in row:
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
    mapping_set: "MappingSet",
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
        name = legacy.get("primaryLabel", "")
        synonym = legacy.get("secondaryLabel", "")
        if name or synonym:
            rows.append({
                "primaryID": legacy["primaryID"],
                "name": name,
                "synonym": synonym,
            })

    with output_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            values = [row.get(col, "") for col in columns]
            f.write("\t".join(values) + "\n")


def write_symbol2prev(
    mapping_set: "MappingSet",
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
        if "primarySymbol" in legacy:
            rows.append({
                "primaryID": legacy["primaryID"],
                "primarySymbol": legacy["primarySymbol"],
                "secondarySymbol": legacy["secondarySymbol"],
                "predicateID": legacy["predicateID"],
                "mapping_cardinality_sec2pri": legacy[
                    "mapping_cardinality_sec2pri"
                ],
                "comment": legacy["comment"],
                "source": legacy["source"],
            })

    with output_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            values = [row.get(col, "") for col in columns]
            f.write("\t".join(values) + "\n")


# =============================================================================
# Multi-output helper
# =============================================================================


def write_outputs(
    mapping_set: "MappingSet",
    output_dir: Path | str,
    output_types: list[str] | None = None,
    datasource_name: str | None = None,
    mapping_date: str | None = None,
) -> dict[str, Path]:
    """Write multiple output files based on requested output types.

    Args:
        mapping_set: The pysec2pri MappingSet to write.
        output_dir: Directory to write output files to.
        output_types: List of output type names to generate.
            If None, generates all available outputs for the datasource.
        datasource_name: Name of the datasource (for filename generation).
            If not provided, uses mapping_set.datasource_name.
        mapping_date: Date string for SSSOM output.

    Returns:
        Dictionary mapping output type names to their file paths.
    """
    from pysec2pri.constants import (
        OutputType,
        get_datasource_config,
        get_output_filename,
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ds_name = datasource_name or mapping_set.datasource_name

    config = get_datasource_config(ds_name)

    # Determine which outputs to generate
    if output_types is None:
        # Default to all available outputs for this datasource
        if config:
            types_to_write = list(config.available_outputs)
        else:
            types_to_write = [OutputType.SSSOM]
    else:
        # Parse the requested output types
        types_to_write = []
        for ot_name in output_types:
            try:
                types_to_write.append(OutputType(ot_name))
            except ValueError:
                # Try to match by name
                for ot in OutputType:
                    if ot.name.lower() == ot_name.lower():
                        types_to_write.append(ot)
                        break

    # Filter to only available outputs for this datasource
    if config:
        types_to_write = [
            t for t in types_to_write if t in config.available_outputs
        ]

    written_files: dict[str, Path] = {}

    for output_type in types_to_write:
        if config:
            filename = get_output_filename(config, output_type)
        else:
            filename = f"{ds_name}_{output_type.value}.tsv"
        output_path = output_dir / filename

        if output_type == OutputType.SSSOM:
            write_sssom(mapping_set, output_path, mapping_date=mapping_date)
        elif output_type == OutputType.PRIMARY_IDS:
            write_primary_ids(mapping_set, output_path)
        elif output_type == OutputType.SEC2PRI:
            write_sec2pri(mapping_set, output_path)
        elif output_type == OutputType.NAME2SYNONYM:
            write_name2synonym(mapping_set, output_path)
        elif output_type == OutputType.SYMBOL2PREV:
            write_symbol2prev(mapping_set, output_path)

        written_files[output_type.value] = output_path

    return written_files
