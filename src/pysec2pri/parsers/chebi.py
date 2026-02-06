"""ChEBI parser for secondary-to-primary identifier mappings.

Supports two formats:
1. TSV flat files (releases >= 245): secondary_ids.tsv, names.tsv, compounds.tsv
2. SDF files (releases < 245): chebi_3_stars.sdf or ChEBI_complete.sdf

Extracts:
1. ID-to-ID mappings: secondary ChEBI IDs -> primary ChEBI IDs
2. Label-to-label mappings: synonyms -> primary names

Uses SSSOM-compliant MappingSet classes with cardinality computation.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from sssom_schema import Mapping

from pysec2pri.parsers.base import (
    BaseParser,
    Sec2PriMappingSet,
)

if TYPE_CHECKING:
    pass

# Threshold version where TSV format was introduced
NEW_FORMAT_VERSION = 245


class ChEBIParser(BaseParser):
    """Parser for ChEBI data files.

    Supports both TSV flat files (>= release 245) and legacy SDF files.
    Extracts secondary-to-primary ChEBI identifier mappings and
    name-to-synonym relationships.

    Returns an IdMappingSet for ID mappings (cardinality computed on IDs)
    and can optionally include synonym mappings via LabelMappingSet.
    """

    datasource_name = "chebi"

    def __init__(
        self,
        version: str | None = None,
        show_progress: bool = True,
        subset: str = "3star",
    ) -> None:
        """Initialize the ChEBI parser.

        Args:
            version: Version/release identifier for the datasource.
            show_progress: Whether to show progress bars during parsing.
            subset: "3star" or "complete" - which compound subset to use.
                For TSV format, filters by stars in compounds.tsv.
                For SDF format, determines which file to download.
        """
        super().__init__(version=version, show_progress=show_progress)
        self.subset = subset

    @property
    def source_url(self) -> str:
        """Get the default download URL from config."""
        return self.get_download_url("sdf") or ""

    def _is_new_format(self) -> bool:
        """Check if we should use new TSV format based on version."""
        if self.version is None:
            return True  # Default to new format for latest
        try:
            return int(self.version) >= NEW_FORMAT_VERSION
        except ValueError:
            return True  # Default to new if version is not numeric

    def parse(
        self,
        input_path: Path | str | None = None,
        *,
        secondary_ids_path: Path | str | None = None,
        compounds_path: Path | str | None = None,
    ) -> Sec2PriMappingSet:
        """Parse ChEBI data into an IdMappingSet.

        For TSV format (>= 245): provide secondary_ids_path and compounds_path
        For SDF format (< 245): provide input_path to the SDF file

        Args:
            input_path: Path to the ChEBI SDF file (legacy format).
            secondary_ids_path: Path to secondary_ids.tsv (new format).
            compounds_path: Path to compounds.tsv for 3-star filtering.

        Returns:
            IdMappingSet with computed cardinalities.
        """
        if secondary_ids_path is not None:
            # New TSV format
            raw_mappings = _parse_secondary_ids_tsv(
                Path(secondary_ids_path),
                compounds_path=Path(compounds_path) if compounds_path else None,
                subset=self.subset,
                show_progress=self.show_progress,
            )
        elif input_path is not None:
            # Legacy SDF format
            raw_mappings, _ = _parse_chebi_sdf_fast(
                Path(input_path),
                show_progress=self.show_progress,
            )
        else:
            raise ValueError("Must provide either input_path (SDF) or secondary_ids_path (TSV)")

        # Build SSSOM Mapping objects for ID mappings
        mappings = self._build_id_mappings(raw_mappings)

        # Create IdMappingSet and compute cardinalities
        return self._create_mapping_set(mappings, mapping_type="id")

    def parse_synonyms(
        self,
        input_path: Path | str | None = None,
        *,
        names_path: Path | str | None = None,
        compounds_path: Path | str | None = None,
    ) -> Sec2PriMappingSet:
        """Parse ChEBI data into a LabelMappingSet for synonyms.

        For TSV format (>= 245): provide names_path and compounds_path
        For SDF format (< 245): provide input_path to the SDF file

        Args:
            input_path: Path to the ChEBI SDF file (legacy format).
            names_path: Path to names.tsv (new format).
            compounds_path: Path to compounds.tsv for 3-star filtering.

        Returns:
            LabelMappingSet with computed cardinalities based on labels.
        """
        if names_path is not None:
            # New TSV format
            raw_mappings = _parse_names_tsv(
                Path(names_path),
                compounds_path=Path(compounds_path) if compounds_path else None,
                subset=self.subset,
                show_progress=self.show_progress,
            )
        elif input_path is not None:
            # Legacy SDF format
            _, raw_mappings = _parse_chebi_sdf_fast(
                Path(input_path),
                show_progress=self.show_progress,
            )
        else:
            raise ValueError("Must provide either input_path (SDF) or names_path (TSV)")

        # Build SSSOM Mapping objects for label mappings
        mappings = self._build_label_mappings(raw_mappings)

        # Create LabelMappingSet and compute cardinalities
        return self._create_mapping_set(mappings, mapping_type="label")

    def _build_id_mappings(self, raw_id_mappings: list[tuple[str, str]]) -> list[Mapping]:
        """Build SSSOM Mapping objects for ID-to-ID mappings.

        Args:
            raw_id_mappings: List of (primary_id, secondary_id) tuples.

        Returns:
            List of sssom_schema.Mapping objects.
        """
        m_meta = self.get_mapping_metadata()
        mappings: list[Mapping] = []

        for primary_id, secondary_id in self._progress(
            raw_id_mappings, desc="Creating ID mappings"
        ):
            # SSSOM: subject_id = primary (current), object_id = secondary (old)
            mapping = Mapping(
                subject_id=primary_id,
                object_id=secondary_id,
                predicate_id=m_meta["predicate_id"],
                mapping_justification=m_meta["mapping_justification"],
                subject_source=m_meta.get("subject_source"),
                object_source=m_meta.get("object_source"),
                mapping_tool=m_meta.get("mapping_tool"),
                confidence=m_meta.get("confidence"),
                license=m_meta.get("license"),
            )
            mappings.append(mapping)

        return mappings

    def _build_label_mappings(self, raw_name_mappings: list[tuple[str, str, str]]) -> list[Mapping]:
        """Build SSSOM Mapping objects for label-to-label (synonym) mappings.

        Args:
            raw_name_mappings: List of (subject_id, primary_name, synonym).

        Returns:
            List of sssom_schema.Mapping objects.
        """
        m_meta = self.get_mapping_metadata()
        mappings: list[Mapping] = []

        for subject_id, primary_name, synonym in self._progress(
            raw_name_mappings, desc="Creating synonym mappings"
        ):
            # For synonyms: subject_label = primary name, object_label = synonym
            mapping = Mapping(
                subject_id=subject_id,
                subject_label=primary_name,
                object_id=subject_id,  # Same entity
                object_label=synonym,
                predicate_id="oboInOwl:hasRelatedSynonym",
                mapping_justification=m_meta["mapping_justification"],
                subject_source=m_meta.get("subject_source"),
                object_source=m_meta.get("object_source"),
                mapping_tool=m_meta.get("mapping_tool"),
                license=m_meta.get("license"),
            )
            mappings.append(mapping)

        return mappings

    def _create_mapping_set(
        self, mappings: list[Mapping], mapping_type: str = "id"
    ) -> Sec2PriMappingSet:
        """Create an IdMappingSet or LabelMappingSet with metadata from config.

        Delegates to BaseParser.create_mapping_set().
        """
        return self.create_mapping_set(mappings, mapping_type)


# =============================================================================
# TSV parsing functions (new format >= 245)
# =============================================================================


def _get_3star_compound_ids(
    compounds_path: Path,
    show_progress: bool = True,
) -> set[int]:
    """Get set of compound IDs with 3 stars from compounds.tsv.

    Args:
        compounds_path: Path to compounds.tsv file.
        show_progress: Whether to show progress.

    Returns:
        Set of compound IDs (as integers) with stars == 3.
    """
    import polars as pl

    if show_progress:
        pass

    df = pl.read_csv(
        compounds_path,
        separator="\t",
        columns=["id", "stars"],
        schema_overrides={"id": pl.Int64, "stars": pl.Int64},
    )

    three_star_ids = set(df.filter(pl.col("stars") == 3)["id"].to_list())

    if show_progress:
        pass

    return three_star_ids


def _parse_secondary_ids_tsv(
    secondary_ids_path: Path,
    compounds_path: Path | None = None,
    subset: str = "3star",
    show_progress: bool = True,
) -> list[tuple[str, str]]:
    """Parse secondary_ids.tsv into (primary_id, secondary_id) tuples.

    Format: compound_id | secondary_id

    Args:
        secondary_ids_path: Path to secondary_ids.tsv file.
        compounds_path: Path to compounds.tsv for 3-star filtering.
        subset: "3star" or "complete".
        show_progress: Whether to show progress.

    Returns:
        List of (primary_curie, secondary_curie) tuples.
    """
    import polars as pl

    if show_progress:
        pass

    df = pl.read_csv(
        secondary_ids_path,
        separator="\t",
        schema_overrides={"compound_id": pl.Int64, "secondary_id": pl.Int64},
    )

    # Filter to 3-star compounds if requested
    if subset == "3star" and compounds_path is not None:
        three_star_ids = _get_3star_compound_ids(compounds_path, show_progress)
        df = df.filter(pl.col("compound_id").is_in(three_star_ids))

    # Build mapping tuples with CHEBI: prefix
    mappings: list[tuple[str, str]] = []
    for row in df.iter_rows():
        primary_id = f"CHEBI:{row[0]}"
        secondary_id = f"CHEBI:{row[1]}"
        mappings.append((primary_id, secondary_id))

    if show_progress:
        pass

    return mappings


def _parse_names_tsv(
    names_path: Path,
    compounds_path: Path | None = None,
    subset: str = "3star",
    show_progress: bool = True,
) -> list[tuple[str, str, str]]:
    """Parse names.tsv into (subject_id, primary_name, synonym) tuples.

    Format: id|compound_id|name|type|status_id|adapted|language_code|ascii_name

    For each compound, the name with the smallest id is considered the primary
    name, and all other names are synonyms.

    Args:
        names_path: Path to names.tsv file.
        compounds_path: Path to compounds.tsv for 3-star filtering.
        subset: "3star" or "complete".
        show_progress: Whether to show progress.

    Returns:
        List of (subject_curie, primary_name, synonym) tuples.
    """
    import polars as pl

    if show_progress:
        pass

    df = pl.read_csv(
        names_path,
        separator="\t",
        columns=["id", "compound_id", "name"],
        schema_overrides={"id": pl.Int64, "compound_id": pl.Int64, "name": pl.Utf8},
    )

    # Filter to 3-star compounds if requested
    if subset == "3star" and compounds_path is not None:
        three_star_ids = _get_3star_compound_ids(compounds_path, show_progress)
        df = df.filter(pl.col("compound_id").is_in(three_star_ids))

    # For each compound, find the primary name (smallest id)
    # and all other names as synonyms
    primary_names = df.group_by("compound_id").agg(
        [
            pl.col("name").sort_by("id").first().alias("primary_name"),
        ]
    )

    # Join to get primary name for each row
    df_with_primary = df.join(primary_names, on="compound_id")

    # Filter to only synonym rows (where name != primary_name)
    synonyms_df = df_with_primary.filter(pl.col("name") != pl.col("primary_name"))

    # Build mapping tuples
    mappings: list[tuple[str, str, str]] = []
    for row in synonyms_df.iter_rows():
        compound_id = row[1]
        name = row[2]
        primary_name = row[3]
        subject_curie = f"CHEBI:{compound_id}"
        mappings.append((subject_curie, primary_name, name))

    if show_progress:
        pass

    return mappings


# =============================================================================
# SDF parsing functions (legacy format < 245)
# =============================================================================


def _read_sdf_lines(input_path: Path) -> list[str]:
    """Read an SDF file into a list of lines."""
    with input_path.open("r", encoding="utf-8") as f:
        return f.readlines()


def _extract_value(lines: list[str], i: int, total_lines: int) -> tuple[str | None, int]:
    """Extract a value from the next line after a tag."""
    i += 1
    if i < total_lines:
        value = lines[i].strip()
        return value, i
    return None, i


def _extract_secondary_ids(
    lines: list[str],
    i: int,
    total_lines: int,
    current_subject_id: str | None,
    raw_id_mappings: list[tuple[str, str]],
) -> int:
    """Extract secondary IDs from the lines."""
    i += 1
    while i < total_lines:
        sec_line = lines[i].strip()
        if sec_line.startswith("CHEBI:"):
            if current_subject_id:
                # Handle semicolon-separated IDs on same line
                for sec_id in sec_line.split(";"):
                    sec_id = sec_id.strip()
                    if sec_id.startswith("CHEBI:"):
                        raw_id_mappings.append((current_subject_id, sec_id))
            i += 1
        else:
            break
    return i


def _extract_synonyms(
    lines: list[str],
    i: int,
    total_lines: int,
    current_subject_id: str | None,
    current_name: str | None,
    raw_name_mappings: list[tuple[str, str, str]],
) -> int:
    """Extract synonyms from the lines."""
    i += 1
    while i < total_lines:
        syn_line = lines[i].strip()
        if syn_line and not syn_line.startswith(">"):
            if current_subject_id and current_name:
                raw_name_mappings.append((current_subject_id, current_name, syn_line))
            i += 1
        else:
            break
    return i


def _process_sdf_line(
    lines: list[str],
    i: int,
    total_lines: int,
    current_subject_id: str | None,
    current_name: str | None,
    raw_id_mappings: list[tuple[str, str]],
    raw_name_mappings: list[tuple[str, str, str]],
) -> tuple[int, str | None, str | None, bool]:
    """Process a single SDF line and return updated state."""
    line = lines[i].strip()

    if "<chebi id>" in line.lower():
        value, i = _extract_value(lines, i, total_lines)
        if value:
            # Handle both "CHEBI:123" and "123" formats
            if value.startswith("CHEBI:"):
                value = f"CHEBI:{value[6:]}"  # Normalize
            else:
                value = f"CHEBI:{value}"
        return i, value, current_name, False

    if "<chebi name>" in line.lower():
        value, i = _extract_value(lines, i, total_lines)
        return i, current_subject_id, value, False

    if "<secondary" in line.lower():
        i = _extract_secondary_ids(lines, i, total_lines, current_subject_id, raw_id_mappings)
        return i, current_subject_id, current_name, True

    if "<synonym" in line.lower():
        i = _extract_synonyms(
            lines,
            i,
            total_lines,
            current_subject_id,
            current_name,
            raw_name_mappings,
        )
        return i, current_subject_id, current_name, True

    if line.startswith("$$$$"):
        return i, None, None, False

    return i, current_subject_id, current_name, False


def _parse_chebi_sdf_fast(
    input_path: Path,
    show_progress: bool = True,
) -> tuple[list[tuple[str, str]], list[tuple[str, str, str]]]:
    """Parse ChEBI SDF file into raw ID and synonym mappings.

    Args:
        input_path: Path to the SDF file.
        show_progress: Whether to display a progress bar.

    Returns:
        Tuple of:
            - list of (primary_id, secondary_id) for ID mappings
            - list of (subject_id, primary_name, synonym) for label mappings
    """
    from tqdm import tqdm

    lines = _read_sdf_lines(input_path)

    raw_id_mappings: list[tuple[str, str]] = []
    raw_name_mappings: list[tuple[str, str, str]] = []

    current_subject_id: str | None = None
    current_name: str | None = None

    i = 0
    total_lines = len(lines)

    pbar = tqdm(total=total_lines, desc="Parsing ChEBI SDF") if show_progress else None

    while i < total_lines:
        i, current_subject_id, current_name, skip = _process_sdf_line(
            lines,
            i,
            total_lines,
            current_subject_id,
            current_name,
            raw_id_mappings,
            raw_name_mappings,
        )
        if not skip:
            i += 1
        if pbar:
            pbar.update(1)

    if pbar:
        pbar.close()

    if show_progress:
        pass

    return raw_id_mappings, raw_name_mappings


__all__ = ["ChEBIParser"]
