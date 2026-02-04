"""ChEBI SDF file parser for secondary-to-primary identifier mappings.

This parser extracts:
1. ID-to-ID mappings: secondary ChEBI IDs -> primary ChEBI IDs
2. Label-to-label mappings: synonyms -> primary names

Uses SSSOM-compliant MappingSet classes with cardinality computation.
"""

from __future__ import annotations

from pathlib import Path

from sssom_schema import Mapping

from pysec2pri.parsers.base import (
    BaseParser,
    Sec2PriMappingSet,
)


class ChEBIParser(BaseParser):
    """Parser for ChEBI SDF files.

    Extracts secondary-to-primary ChEBI identifier mappings and
    name-to-synonym relationships from ChEBI SDF format files.

    Returns an IdMappingSet for ID mappings (cardinality computed on IDs)
    and can optionally include synonym mappings via LabelMappingSet.
    """

    datasource_name = "chebi"

    def __init__(self, version: str | None = None, show_progress: bool = True) -> None:
        """Initialize the ChEBI parser.

        Args:
            version: Version/release identifier for the datasource.
            show_progress: Whether to show progress bars during parsing.
        """
        super().__init__(version=version, show_progress=show_progress)

    @property
    def source_url(self) -> str:
        """Get the default SDF download URL from config."""
        return self.get_download_url("sdf") or ""

    def parse(self, input_path: Path | str | None) -> Sec2PriMappingSet:
        """Parse a ChEBI SDF file into an IdMappingSet.

        Args:
            input_path: Path to the ChEBI SDF file.

        Returns:
            IdMappingSet with computed cardinalities based on subject_id/object_id.
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        input_path = Path(str(input_path))
        lines = _read_sdf_lines(input_path)
        raw_id_mappings, _ = _parse_chebi_sdf(lines, show_progress=self.show_progress)

        # Build SSSOM Mapping objects for ID mappings
        mappings = self._build_id_mappings(raw_id_mappings)

        # Create IdMappingSet and compute cardinalities
        mapping_set = self._create_mapping_set(mappings, mapping_type="id")
        return mapping_set

    def parse_synonyms(self, input_path: Path | str | None) -> Sec2PriMappingSet:
        """Parse a ChEBI SDF file into a LabelMappingSet for synonyms.

        Args:
            input_path: Path to the ChEBI SDF file.

        Returns:
            LabelMappingSet with computed cardinalities based on labels.
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        input_path = Path(str(input_path))
        lines = _read_sdf_lines(input_path)
        _, raw_name_mappings = _parse_chebi_sdf(lines, show_progress=self.show_progress)

        # Build SSSOM Mapping objects for label mappings
        mappings = self._build_label_mappings(raw_name_mappings)

        # Create LabelMappingSet and compute cardinalities
        mapping_set = self._create_mapping_set(mappings, mapping_type="label")
        return mapping_set

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
            raw_name_mappings: List of (subject_id, primary_name, synonym) tuples.

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
# Private helper functions for SDF parsing
# =============================================================================


def _read_sdf_lines(input_path: Path) -> list[str]:
    """Read an SDF file into a list of lines.

    Args:
        input_path: Path to the SDF file.

    Returns:
        List of lines including newlines.
    """
    with input_path.open("r", encoding="utf-8") as f:
        return f.readlines()


def _extract_value(lines: list[str], i: int, total_lines: int) -> tuple[str | None, int]:
    """Extract a value from the next line after a tag.

    Args:
        lines: Lines from the SDF file.
        i: Current index in the lines.
        total_lines: Total number of lines.

    Returns:
        Tuple of the extracted value (or None) and the updated index.
    """
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
    """Extract secondary IDs from the lines.

    Args:
        lines: Lines from the SDF file.
        i: Current index in the lines.
        total_lines: Total number of lines.
        current_subject_id: Current primary ID.
        raw_id_mappings: List to store raw ID mappings.

    Returns:
        Updated index.
    """
    i += 1
    while i < total_lines:
        sec_line = lines[i].strip()
        if sec_line.startswith("CHEBI:"):
            if current_subject_id:
                for sec_id in sec_line.split(";"):
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
    """Extract synonyms from the lines.

    Args:
        lines: Lines from the SDF file.
        i: Current index in the lines.
        total_lines: Total number of lines.
        current_subject_id: Current primary ID.
        current_name: Current name.
        raw_name_mappings: List to store raw name mappings.

    Returns:
        Updated index.
    """
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
    """Process a single SDF line and return updated state.

    Args:
        lines: Lines from the SDF file.
        i: Current index in the lines.
        total_lines: Total number of lines.
        current_subject_id: Current primary ID.
        current_name: Current name.
        raw_id_mappings: List to store raw ID mappings.
        raw_name_mappings: List to store raw name mappings.

    Returns:
        Tuple of (new index, subject_id, name, skip_increment).
    """
    line = lines[i].strip()

    if "<chebi id>" in line.lower():
        value, i = _extract_value(lines, i, total_lines)
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


def _parse_chebi_sdf(
    lines: list[str],
    *,
    show_progress: bool,
) -> tuple[list[tuple[str, str]], list[tuple[str, str, str]]]:
    """Parse ChEBI SDF lines into raw ID and synonym mappings.

    Args:
        lines: Lines from the SDF file.
        show_progress: Whether to display a progress bar.

    Returns:
        Tuple of:
            - list of (primary_id, secondary_id) for ID mappings
            - list of (subject_id, primary_name, synonym) for label mappings
    """
    from tqdm import tqdm

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

    return raw_id_mappings, raw_name_mappings


__all__ = ["ChEBIParser"]
