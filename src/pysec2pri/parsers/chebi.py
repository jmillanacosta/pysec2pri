"""ChEBI SDF file parser for secondary-to-primary identifier mappings."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pysec2pri.models import ChEBIMapping, MappingCardinality, MappingSet
from pysec2pri.parsers.base import BaseParser


# Utility for progress bar
class _Progressor:  # TODO move to base
    """Class for tqdm progress bar."""

    def __init__(self, show_progress: bool):
        """Init."""
        self.show_progress = show_progress

    def progress(self, iterable: Any, desc: Any = None, total: Any = None) -> Any:
        """Wrap tqdm."""
        if self.show_progress:
            from tqdm import tqdm

            return tqdm(iterable, desc=desc, total=total)
        return iterable


class ChEBIParser(BaseParser):
    """Parser for ChEBI SDF files.

    Extracts secondary-to-primary ChEBI identifier mappings and
    name-to-synonym relationships from ChEBI SDF format files.
    """

    datasource_name = "ChEBI"
    default_source_url = "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_complete_3star.sdf"

    def __init__(self, version: str | None = None, show_progress: bool = True) -> None:
        """Initialize the ChEBI parser.

        Args:
            version: Version/release identifier for the datasource.
            show_progress: Whether to show progress bars during parsing.
        """
        super().__init__(version=version, show_progress=show_progress)

    def parse(self, input_path: Path | str | None) -> MappingSet:
        """Parse a ChEBI SDF file into a MappingSet."""
        if input_path is None:
            raise ValueError("input_path must not be None")
        input_path = Path(str(input_path))
        lines = read_sdf_lines(input_path)
        raw_id_mappings, raw_name_mappings = parse_chebi_sdf(
            lines, show_progress=self.show_progress
        )
        mappings = []
        mappings.extend(
            build_id_mappings(
                raw_id_mappings,
                source_url=self.default_source_url,
                show_progress=self.show_progress,
            )
        )
        mappings.extend(
            build_synonym_mappings(
                raw_name_mappings,
                source_url=self.default_source_url,
                show_progress=self.show_progress,
            )
        )
        return build_mapping_set(
            mappings=mappings,
            datasource_name=self.datasource_name,
            version=self.version if self.version is not None else "",
        )


def read_sdf_lines(input_path: Path) -> list[str]:
    """Read an SDF file into a list of lines.

    Args:
        input_path: Path to the SDF file.

    Returns:
        List of lines including newlines.
    """
    with input_path.open("r", encoding="utf-8") as f:
        return f.readlines()


def extract_subject_id(
    lines: list[str], i: int, total_lines: int, pbar: object = None
) -> tuple[str | None, int]:
    """Extract the primary ID from the lines.

    Args:
        lines: Lines from the SDF file.
        i: Current index in the lines.
        total_lines: Total number of lines.
        pbar: Progress bar instance.

    Returns:
        Tuple of the primary ID (or None) and the updated index.
    """
    i += 1
    if i < total_lines:
        subject_id = lines[i].strip()
        return subject_id, i
    return None, i


def extract_object_ids(
    lines: list[str],
    i: int,
    total_lines: int,
    current_subject_id: str | None,
    raw_id_mappings: list[tuple[str, str]],
    pbar: object = None,
) -> int:
    """Extract secondary IDs from the lines.

    Args:
        lines: Lines from the SDF file.
        i: Current index in the lines.
        total_lines: Total number of lines.
        current_subject_id: Current primary ID.
        raw_id_mappings: List to store raw ID mappings.
        pbar: Progress bar instance.

    Returns:
        Updated index.
    """
    i += 1
    while i < total_lines:
        sec_line = lines[i].strip()
        if sec_line.startswith("CHEBI:"):
            if current_subject_id:
                raw_id_mappings.append((current_subject_id, sec_line))
            i += 1
        else:
            break
    return i


def extract_synonyms(
    lines: list[str],
    i: int,
    total_lines: int,
    current_subject_id: str | None,
    current_name: str | None,
    raw_name_mappings: list[tuple[str, str, str]],
    pbar: object = None,
) -> int:
    """Extract synonyms from the lines.

    Args:
        lines: Lines from the SDF file.
        i: Current index in the lines.
        total_lines: Total number of lines.
        current_subject_id: Current primary ID.
        current_name: Current name.
        raw_name_mappings: List to store raw name mappings.
        pbar: Progress bar instance.

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


def parse_chebi_sdf(
    lines: list[str],
    *,
    show_progress: bool,
) -> tuple[list[tuple[str, str]], list[tuple[str, str, str]]]:
    """Parse ChEBI SDF lines into raw ID and synonym mappings.

    Deprecated, now trying out polars (parse()).

    Args:
        lines: Lines from the SDF file.
        show_progress: Whether to display a progress bar.

    Returns:
        Tuple of:
            - list of (subject_id, object_id)
            - list of (subject_id, primary_name, synonym)
    """
    raw_id_mappings: list[tuple[str, str]] = []
    raw_name_mappings: list[tuple[str, str, str]] = []

    current_subject_id: str | None = None
    current_name: str | None = None

    i = 0
    total_lines = len(lines)
    progressor = _Progressor(show_progress)
    pbar = (
        progressor.progress(range(total_lines), desc="Parsing ChEBI SDF", total=total_lines)
        if show_progress
        else None
    )

    while i < total_lines:
        line = lines[i].strip()

        # ...existing code...

        if "<chebi id>" in line.lower():
            current_subject_id, i = extract_subject_id(lines, i, total_lines, pbar)

        elif "<chebi name>" in line.lower():
            current_name, i = extract_subject_id(lines, i, total_lines, pbar)

        elif "<secondary" in line.lower():
            i = extract_object_ids(lines, i, total_lines, current_subject_id, raw_id_mappings, pbar)
            continue

        elif "<synonym" in line.lower():
            i = extract_synonyms(
                lines,
                i,
                total_lines,
                current_subject_id,
                current_name,
                raw_name_mappings,
                pbar,
            )
            continue

        elif line.startswith("$$$$"):
            current_subject_id = None
            current_name = None

        i += 1

    # ...existing code...

    return raw_id_mappings, raw_name_mappings


def build_id_mappings(
    raw_id_mappings: list[tuple[str, str]],
    *,
    source_url: str,
    show_progress: bool,
) -> list[ChEBIMapping]:
    """Create ChEBIMapping objects for ID mappings.

    Args:
        raw_id_mappings: List of (subject_id, object_id).
        source_url: Source URL to attach to mappings.
        show_progress: Whether to show progress.

    Returns:
        List of ChEBIMapping objects with cardinality populated.
    """
    cardinality_map = BaseParser.compute_cardinality(raw_id_mappings)
    progressor = _Progressor(show_progress)
    iterator = progressor.progress(raw_id_mappings, desc="Creating ID mappings")

    mappings: list[ChEBIMapping] = []

    for subject_id, object_id in iterator:
        cardinality = cardinality_map.get(
            (subject_id, object_id),
            MappingCardinality.ONE_TO_ONE,
        )
        # Use IAO:0100001 (term replaced by) for clear replacements
        predicate = "IAO:0100001"
        if cardinality in (
            MappingCardinality.ONE_TO_MANY,
            MappingCardinality.MANY_TO_MANY,
            MappingCardinality.ONE_TO_ZERO,
        ):
            predicate = "oboInOwl:consider"
        mappings.append(
            ChEBIMapping(
                subject_id=subject_id,
                predicate_id=predicate,
                object_id=object_id,
                mapping_cardinality=cardinality,
                source_url=source_url,
            )
        )

    return mappings


def build_synonym_mappings(
    raw_name_mappings: list[tuple[str, str, str]],
    *,
    source_url: str,
    show_progress: bool,
) -> list[ChEBIMapping]:
    """Create ChEBIMapping objects for name/synonym mappings.

    Args:
        raw_name_mappings: List of (subject_id, name, synonym).
        source_url: Source URL to attach to mappings.
        show_progress: Whether to show progress.

    Returns:
        List of ChEBIMapping objects without cardinality.
    """
    progressor = _Progressor(show_progress)
    iterator = progressor.progress(
        raw_name_mappings,
        desc="Creating synonym mappings",
    )

    return [
        ChEBIMapping(
            subject_id=subject_id,
            predicate_id="oboInOwl:hasRelatedSynonym",
            object_id=synonym,
            subject_label=name,
            object_label=synonym,
            source_url=source_url,
        )
        for subject_id, name, synonym in iterator
    ]


def build_mapping_set(
    *,
    mappings: list[ChEBIMapping],
    datasource_name: str,
    version: str,
) -> MappingSet:
    """Assemble a MappingSet for ChEBI mappings.

    Args:
        mappings: All mapping objects.
        datasource_name: Datasource name.
        version: Datasource version.

    Returns:
        Fully populated MappingSet.
    """
    return MappingSet(
        mappings=mappings,
        datasource_name=datasource_name,
        version=version,
        mapping_set_id="omicsfixid_chebi_01",
        mapping_set_description=(
            "Secondary to primary ID mappings for ChEBI database, "
            "generated for the omicsFixID project."
        ),
        comment=(
            "object_label represents a synonym name by which the subject_label is also known."
        ),
        curie_map={
            "CHEBI:": "http://purl.obolibrary.org/obo/CHEBI_",
            "IAO:": "http://purl.obolibrary.org/obo/IAO_",
            "oboInOwl:": "http://www.geneontology.org/formats/oboInOwl#",
        },
    )


__all__ = ["ChEBIParser"]
