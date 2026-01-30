"""HMDB XML file parser for secondary-to-primary identifier mappings.

Uses lxml iterparse for memory-efficient streaming of large XML files.
"""

from __future__ import annotations

import zipfile
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path
from typing import BinaryIO

from lxml import etree
from tqdm import tqdm

from pysec2pri.models import (
    HMDBMapping,
    MappingCardinality,
    MappingSet,
)
from pysec2pri.parsers.base import BaseParser

# HMDB XML namespace
HMDB_NS = "http://www.hmdb.ca"
HMDB_NS_MAP = {"hmdb": HMDB_NS}

# Element tags we're interested in (with namespace)
TAG_METABOLITE = f"{{{HMDB_NS}}}metabolite"
TAG_ACCESSION = f"{{{HMDB_NS}}}accession"
TAG_SECONDARY_ACCESSIONS = f"{{{HMDB_NS}}}secondary_accessions"
TAG_NAME = f"{{{HMDB_NS}}}name"
TAG_SYNONYMS = f"{{{HMDB_NS}}}synonyms"
TAG_SYNONYM = f"{{{HMDB_NS}}}synonym"


def _get_text(elem: etree._Element | None) -> str | None:
    """Extract stripped text from an element, returning None if empty."""
    if elem is None or elem.text is None:
        return None
    text = elem.text.strip()
    return text or None


def _get_primary_accession(elem: etree._Element) -> str | None:
    """Get primary accession from a metabolite element."""
    return _get_text(elem.find(TAG_ACCESSION))


def _iter_secondary_accessions(elem: etree._Element) -> Iterable[str]:
    """Iterate over secondary accessions in a metabolite element."""
    container = elem.find(TAG_SECONDARY_ACCESSIONS)
    if container is None:
        return
    for child in container:
        text = _get_text(child)
        if text:
            yield text


def _get_primary_name(elem: etree._Element) -> str | None:
    """Get primary name from a metabolite element."""
    return _get_text(elem.find(TAG_NAME))


def _iter_synonyms(elem: etree._Element) -> Iterable[str]:
    """Iterate over synonyms in a metabolite element."""
    container = elem.find(TAG_SYNONYMS)
    if container is None:
        return
    for syn_elem in container.findall(TAG_SYNONYM):
        text = _get_text(syn_elem)
        if text:
            yield text


class HMDBParser(BaseParser):
    """Parser for HMDB XML files using lxml iterparse.

    Extracts secondary-to-primary HMDB identifier mappings and
    name-to-synonym relationships from HMDB XML format files.

    Uses lxml's iterparse for memory-efficient streaming of large XML
    files without loading the entire file into memory. Each metabolite
    element is processed and then cleared from memory.
    """

    datasource_name = "HMDB"
    default_source_url = "https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip"

    def parse(self, input_path: Path | str | None) -> MappingSet:
        """Parse HMDB XML files.

        Args:
            input_path: Path to a ZIP file containing hmdb_metabolites.xml,
                        or path to the XML file directly.

        Returns:
            MappingSet with HMDB identifier and synonym mappings.
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        input_path = Path(str(input_path))

        raw_id_mappings: list[tuple[str, str]] = []
        raw_name_mappings: list[tuple[str, str, str]] = []

        if input_path.suffix == ".zip":
            self._parse_zip(input_path, raw_id_mappings, raw_name_mappings)
        else:
            # Direct XML file
            metabolite_count = None
            if self.show_progress:
                metabolite_count = self._count_metabolites(input_path)

            with input_path.open("rb") as f:
                self._parse_xml_stream(f, raw_id_mappings, raw_name_mappings, metabolite_count)

        # Compute cardinality
        cardinality_map = self._compute_cardinality(raw_id_mappings)

        # Create mapping objects
        mappings: list[HMDBMapping] = []

        # Add ID mappings with progress
        id_iter = raw_id_mappings
        if self.show_progress:
            id_iter = tqdm(raw_id_mappings, desc="Creating ID mappings")

        for subject_id, object_id in id_iter:
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
            mapping = HMDBMapping(
                subject_id=subject_id,
                predicate_id=predicate,
                object_id=object_id,
                mapping_cardinality=cardinality,
                source_url=self.default_source_url,
            )
            mappings.append(mapping)

        # Add name/synonym mappings with progress
        name_iter = raw_name_mappings
        if self.show_progress:
            name_iter = tqdm(raw_name_mappings, desc="Creating synonym mappings")

        for subject_id, name, synonym in name_iter:
            mapping = HMDBMapping(
                subject_id=subject_id,
                predicate_id="oboInOwl:hasRelatedSynonym",
                object_id=synonym,
                subject_label=name,
                object_label=synonym,
                source_url=self.default_source_url,
            )
            mappings.append(mapping)

        mapping_set = MappingSet(
            mappings=mappings,
            datasource_name=self.datasource_name,
            version=self.version,
            mapping_set_id="omicsfixid_hmdb_01",
            mapping_set_description=(
                "Secondary to primary ID mappings for HMDB database, "
                "generated for the omicsFixID project."
            ),
            comment=(
                "object_label represents a synonym name by which the subject_label is also known."
            ),
            curie_map={
                "HMDB:": "http://www.hmdb.ca/metabolites/",
                "IAO:": "http://purl.obolibrary.org/obo/IAO_",
                "oboInOwl:": "http://www.geneontology.org/formats/oboInOwl#",
            },
        )

        return mapping_set

    def _parse_zip(
        self,
        zip_path: Path,
        id_mappings: list[tuple[str, str]],
        name_mappings: list[tuple[str, str, str]],
    ) -> None:
        """Parse HMDB from a ZIP archive containing hmdb_metabolites.xml."""
        with zipfile.ZipFile(zip_path, "r") as zf:
            # Find the main metabolites XML file
            xml_name = None
            for name in zf.namelist():
                if name.endswith("hmdb_metabolites.xml"):
                    xml_name = name
                    break

            if xml_name is None:
                msg = f"No hmdb_metabolites.xml found in {zip_path}"
                raise ValueError(msg)

            # Count metabolites for progress bar
            metabolite_count = None
            if self.show_progress:
                with zf.open(xml_name) as f:
                    import typing

                    metabolite_count = self._count_metabolites_stream(typing.cast("BinaryIO", f))

            # Parse the XML stream
            with zf.open(xml_name) as f:
                import typing

                self._parse_xml_stream(
                    typing.cast("BinaryIO", f), id_mappings, name_mappings, metabolite_count
                )

    def _count_metabolites(self, file_path: Path) -> int:
        """Count metabolite elements in XML file."""
        with file_path.open("rb") as f:
            return self._count_metabolites_stream(f)

    def _count_metabolites_stream(self, stream: BinaryIO) -> int:
        """Count metabolite elements by fast scanning for closing tags."""
        count = 0
        # Scan for </metabolite> tags
        for line in stream:
            if b"</metabolite>" in line:
                count += 1
        return count

    def _parse_xml_stream(
        self,
        stream: BinaryIO,
        id_mappings: list[tuple[str, str]],
        name_mappings: list[tuple[str, str, str]],
        total: int | None = None,
    ) -> None:
        """Parse XML stream using iterparse for memory efficiency.

        Uses lxml's iterparse to process one metabolite at a time,
        clearing processed elements from memory immediately.
        """
        # Create progress bar if enabled
        pbar = None
        if self.show_progress:
            pbar = tqdm(
                total=total,
                desc="Parsing HMDB metabolites",  # unit=" mol"
            )

        # Use iterparse with events for metabolite elements
        context = etree.iterparse(
            stream,
            events=("end",),
            tag=TAG_METABOLITE,
            recover=True,  # Continue on malformed XML
        )

        for _event, elem in context:
            self._process_metabolite(elem, id_mappings, name_mappings)

            # Clear element and its ancestors to free memory
            # This is critical for processing large files
            elem.clear()
            # Also clear preceding siblings that are no longer needed
            while elem.getprevious() is not None:
                del elem.getparent()[0]

            if pbar:
                pbar.update(1)

        if pbar:
            pbar.close()

    def _process_metabolite(
        self,
        elem: etree._Element,
        id_mappings: list[tuple[str, str]],
        name_mappings: list[tuple[str, str, str]],
    ) -> None:
        """Extract mappings from a single metabolite element."""
        subject_id = _get_primary_accession(elem)
        if not subject_id:
            return

        for object_id in _iter_secondary_accessions(elem):
            id_mappings.append((subject_id, object_id))

        primary_name = _get_primary_name(elem)
        if not primary_name:
            return

        for synonym in _iter_synonyms(elem):
            name_mappings.append((subject_id, primary_name, synonym))

    def _compute_cardinality(
        self,
        mappings: list[tuple[str, str]],
    ) -> dict[tuple[str, str], MappingCardinality]:
        """Compute cardinality for each mapping pair."""
        # Count how many secondaries each primary maps to
        primary_counts: dict[str, int] = defaultdict(int)
        # Count how many primaries each secondary maps to
        secondary_counts: dict[str, int] = defaultdict(int)

        for primary, secondary in mappings:
            primary_counts[primary] += 1
            secondary_counts[secondary] += 1

        # Assign cardinality to each pair
        result: dict[tuple[str, str], MappingCardinality] = {}
        for primary, secondary in mappings:
            p_count = primary_counts[primary]
            s_count = secondary_counts[secondary]

            if p_count == 1 and s_count == 1:
                cardinality = MappingCardinality.ONE_TO_ONE
            elif p_count == 1:
                cardinality = MappingCardinality.MANY_TO_ONE
            elif s_count == 1:
                cardinality = MappingCardinality.ONE_TO_MANY
            else:
                cardinality = MappingCardinality.MANY_TO_MANY

            result[(primary, secondary)] = cardinality

        return result


__all__ = ["HMDBParser"]
