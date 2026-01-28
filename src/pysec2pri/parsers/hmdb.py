"""HMDB XML file parser for secondary-to-primary identifier mappings.

Uses lxml iterparse for memory-efficient streaming of large XML files.
"""

from __future__ import annotations

import zipfile
from collections import defaultdict
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


class HMDBParser(BaseParser):
    """Parser for HMDB XML files using lxml iterparse.

    Extracts secondary-to-primary HMDB identifier mappings and
    name-to-synonym relationships from HMDB XML format files.

    Uses lxml's iterparse for memory-efficient streaming of large XML
    files without loading the entire file into memory. Each metabolite
    element is processed and then cleared from memory.
    """

    datasource_name = "HMDB"
    default_source_url = (
        "https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip"
    )

    def parse(self, input_path: Path | str) -> MappingSet:
        """Parse HMDB XML files.

        Args:
            input_path: Path to a ZIP file containing hmdb_metabolites.xml,
                        or path to the XML file directly.

        Returns:
            MappingSet with HMDB identifier and synonym mappings.
        """
        input_path = Path(input_path)

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
                self._parse_xml_stream(
                    f, raw_id_mappings, raw_name_mappings, metabolite_count
                )

        # Compute cardinality
        cardinality_map = self._compute_cardinality(raw_id_mappings)

        # Create mapping objects
        mappings: list[HMDBMapping] = []

        # Add ID mappings with progress
        id_iter = raw_id_mappings
        if self.show_progress:
            id_iter = tqdm(raw_id_mappings, desc="Creating ID mappings")

        for primary_id, secondary_id in id_iter:
            cardinality = cardinality_map.get(
                (primary_id, secondary_id),
                MappingCardinality.ONE_TO_ONE,
            )
            mapping = HMDBMapping(
                primary_id=primary_id,
                secondary_id=secondary_id,
                mapping_cardinality=cardinality,
                source_url=self.default_source_url,
            )
            mappings.append(mapping)

        # Add name/synonym mappings with progress
        name_iter = raw_name_mappings
        if self.show_progress:
            name_iter = tqdm(
                raw_name_mappings, desc="Creating synonym mappings"
            )

        for primary_id, name, synonym in name_iter:
            mapping = HMDBMapping(
                primary_id=primary_id,
                secondary_id=None,
                primary_label=name,
                secondary_label=synonym,
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
                "object_label represents a synonym name by which "
                "the subject_label is also known."
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
                    metabolite_count = self._count_metabolites_stream(f)

            # Parse the XML stream
            with zf.open(xml_name) as f:
                self._parse_xml_stream(
                    f, id_mappings, name_mappings, metabolite_count
                )

    def _count_metabolites(self, file_path: Path) -> int:
        """Count metabolite elements in XML file."""
        with file_path.open("rb") as f:
            return self._count_metabolites_stream(f)

    def _count_metabolites_stream(self, stream: BinaryIO) -> int:
        """Count metabolite elements by fast scanning for closing tags."""
        count = 0
        # Fast byte-level scan for </metabolite> tags
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
                total=total, desc="Parsing HMDB metabolites",  # unit=" mol"
            )

        # Use iterparse with events for metabolite elements
        context = etree.iterparse(
            stream,
            events=("end",),
            tag=TAG_METABOLITE,
            recover=True,  # Continue on malformed XML
        )

        for event, elem in context:
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
        # Get primary accession (first direct child accession)
        accession_elem = elem.find(TAG_ACCESSION)
        if accession_elem is None or not accession_elem.text:
            return

        primary_id = accession_elem.text.strip()
        if not primary_id:
            return

        # Get secondary accessions
        sec_container = elem.find(TAG_SECONDARY_ACCESSIONS)
        if sec_container is not None:
            for sec_elem in sec_container:
                if sec_elem.text:
                    secondary_id = sec_elem.text.strip()
                    if secondary_id:
                        id_mappings.append((primary_id, secondary_id))

        # Get name
        name_elem = elem.find(TAG_NAME)
        primary_name = None
        if name_elem is not None and name_elem.text:
            primary_name = name_elem.text.strip()

        # Get synonyms
        if primary_name:
            syn_container = elem.find(TAG_SYNONYMS)
            if syn_container is not None:
                for syn_elem in syn_container.findall(TAG_SYNONYM):
                    if syn_elem.text:
                        synonym = syn_elem.text.strip()
                        if synonym:
                            name_mappings.append(
                                (primary_id, primary_name, synonym)
                            )

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
