"""HMDB XML file parser for secondary-to-primary identifier mappings.

This parser extracts ID-to-ID mappings:
- Secondary accessions -> primary accessions (from hmdb_metabolites.xml)

Uses SSSOM-compliant IdMappingSet with cardinality computation.
"""

from __future__ import annotations

import gzip
import zipfile
from pathlib import Path
from typing import TYPE_CHECKING, Any

import defusedxml.ElementTree as DefusedET
from sssom_schema import Mapping

from pysec2pri.logging import logger
from pysec2pri.parsers.base import (
    BaseParser,
    Sec2PriMappingSet,
)

if TYPE_CHECKING:
    from xml.etree.ElementTree import Element

# HMDB XML namespace
HMDB_NS = {"hmdb": "http://www.hmdb.ca"}


class HMDBParser(BaseParser):
    """Parser for HMDB XML files.

    Extracts secondary-to-primary HMDB accession mappings from
    the hmdb_metabolites.xml file.

    Returns IdMappingSet for all mappings (HMDB only has ID mappings).
    """

    datasource_name = "hmdb"

    @property
    def source_url(self) -> str:
        """Get the metabolites download URL from config."""
        return self.get_download_url("metabolites") or ""

    def parse(self, input_path: Path | str | None) -> Sec2PriMappingSet:
        """Parse HMDB metabolites XML file into an IdMappingSet.

        Args:
            input_path: Path to hmdb_metabolites.xml (or .zip/.gz).

        Returns:
            IdMappingSet with computed cardinalities based on IDs.
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        input_path = Path(input_path)

        # Parse XML file for ID mappings
        mappings = self._parse_metabolites_xml(input_path)

        # Create IdMappingSet and compute cardinalities
        mapping_set = self._create_mapping_set(mappings)
        return mapping_set

    def _parse_metabolites_xml(self, file_path: Path) -> list[Mapping]:
        """Parse HMDB metabolites XML for secondary -> primary mappings.

        Args:
            file_path: Path to the XML file (can be .zip or .gz).

        Returns:
            List of SSSOM Mapping objects.
        """
        # Handle compressed files
        xml_content = self._read_xml_content(file_path)
        if xml_content is None:
            return []

        m_meta = self.get_mapping_metadata()
        mappings: list[Mapping] = []

        # Parse XML iteratively to handle large files
        try:
            context = DefusedET.iterparse(xml_content, events=("end",))
            for _event, elem in self._progress(context, desc="Parsing HMDB XML"):
                # Look for metabolite elements
                tag = elem.tag.replace(f"{{{HMDB_NS['hmdb']}}}", "")
                if tag == "metabolite":
                    metabolite_mappings = self._process_metabolite(elem, m_meta)
                    mappings.extend(metabolite_mappings)
                    elem.clear()  # Free memory
        except DefusedET.ParseError:
            # Fall back to simpler parsing for test files
            mappings = self._parse_simple_xml(file_path, m_meta)

        return mappings

    def _process_metabolite(self, elem: Element, m_meta: dict[str, Any]) -> list[Mapping]:
        """Process a single metabolite element.

        Args:
            elem: The metabolite XML element.
            m_meta: Mapping metadata from config.

        Returns:
            List of mappings for this metabolite.
        """
        mappings: list[Mapping] = []

        # Get primary accession
        accession_elem = elem.find("hmdb:accession", HMDB_NS)
        if accession_elem is None:
            accession_elem = elem.find("accession")
        if accession_elem is None or not accession_elem.text:
            return []

        primary_id = accession_elem.text.strip()

        # Get name for label
        name_elem = elem.find("hmdb:name", HMDB_NS)
        if name_elem is None:
            name_elem = elem.find("name")
        primary_label = name_elem.text.strip() if name_elem is not None and name_elem.text else ""

        # Get secondary accessions
        sec_accessions = elem.find("hmdb:secondary_accessions", HMDB_NS)
        if sec_accessions is None:
            sec_accessions = elem.find("secondary_accessions")

        if sec_accessions is not None:
            for sec_elem in sec_accessions:
                if sec_elem.text:
                    secondary_id = sec_elem.text.strip()
                    # SSSOM: subject = primary (current), object = secondary (old)
                    mapping = Mapping(
                        subject_id=f"HMDB:{primary_id}",
                        subject_label=primary_label,
                        object_id=f"HMDB:{secondary_id}",
                        predicate_id=m_meta["predicate_id"],
                        mapping_justification=m_meta["mapping_justification"],
                        subject_source=m_meta.get("subject_source"),
                        object_source=m_meta.get("object_source"),
                        mapping_tool=m_meta.get("mapping_tool"),
                        license=m_meta.get("license"),
                    )
                    mappings.append(mapping)

        return mappings

    def _parse_simple_xml(self, file_path: Path, m_meta: dict[str, Any]) -> list[Mapping]:
        """Parse XML for test files without namespace.

        Args:
            file_path: Path to the XML file.
            m_meta: Mapping metadata from config.

        Returns:
            List of SSSOM Mapping objects.
        """
        mappings: list[Mapping] = []

        try:
            tree = DefusedET.parse(file_path)
            root = tree.getroot()

            for metabolite in root.findall(".//metabolite"):
                accession = metabolite.find("accession")
                name = metabolite.find("name")
                sec_accessions = metabolite.find("secondary_accessions")

                if accession is None or not accession.text:
                    continue

                primary_id = accession.text.strip()
                primary_label = name.text.strip() if name is not None and name.text else ""

                if sec_accessions is not None:
                    for sec_elem in sec_accessions:
                        if sec_elem.text:
                            secondary_id = sec_elem.text.strip()
                            mapping = Mapping(
                                subject_id=f"HMDB:{primary_id}",
                                subject_label=primary_label,
                                object_id=f"HMDB:{secondary_id}",
                                predicate_id=m_meta["predicate_id"],
                                mapping_justification=m_meta["mapping_justification"],
                                subject_source=m_meta.get("subject_source"),
                                object_source=m_meta.get("object_source"),
                                mapping_tool=m_meta.get("mapping_tool"),
                                license=m_meta.get("license"),
                            )
                            mappings.append(mapping)
        except Exception as e:
            logger.warning("Failed to parse HMDB XML: %s", e)

        return mappings

    def _read_xml_content(self, file_path: Path) -> Any:
        """Read XML content, handling compressed files.

        Args:
            file_path: Path to the file.

        Returns:
            File-like object for parsing, or None if failed.
        """
        if file_path.suffix == ".zip":
            try:
                with zipfile.ZipFile(file_path, "r") as zf:
                    xml_files = [n for n in zf.namelist() if n.endswith(".xml")]
                    if xml_files:
                        return zf.open(xml_files[0])
            except Exception as e:
                logger.warning("Failed to read zip file %s: %s", file_path, e)
                return None
        elif file_path.suffix == ".gz":
            try:
                return gzip.open(file_path, "rt", encoding="utf-8")
            except Exception as e:
                logger.warning("Failed to read gzip file %s: %s", file_path, e)
                return None
        else:
            try:
                return file_path.open("r", encoding="utf-8")
            except Exception as e:
                logger.warning("Failed to read file %s: %s", file_path, e)
                return None

    def _create_mapping_set(
        self, mappings: list[Mapping], mapping_type: str = "id"
    ) -> Sec2PriMappingSet:
        """Delegate to base class method."""
        return self.create_mapping_set(mappings, mapping_type)


__all__ = ["HMDBParser"]
