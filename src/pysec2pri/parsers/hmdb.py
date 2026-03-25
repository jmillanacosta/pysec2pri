"""HMDB XML file parser for secondary-to-primary identifier mappings.

This parser extracts ID-to-ID mappings from:
- hmdb_metabolites.xml  → HMDB0... accessions
- hmdb_proteins.xml     → HMDBP... accessions

Uses SSSOM-compliant IdMappingSet with cardinality computation.
"""

from __future__ import annotations

import gzip
import re
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

# Pattern for bare numeric legacy protein IDs (e.g. "5229")
_BARE_NUM_RE = re.compile(r"^\d+$")


class HMDBParser(BaseParser):
    """Parser for HMDB XML files (metabolites and proteins).

    Extracts secondary-to-primary HMDB accession mappings from
    ``hmdb_metabolites.xml`` and/or ``hmdb_proteins.xml``.
    """

    datasource_name = "hmdb"

    @property
    def source_url(self) -> str:
        """Get the metabolites download URL from config."""
        return self.get_download_url("metabolites") or ""

    def parse(self, input_path: Path | str | None) -> Sec2PriMappingSet:
        """Parse HMDB metabolites XML file.

        Args:
            input_path: Path to ``hmdb_metabolites.xml`` (or ``.zip``/``.gz``).

        Returns:
            IdMappingSet for metabolite accessions.
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        return self._parse_xml(
            Path(input_path),
            element_tag="metabolite",
            prefix="HMDB",
            desc="Parsing HMDB metabolites XML",
        )

    def parse_proteins(self, input_path: Path | str) -> Sec2PriMappingSet:
        """Parse HMDB proteins XML file.

        Primary accessions have the form ``HMDBP00001``.
        Secondary accessions may be bare numbers (legacy format, e.g.
        ``5229``) or full ``HMDBP`` accessions; both are normalised to
        ``HMDBP:HMDBP<zero-padded>`` using the same prefix logic as
        the metabolites parser.

        Args:
            input_path: Path to ``hmdb_proteins.xml`` (or ``.zip``/``.gz``).

        Returns:
            IdMappingSet for protein accessions.
        """
        return self._parse_xml(
            Path(input_path),
            element_tag="protein",
            prefix="HMDB",
            desc="Parsing HMDB proteins XML",
            secondary_normaliser=self._normalise_protein_secondary,
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _normalise_protein_secondary(raw: str) -> str:
        """Return a CURIE-ready protein accession from a raw secondary value.

        Bare numeric values (legacy IDs) are zero-padded to 5 digits and
        prefixed with ``HMDBP`` so they match the primary ``HMDBP`` pattern.
        Already-prefixed values are returned unchanged.
        """
        raw = raw.strip()
        if _BARE_NUM_RE.match(raw):
            return f"HMDBP{int(raw):05d}"
        return raw

    def _parse_xml(
        self,
        file_path: Path,
        element_tag: str,
        prefix: str,
        desc: str,
        secondary_normaliser: Any = None,
    ) -> Sec2PriMappingSet:
        """Generic HMDB XML parser for metabolites and proteins.

        Args:
            file_path: Path to the XML file (may be ``.zip`` or ``.gz``).
            element_tag: Top-level record element name (``"metabolite"`` or
                ``"protein"``).
            prefix: CURIE prefix (``"HMDB"`` for both).
            desc: Progress-bar description string.
            secondary_normaliser: Optional callable ``(raw_str) -> str``
                applied to each raw secondary accession before prefixing.

        Returns:
            Sec2PriMappingSet with computed cardinalities.
        """
        xml_content = self._read_xml_content(file_path)
        if xml_content is None:
            return self._create_mapping_set([])

        m_meta = self.get_mapping_metadata()
        mappings: list[Mapping] = []

        try:
            context = DefusedET.iterparse(xml_content, events=("end",))
            for _event, elem in self._progress(context, desc=desc):
                tag = elem.tag.replace(f"{{{HMDB_NS['hmdb']}}}", "")
                if tag == element_tag:
                    mappings.extend(
                        self._process_record(
                            elem, prefix, m_meta, secondary_normaliser
                        )
                    )
                    elem.clear()
        except DefusedET.ParseError:
            mappings = self._parse_simple_xml(
                file_path, element_tag, prefix, m_meta, secondary_normaliser
            )

        return self._create_mapping_set(mappings)

    def _process_record(
        self,
        elem: Element,
        prefix: str,
        m_meta: dict[str, Any],
        secondary_normaliser: Any,
    ) -> list[Mapping]:
        """Extract mappings from a single record element."""
        # Primary accession
        accession_elem = elem.find("hmdb:accession", HMDB_NS)
        if accession_elem is None:
            accession_elem = elem.find("accession")
        if accession_elem is None or not accession_elem.text:
            return []
        primary_raw = accession_elem.text.strip()
        primary_id = f"{prefix}:{primary_raw}"

        # Label (name)
        name_elem = elem.find("hmdb:name", HMDB_NS)
        if name_elem is None:
            name_elem = elem.find("name")
        primary_label = (
            name_elem.text.strip()
            if name_elem is not None and name_elem.text
            else ""
        )

        # Secondary accessions
        sec_block = elem.find("hmdb:secondary_accessions", HMDB_NS)
        if sec_block is None:
            sec_block = elem.find("secondary_accessions")
        if sec_block is None:
            return []

        mappings: list[Mapping] = []
        for sec_elem in sec_block:
            if not sec_elem.text:
                continue
            raw_sec = sec_elem.text.strip()
            if secondary_normaliser is not None:
                raw_sec = secondary_normaliser(raw_sec)
            secondary_id = f"{prefix}:{raw_sec}"
            mappings.append(
                Mapping(
                    subject_id=primary_id,
                    subject_label=primary_label,
                    object_id=secondary_id,
                    predicate_id=m_meta["predicate_id"],
                    predicate_label=m_meta.get("predicate_label"),
                    mapping_justification=m_meta["mapping_justification"],
                    subject_source=m_meta.get("subject_source"),
                    object_source=m_meta.get("object_source"),
                    mapping_tool=m_meta.get("mapping_tool"),
                    license=m_meta.get("license"),
                )
            )
        return mappings

    def _parse_simple_xml(
        self,
        file_path: Path,
        element_tag: str,
        prefix: str,
        m_meta: dict[str, Any],
        secondary_normaliser: Any,
    ) -> list[Mapping]:
        """Fallback parser for test files without namespace."""
        mappings: list[Mapping] = []
        try:
            tree = DefusedET.parse(file_path)
            root = tree.getroot()
            for record in root.findall(f".//{element_tag}"):
                mappings.extend(
                    self._process_record(
                        record, prefix, m_meta, secondary_normaliser
                    )
                )
        except Exception as e:
            logger.warning("Failed to parse HMDB XML: %s", e)
        return mappings

    def _read_xml_content(self, file_path: Path) -> Any:
        """Read XML content, handling compressed files."""
        if file_path.suffix == ".zip":
            try:
                with zipfile.ZipFile(file_path, "r") as zf:
                    xml_files = [
                        n for n in zf.namelist() if n.endswith(".xml")
                    ]
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
