"""HMDB XML file parser for secondary-to-primary identifier mappings.

This parser extracts ID-to-ID mappings from:
- hmdb_metabolites.xml  -> HMDB0... accessions  (HMDBMetaboliteParser)
- hmdb_proteins.xml     -> HMDBP... accessions  (HMDBProteinParser)

Each subclass reads its own YAML config:
- hmdb_metabolites.yaml
- hmdb_proteins.yaml

Uses SSSOM-compliant IdMappingSet with cardinality computation.
"""

from __future__ import annotations

import gzip
import io
import re
import zipfile
from collections.abc import Callable
from pathlib import Path
from typing import IO, TYPE_CHECKING

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

_BARE_NUM_RE = re.compile(r"^\d+$")


class HMDBParser(BaseParser):
    """Shared XML-parser for HMDB metabolite and protein files.

    Use :class:`HMDBMetaboliteParser` or :class:`HMDBProteinParser`.
    """

    # Subclasses must declare their own datasource_name so each reads
    # its own YAML config file (hmdb_metabolites.yaml / hmdb_proteins.yaml).
    datasource_name: str  # must be set by subclass

    # Internal helpers

    def _extract_version_from_file(self, file_path: Path) -> str | None:
        """Extract the HMDB dataset version from the first <version> element."""
        import re as _re

        version_re = _re.compile(r"<version>\s*([^<]+)\s*</version>")
        read_bytes = 2048

        try:
            if file_path.suffix == ".zip":
                with zipfile.ZipFile(file_path, "r") as zf:
                    xml_files = [n for n in zf.namelist() if n.endswith(".xml")]
                    if not xml_files:
                        return None
                    with zf.open(xml_files[0]) as fh:
                        head = fh.read(read_bytes).decode("utf-8", errors="replace")
            elif file_path.suffix == ".gz":
                with gzip.open(file_path, "rt", encoding="utf-8") as fh:
                    head = fh.read(read_bytes)
            else:
                with file_path.open("r", encoding="utf-8") as fh:
                    head = fh.read(read_bytes)
        except Exception as e:
            logger.debug("Could not read HMDB version from %s: %s", file_path, e)
            return None

        m = version_re.search(head)
        return m.group(1).strip() if m else None

    def _parse_xml(
        self,
        file_path: Path,
        element_tag: str,
        prefix: str,
        desc: str,
        secondary_normaliser: Callable[[str], str] | None = None,
    ) -> Sec2PriMappingSet:
        """Parse an HMDB XML file and return a mapping set.

        Args:
            file_path: Path to the XML file (plain, ``.zip``, or ``.gz``).
            element_tag: Top-level record element name (``"metabolite"`` or
                ``"protein"``).
            prefix: CURIE prefix (``"HMDB"`` for both record types).
            desc: Progress-bar description string.
            secondary_normaliser: Optional callable ``(raw_str) -> str``
                applied to each raw secondary accession before prefixing.

        Returns:
            :class:`Sec2PriMappingSet` with computed cardinalities and
            ``_primary_ids`` populated.
        """
        xml_content = self._read_xml_content(file_path)
        if xml_content is None:
            return self._create_mapping_set([])

        self._resolve_version(file_path)
        m_meta = self.get_mapping_metadata()
        fixed: dict[str, str | None] = {
            "predicate_id": m_meta["predicate_id"],
            "predicate_label": m_meta.get("predicate_label"),
            "mapping_justification": m_meta["mapping_justification"],
            "subject_source": m_meta.get("subject_source"),
            "object_source": m_meta.get("object_source"),
            "mapping_tool": m_meta.get("mapping_tool"),
            "license": m_meta.get("license"),
        }
        rows_data: list[dict[str, str]] = []
        primary_ids_found: set[str] = set()

        try:
            context = DefusedET.iterparse(xml_content, events=("end",))
            for _event, elem in self._progress(context, desc=desc):
                tag = elem.tag.replace(f"{{{HMDB_NS['hmdb']}}}", "")
                if tag == element_tag:
                    accession_elem = elem.find("hmdb:accession", HMDB_NS)
                    if accession_elem is None:
                        accession_elem = elem.find("accession")
                    if accession_elem is not None and accession_elem.text:
                        primary_ids_found.add(f"{prefix}:{accession_elem.text.strip()}")
                    rows_data.extend(self._process_record(elem, prefix, secondary_normaliser))
                    elem.clear()
        except DefusedET.ParseError:
            rows_data = self._parse_simple_xml(file_path, element_tag, prefix, secondary_normaliser)

        mappings = self._build_mappings(
            rows_data, fixed, desc="Building HMDB mappings", total=len(rows_data)
        )
        ms = self._create_mapping_set(mappings)
        if primary_ids_found:
            object.__setattr__(ms, "_primary_ids", primary_ids_found)
        return ms

    def _process_record(
        self,
        elem: Element,
        prefix: str,
        secondary_normaliser: Callable[[str], str] | None,
    ) -> list[dict[str, str]]:
        """Extract row dicts from a single record element."""
        # Primary accession
        accession_elem = elem.find("hmdb:accession", HMDB_NS)
        if accession_elem is None:
            accession_elem = elem.find("accession")
        if accession_elem is None or not accession_elem.text:
            return []
        primary_raw = accession_elem.text.strip()
        primary_id = f"{prefix}:{primary_raw}"

        # Label name
        name_elem = elem.find("hmdb:name", HMDB_NS)
        if name_elem is None:
            name_elem = elem.find("name")
        primary_label = name_elem.text.strip() if name_elem is not None and name_elem.text else ""

        # Secondary accessions
        sec_block = elem.find("hmdb:secondary_accessions", HMDB_NS)
        if sec_block is None:
            sec_block = elem.find("secondary_accessions")
        if sec_block is None:
            return []

        rows: list[dict[str, str]] = []
        for sec_elem in sec_block:
            if not sec_elem.text:
                continue
            raw_sec = sec_elem.text.strip()
            if secondary_normaliser is not None:
                raw_sec = secondary_normaliser(raw_sec)
            secondary_id = f"{prefix}:{raw_sec}"
            rows.append(
                {
                    "subject_id": secondary_id,
                    "subject_label": "",
                    "object_id": primary_id,
                    "object_label": primary_label,
                }
            )
        return rows

    def _parse_simple_xml(
        self,
        file_path: Path,
        element_tag: str,
        prefix: str,
        secondary_normaliser: Callable[[str], str] | None,
    ) -> list[dict[str, str]]:
        """Fallback parser for test files without a namespace declaration."""
        rows_data: list[dict[str, str]] = []
        try:
            tree = DefusedET.parse(file_path)
            root = tree.getroot()
            for record in root.findall(f".//{element_tag}"):
                rows_data.extend(self._process_record(record, prefix, secondary_normaliser))
        except Exception as e:
            logger.warning("Failed to parse HMDB XML: %s", e)
        return rows_data

    def _read_xml_content(self, file_path: Path) -> IO[bytes] | IO[str] | None:
        """Open an XML file, transparently handling ``.zip`` and ``.gz`` wrappers.

        Args:
            file_path: Path to the XML file (plain, ``.zip``, or ``.gz``).

        Returns:
            An open file-like object ready for parsing, or ``None`` on failure.
        """
        if file_path.suffix == ".zip":
            try:
                with zipfile.ZipFile(file_path, "r") as zf:
                    xml_files = [n for n in zf.namelist() if n.endswith(".xml")]
                    if xml_files:
                        return io.BytesIO(zf.read(xml_files[0]))
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
        """Delegate to the base-class factory."""
        return self.create_mapping_set(mappings, mapping_type)

    def parse_primary_ids(
        self,
        metabolites_path: Path | str | None = None,
        proteins_path: Path | str | None = None,
    ) -> Sec2PriMappingSet:
        """Return a mapping set containing the full list of current HMDB primary IDs.

        Reads one or both of ``hmdb_metabolites.xml`` and ``hmdb_proteins.xml``
        and collects all primary accession numbers.  The returned mapping set
        has an empty ``mappings`` list; ``_primary_ids`` is populated with every
        current ``HMDB:<acc>`` CURIE.

        Args:
            metabolites_path: Path to ``hmdb_metabolites.xml`` (or zip/gz).
            proteins_path: Path to ``hmdb_proteins.xml`` (or zip/gz).

        Returns:
            :class:`~pysec2pri.parsers.base.IdMappingSet` with ``_primary_ids``
            populated.  At least one of the two path arguments must be supplied.
        """
        if metabolites_path is None and proteins_path is None:
            raise ValueError("At least one of metabolites or proteins_path must be supplied.")

        primary_ids: set[str] = set()

        if metabolites_path is not None:
            ms_m = self.parse(metabolites_path)
            primary_ids |= object.__getattribute__(ms_m, "_primary_ids")

        if proteins_path is not None:
            ms_p = self.parse(proteins_path)
            primary_ids |= object.__getattribute__(ms_p, "_primary_ids")

        ms = self._create_mapping_set([], mapping_type="id")
        object.__setattr__(ms, "_primary_ids", primary_ids)
        return ms


# Concrete parsers


class HMDBMetaboliteParser(HMDBParser):
    """Parser for ``hmdb_metabolites.xml``.

    Reads configuration from ``hmdb_metabolites.yaml``.
    Primary accessions have the form ``HMDB0…`` (e.g. ``HMDB0000001``).
    """

    datasource_name = "hmdb_metabolites"

    @property
    def source_url(self) -> str:
        """Metabolites download URL from ``hmdb_metabolites.yaml``."""
        return self.get_download_url("metabolites") or ""

    def parse(self, input_path: Path | str | None) -> Sec2PriMappingSet:
        """Parse ``hmdb_metabolites.xml`` (or ``.zip`` / ``.gz``).

        Args:
            input_path: Path to the metabolites XML file.

        Returns:
            :class:`Sec2PriMappingSet` for metabolite accessions.
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        return self._parse_xml(
            Path(input_path),
            element_tag="metabolite",
            prefix="HMDB",
            desc="Parsing HMDB metabolites XML",
        )


class HMDBProteinParser(HMDBParser):
    """Parser for ``hmdb_proteins.xml``.

    Reads configuration from ``hmdb_proteins.yaml``.
    Primary accessions have the form ``HMDBP…`` (e.g. ``HMDBP00001``).
    Secondary accessions may be bare numbers (legacy format, e.g. ``5229``)
    or full ``HMDBP`` accessions; both are normalised by
    :meth:`_normalise_protein_secondary`.
    """

    datasource_name = "hmdb_proteins"

    @property
    def source_url(self) -> str:
        """Proteins download URL from ``hmdb_proteins.yaml``."""
        return self.get_download_url("proteins") or ""

    @staticmethod
    def _normalise_protein_secondary(raw: str) -> str:
        """Normalise a raw secondary protein accession to a ``HMDBP``-prefixed form.

        Bare numeric values (legacy IDs) are zero-padded to five digits and
        prefixed with ``HMDBP`` to match the primary accession pattern.
        Already-prefixed values are returned unchanged.

        Args:
            raw: Raw secondary accession string (e.g. ``"5229"`` or
                ``"HMDBP05229"``).

        Returns:
            Normalised accession string (e.g. ``"HMDBP05229"``).
        """
        raw = raw.strip()
        if _BARE_NUM_RE.match(raw):
            return f"HMDBP{int(raw):05d}"
        return raw

    def parse(self, input_path: Path | str | None) -> Sec2PriMappingSet:
        """Parse ``hmdb_proteins.xml`` (or ``.zip`` / ``.gz``).

        Args:
            input_path: Path to the proteins XML file.

        Returns:
            :class:`Sec2PriMappingSet` for protein accessions.
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        return self._parse_xml(
            Path(input_path),
            element_tag="protein",
            prefix="HMDB",
            desc="Parsing HMDB proteins XML",
            secondary_normaliser=self._normalise_protein_secondary,
        )


__all__ = [
    "HMDBMetaboliteParser",
    "HMDBParser",
    "HMDBProteinParser",
]
