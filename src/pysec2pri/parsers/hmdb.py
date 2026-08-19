"""HMDB XML file parser for secondary-to-primary identifier mappings.

This parser extracts:
1. ID-to-ID mappings from:
   - hmdb_metabolites.xml  -> HMDB0... accessions  (HMDBMetaboliteParser)
   - hmdb_proteins.xml     -> HMDBP... accessions  (HMDBProteinParser)
2. Label-to-label mappings: synonyms -> primary names

Each subclass reads a YAML config:
- hmdb_metabolites.yaml
- hmdb_proteins.yaml

Uses SSSOM-compliant IdMappingSet and LabelMappingSet with cardinality
computation.
"""

from __future__ import annotations

import gzip
import io
import re
import zipfile
from pathlib import Path
from typing import IO, TYPE_CHECKING, ClassVar

import defusedxml.ElementTree as DefusedET

from pysec2pri.logging import logger
from pysec2pri.parsers.base import (
    BaseMappingSet,
    BaseParser,
)

if TYPE_CHECKING:
    from collections.abc import Callable
    from xml.etree.ElementTree import Element

# HMDB XML namespace
HMDB_NS = {"hmdb": "http://www.hmdb.ca"}

_BARE_NUM_RE = re.compile(r"^\d+$")


class HMDBParser(BaseParser):
    """Shared XML-parser for HMDB metabolite and protein files.

    Use :class:`HMDBMetaboliteParser` or :class:`HMDBProteinParser`.
    """

    # Subclasses declare their own datasource_name so each reads
    # their YAML config file (hmdb_metabolites.yaml / hmdb_proteins.yaml)
    datasource_name: str  # must be set by subclass
    element_tag: ClassVar[str]
    prefix: ClassVar[str]
    species: ClassVar[str] = "9606"

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

    @staticmethod
    def _accession_and_name(elem: Element, prefix: str) -> tuple[str | None, str]:
        """Return ``(primary_id, primary_name)`` for a record element."""
        accession_elem = elem.find("hmdb:accession", HMDB_NS)
        if accession_elem is None:
            accession_elem = elem.find("accession")
        if accession_elem is None or not accession_elem.text:
            return None, ""
        primary_id = f"{prefix}:{accession_elem.text.strip()}"

        name_elem = elem.find("hmdb:name", HMDB_NS)
        if name_elem is None:
            name_elem = elem.find("name")
        primary_name = name_elem.text.strip() if name_elem is not None and name_elem.text else ""
        return primary_id, primary_name

    def _iter_records(
        self,
        file_path: Path,
        element_tag: str,
        process_fn: Callable[[Element], list[dict[str, str]]],
        desc: str,
    ) -> list[dict[str, str]]:
        """Walk every ``element_tag`` record in an HMDB XML file, building rows.

        Falls back to a plain DOM parse for test files without a namespace
        declaration (``DefusedET.iterparse`` requires one).

        Args:
            file_path: Path to the XML file (plain, ``.zip``, or ``.gz``).
            element_tag: Top-level record element name (``"metabolite"`` or
                ``"protein"``).
            process_fn: Called with each record element; returns its row dicts.
            desc: Progress-bar description string.

        Returns:
            Flattened list of row dicts across every record.
        """
        xml_content = self._read_xml_content(file_path)
        if xml_content is None:
            return []

        rows_data: list[dict[str, str]] = []
        try:
            with xml_content:
                context = DefusedET.iterparse(xml_content, events=("end",))
                for _event, elem in self._progress(context, desc=desc):
                    tag = elem.tag.replace(f"{{{HMDB_NS['hmdb']}}}", "")
                    if tag == element_tag:
                        rows_data.extend(process_fn(elem))
                        elem.clear()
        except DefusedET.ParseError:
            try:
                root = DefusedET.parse(file_path).getroot()
                if root is not None:
                    for record in root.findall(f".//{element_tag}"):
                        rows_data.extend(process_fn(record))
            except Exception as e:
                logger.warning("Failed to parse HMDB XML: %s", e)
        return rows_data

    def _finalize_mapping_set(
        self,
        rows_data: list[dict[str, str]],
        fixed: dict[str, str | None],
        *,
        build_desc: str,
        mapping_type: str = "id",
        primary_ids: set[str] | None = None,
        primary_labels: dict[str, set[str]] | None = None,
    ) -> BaseMappingSet:
        """Build ``Mapping`` objects from row dicts and wrap them in a mapping set."""
        mappings = self._build_mappings(rows_data, fixed, desc=build_desc, total=len(rows_data))
        return self.create_mapping_set(
            mappings,
            mapping_type=mapping_type,
            primary_ids=primary_ids,
            primary_labels=primary_labels,
        )

    def _parse_xml(
        self,
        file_path: Path,
        element_tag: str,
        prefix: str,
        desc: str,
    ) -> BaseMappingSet:
        """Parse an HMDB XML file and return an ID mapping set.

        Args:
            file_path: Path to the XML file (plain, ``.zip``, or ``.gz``).
            element_tag: Top-level record element name (``"metabolite"`` or
                ``"protein"``).
            prefix: CURIE prefix (``"HMDB"`` for both record types).
            desc: Progress-bar description string.

        Returns:
            :class:`BaseMappingSet` with computed cardinalities and
            ``_primary_ids`` populated.
        """
        self._resolve_version(file_path)
        m_meta = self.get_mapping_metadata()
        fixed: dict[str, str | None] = {
            **self._fixed_mapping_fields(),
            "predicate_id": m_meta["predicate_id"],
            "predicate_label": m_meta.get("predicate_label"),
        }
        primary_ids_found: set[str] = set()
        rows_data = self._iter_records(
            file_path,
            element_tag,
            lambda elem: self._process_record(elem, prefix, primary_ids_found),
            desc,
        )
        return self._finalize_mapping_set(
            rows_data,
            fixed,
            build_desc="Building HMDB mappings",
            primary_ids=primary_ids_found or None,
        )

    def _process_record(
        self,
        elem: Element,
        prefix: str,
        primary_ids: set[str],
    ) -> list[dict[str, str]]:
        """Extract secondary-to-primary row dicts from a single record element.

        Bare-numeric ``secondary_accessions`` entries are skipped: in the
        proteins file these are HMDB's internal database row ids (not legacy
        accessions), and a bare integer is never a valid HMDB accession for
        either record type. See issue #44.
        """
        primary_id, primary_label = self._accession_and_name(elem, prefix)
        if primary_id is None:
            return []
        primary_ids.add(primary_id)

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
            if _BARE_NUM_RE.match(raw_sec):
                continue
            rows.append(
                {
                    "subject_id": f"{prefix}:{raw_sec}",
                    "subject_label": "",
                    "object_id": primary_id,
                    "object_label": primary_label,
                }
            )
        return rows

    def _parse_synonyms_xml(
        self,
        file_path: Path,
        element_tag: str,
        prefix: str,
        desc: str,
    ) -> BaseMappingSet:
        """Parse an HMDB XML file into a synonym-to-name label mapping set.

        Args:
            file_path: Path to the XML file (plain, ``.zip``, or ``.gz``).
            element_tag: Top-level record element name (``"metabolite"`` or
                ``"protein"``).
            prefix: CURIE prefix (``"HMDB"`` for both record types).
            desc: Progress-bar description string.

        Returns:
            :class:`BaseMappingSet` with computed cardinalities and
            ``_primary_labels`` populated.
        """
        self._resolve_version(file_path)
        fixed = self._fixed_mapping_fields()
        primary_labels_found: dict[str, set[str]] = {}
        rows_data = self._iter_records(
            file_path,
            element_tag,
            lambda elem: self._process_synonym_record(elem, prefix, primary_labels_found),
            desc,
        )
        return self._finalize_mapping_set(
            rows_data,
            fixed,
            build_desc="Building HMDB synonym mappings",
            mapping_type="label",
            primary_labels=primary_labels_found or None,
        )

    def _process_synonym_record(
        self,
        elem: Element,
        prefix: str,
        primary_labels: dict[str, set[str]],
    ) -> list[dict[str, str]]:
        """Extract synonym label rows from a single record element.

        Also records the record's primary name in ``primary_labels`` even
        when it carries no synonyms, so ``to_pri_labels()`` covers the full
        entity set.
        """
        primary_id, primary_label = self._accession_and_name(elem, prefix)
        if primary_id is None or not primary_label:
            return []
        primary_labels.setdefault(primary_label, set()).add(primary_id)

        syn_block = elem.find("hmdb:synonyms", HMDB_NS)
        if syn_block is None:
            syn_block = elem.find("synonyms")
        if syn_block is None:
            return []

        rows: list[dict[str, str]] = []
        for syn_elem in syn_block:
            if not syn_elem.text:
                continue
            synonym = syn_elem.text.strip()
            if not synonym or synonym == primary_label:
                continue
            rows.append(
                {
                    "object_id": primary_id,
                    "subject_label": synonym,
                    "subject_type": "rdfs literal",
                    "object_label": primary_label,
                    "_label_type": "alias",
                }
            )
        return rows

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

    def parse_primary_ids(
        self,
        metabolites_path: Path | str | None = None,
        proteins_path: Path | str | None = None,
    ) -> BaseMappingSet:
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

        return self.create_mapping_set([], mapping_type="id", primary_ids=primary_ids)

    def parse_primary_labels(
        self,
        metabolites_path: Path | str | None = None,
        proteins_path: Path | str | None = None,
    ) -> BaseMappingSet:
        """Return a mapping set containing the full list of current HMDB primary names.

        Args:
            metabolites_path: Path to ``hmdb_metabolites.xml`` (or zip/gz).
            proteins_path: Path to ``hmdb_proteins.xml`` (or zip/gz).

        Returns:
            :class:`~pysec2pri.parsers.base.LabelMappingSet` with
            ``_primary_labels`` populated. At least one of the two path
            arguments must be supplied.
        """
        if metabolites_path is None and proteins_path is None:
            raise ValueError("At least one of metabolites or proteins_path must be supplied.")

        primary_labels: dict[str, set[str]] = {}

        if metabolites_path is not None:
            ms_m = self.parse_synonyms(metabolites_path)
            for label, ids in object.__getattribute__(ms_m, "_primary_labels").items():
                primary_labels.setdefault(label, set()).update(ids)

        if proteins_path is not None:
            ms_p = self.parse_synonyms(proteins_path)
            for label, ids in object.__getattribute__(ms_p, "_primary_labels").items():
                primary_labels.setdefault(label, set()).update(ids)

        return self.create_mapping_set([], mapping_type="label", primary_labels=primary_labels)

    def parse(self, input_path: Path | str | None) -> BaseMappingSet:
        """Parse the datasource's XML file (or ``.zip`` / ``.gz``) into an IdMappingSet.

        Args:
            input_path: Path to ``hmdb_metabolites.xml`` or ``hmdb_proteins.xml``.

        Returns:
            :class:`BaseMappingSet` for this subclass's accessions.
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        return self._parse_xml(
            Path(input_path),
            element_tag=self.element_tag,
            prefix=self.prefix,
            desc=f"Parsing HMDB {self.element_tag}s XML",
        )

    def parse_synonyms(self, input_path: Path | str | None) -> BaseMappingSet:
        """Parse the datasource's XML file (or ``.zip`` / ``.gz``) into a LabelMappingSet.

        Args:
            input_path: Path to ``hmdb_metabolites.xml`` or ``hmdb_proteins.xml``.

        Returns:
            :class:`BaseMappingSet` for synonym-to-name mappings.
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        return self._parse_synonyms_xml(
            Path(input_path),
            element_tag=self.element_tag,
            prefix=self.prefix,
            desc=f"Parsing HMDB {self.element_tag}s XML for synonyms",
        )


# Concrete parsers


class HMDBMetaboliteParser(HMDBParser):
    """Parser for ``hmdb_metabolites.xml``.

    Reads configuration from ``hmdb_metabolites.yaml``.
    Primary accessions have the form ``HMDB0…`` (e.g. ``HMDB0000001``).
    """

    datasource_name = "hmdb_metabolites"
    element_tag = "metabolite"
    prefix = "HMDB"


class HMDBProteinParser(HMDBParser):
    """Parser for ``hmdb_proteins.xml``.

    Reads configuration from ``hmdb_proteins.yaml``.
    Primary accessions have the form ``HMDBP…`` (e.g. ``HMDBP00001``).
    """

    datasource_name = "hmdb_proteins"
    element_tag = "protein"
    prefix = "HMDBP"


__all__ = [
    "HMDBMetaboliteParser",
    "HMDBParser",
    "HMDBProteinParser",
]
