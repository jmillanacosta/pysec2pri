"""UniProt file parser for secondary-to-primary identifier mappings.

This parser extracts ID-to-ID mappings:
- Secondary accessions -> primary accessions (from sec_ac.txt)
- Deleted accessions -> sssom:NoTermFound (from delac_sp.txt)

Uses SSSOM-compliant IdMappingSet with cardinality computation.
"""

from __future__ import annotations

import re
from pathlib import Path

import polars as pl
from sssom_schema import Mapping

from pysec2pri.parsers.base import (
    WITHDRAWN_ENTRY,
    WITHDRAWN_ENTRY_LABEL,
    BaseParser,
    Sec2PriMappingSet,
)


class UniProtParser(BaseParser):
    """Parser for UniProt files using Polars.

    Extracts secondary-to-primary UniProt accession mappings from
    sec_ac.txt (secondary accessions) and delac_sp.txt (deleted accessions).

    Returns IdMappingSet for all mappings (UniProt only has ID mappings).
    """

    datasource_name = "uniprot"

    @property
    def sec_ac_url(self) -> str:
        """Get the sec_ac.txt download URL from config."""
        return self.get_download_url("sec_ac") or ""

    @property
    def delac_url(self) -> str:
        """Get the delac_sp.txt download URL from config."""
        return self.get_download_url("delac_sp") or ""

    def parse(
        self,
        input_path: Path | str | None = None,
        delac_path: Path | str | None = None,
    ) -> Sec2PriMappingSet:
        """Parse UniProt mapping files into an IdMappingSet.

        Args:
            input_path: Path to sec_ac.txt (secondary accessions file).
            delac_path: Path to delac_sp.txt (deleted accessions file).

        Returns:
            IdMappingSet with computed cardinalities based on IDs.
        """
        mappings: list[Mapping] = []

        # Parse secondary accessions file
        if input_path is not None:
            sec_mappings = self._parse_sec_ac(Path(input_path))
            mappings.extend(sec_mappings)

        # Parse deleted accessions if provided
        if delac_path is not None:
            del_mappings = self._parse_delac(Path(delac_path))
            mappings.extend(del_mappings)

        # Create IdMappingSet and compute cardinalities
        mapping_set = self._create_mapping_set(mappings)
        return mapping_set

    def _parse_sec_ac(self, file_path: Path) -> list[Mapping]:
        """Parse sec_ac.txt for secondary -> primary accession mappings.

        Args:
            file_path: Path to sec_ac.txt file.

        Returns:
            List of SSSOM Mapping objects.
        """
        # Parse the file into a DataFrame
        lines = list(self._iter_sec_ac_data_lines(file_path))
        df = self._sec_ac_lines_to_df(lines)

        if df.is_empty():
            return []

        m_meta = self.get_mapping_metadata()
        mappings: list[Mapping] = []

        rows = list(df.iter_rows(named=True))
        for row in self._progress(rows, desc="Processing sec_ac"):
            secondary_id = row.get("object_id")
            primary_id = row.get("subject_id")

            if not secondary_id or not primary_id:
                continue

            # SSSOM: subject = primary (current), object = secondary (old)
            mapping = Mapping(
                subject_id=f"UniProtKB:{primary_id}",
                object_id=f"UniProtKB:{secondary_id}",
                predicate_id=m_meta["predicate_id"],
                mapping_justification=m_meta["mapping_justification"],
                subject_source=m_meta.get("subject_source"),
                object_source=m_meta.get("object_source"),
                mapping_tool=m_meta.get("mapping_tool"),
                license=m_meta.get("license"),
            )
            mappings.append(mapping)

        return mappings

    def _parse_delac(self, file_path: Path) -> list[Mapping]:
        """Parse delac_sp.txt for deleted accession mappings.

        Deleted accessions map to sssom:NoTermFound (1:0 cardinality).

        Args:
            file_path: Path to delac_sp.txt file.

        Returns:
            List of SSSOM Mapping objects.
        """
        m_meta = self.get_mapping_metadata()
        mappings: list[Mapping] = []

        with file_path.open("r", encoding="utf-8") as f:
            lines = [line.strip() for line in f if line.strip()]

        # Skip header lines (usually starts with "_" or empty)
        data_lines = [
            line for line in lines if line and not line.startswith("_") and not line.startswith("-")
        ]

        for line in self._progress(data_lines, desc="Processing delac"):
            # Each line is a deleted accession
            deleted_id = line.strip()
            if not deleted_id or not re.match(r"^[A-Z0-9]+$", deleted_id):
                continue

            # Deleted = withdrawn with no replacement
            mapping = Mapping(
                subject_id=WITHDRAWN_ENTRY,
                object_id=f"UniProtKB:{deleted_id}",
                subject_label=WITHDRAWN_ENTRY_LABEL,
                predicate_id="oboInOwl:consider",
                mapping_justification=m_meta["mapping_justification"],
                subject_source=m_meta.get("subject_source"),
                object_source=m_meta.get("object_source"),
                mapping_tool=m_meta.get("mapping_tool"),
                license=m_meta.get("license"),
                comment="Deleted accession with no replacement.",
            )
            mappings.append(mapping)

        return mappings

    def _iter_sec_ac_data_lines(self, file_path: Path) -> list[str]:
        """Yield data lines from sec_ac.txt, skipping headers."""
        pattern = re.compile(r"^[A-Z0-9]+\s+[A-Z0-9]+$")
        lines = []

        with file_path.open("r", encoding="utf-8") as f:
            in_data = False
            for raw_line in f:
                line = raw_line.strip()
                if not in_data:
                    if pattern.match(line):
                        in_data = True
                    else:
                        continue

                if not line or line.startswith("_"):
                    continue

                lines.append(line)

        return lines

    def _sec_ac_lines_to_df(self, lines: list[str]) -> pl.DataFrame:
        """Convert parsed sec_ac lines into a Polars DataFrame."""
        rows = []
        for line in lines:
            parts = line.split()
            if len(parts) >= 2:
                rows.append(
                    {
                        "object_id": parts[0],
                        "subject_id": parts[1],
                    }
                )

        if not rows:
            return pl.DataFrame(schema={"object_id": pl.Utf8, "subject_id": pl.Utf8})

        return pl.DataFrame(rows)

    def _create_mapping_set(
        self, mappings: list[Mapping], mapping_type: str = "id"
    ) -> Sec2PriMappingSet:
        """Delegate to base class method."""
        return self.create_mapping_set(mappings, mapping_type)


__all__ = ["UniProtParser"]
