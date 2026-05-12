"""UniProt file parser for secondary-to-primary identifier mappings.

This parser extracts ID-to-ID mappings:
- Secondary accessions -> primary accessions (from sec_ac.txt)
- Deleted accessions -> sssom:NoTermFound (from delac_sp.txt)

Uses SSSOM-compliant IdMappingSet with cardinality computation.
"""

from __future__ import annotations

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
        # Resolve version
        self._resolve_version(
            Path(input_path)
            if input_path is not None
            else (Path(delac_path) if delac_path is not None else None)
        )

        mappings: list[Mapping] = []

        if input_path is not None:
            mappings.extend(self._parse_sec_ac(Path(input_path)))

        if delac_path is not None:
            mappings.extend(self._parse_delac(Path(delac_path)))

        return self._create_mapping_set(mappings)

    def _parse_sec_ac(self, file_path: Path) -> list[Mapping]:
        """Parse sec_ac.txt for secondary -> primary accession mappings.

        Args:
            file_path: Path to sec_ac.txt file.

        Returns:
            List of SSSOM Mapping objects.
        """
        # Count lines until (and including) the separator row.
        skip_rows = 0
        with file_path.open("r", encoding="utf-8") as f:
            for raw_line in f:
                skip_rows += 1
                if raw_line.startswith("_"):
                    break  # next line is first data row

        df = (
            pl.scan_csv(
                file_path,
                has_header=False,
                skip_rows=skip_rows,
                separator="\n",
                new_columns=["line"],
                infer_schema_length=0,
                quote_char=None,
            )
            .filter(
                pl.col("line").str.len_chars() > 0,
                ~pl.col("line").str.starts_with("-"),
                ~pl.col("line").str.starts_with("_"),
            )
            .with_columns(
                pl.col("line")
                .str.split_exact(" ", 1)
                .struct.field("field_0")
                .str.strip_chars()
                .alias("subject_id"),
                pl.col("line").str.split(" ").list.last().str.strip_chars().alias("object_id"),
            )
            .filter(
                pl.col("subject_id").str.len_chars() > 0,
                pl.col("object_id").str.len_chars() > 0,
                pl.col("subject_id") != pl.col("object_id"),
            )
            .with_columns(
                (pl.lit("UniProtKB:") + pl.col("subject_id")).alias("subject_id"),
                (pl.lit("UniProtKB:") + pl.col("object_id")).alias("object_id"),
            )
            .select(["subject_id", "object_id"])
            .collect()
        )

        if df.is_empty():
            return []

        m_meta = self.get_mapping_metadata()
        fixed = {
            "predicate_id": m_meta["predicate_id"],
            "predicate_label": m_meta.get("predicate_label"),
            "mapping_justification": m_meta["mapping_justification"],
            "subject_source": m_meta.get("subject_source"),
            "object_source": m_meta.get("object_source"),
            "mapping_tool": m_meta.get("mapping_tool"),
            "license": m_meta.get("license"),
        }
        rows = df.select(["subject_id", "object_id"]).to_dicts()
        return self._build_mappings(rows, fixed, desc="Processing sec_ac", total=len(rows))

    def _parse_delac(self, file_path: Path) -> list[Mapping]:
        """Parse delac_sp.txt for deleted accession mappings.

        Deleted accessions map to sssom:NoTermFound (1:0 cardinality).

        Args:
            file_path: Path to delac_sp.txt file.

        Returns:
            List of SSSOM Mapping objects.
        """
        skip_rows = 0
        with file_path.open("r", encoding="utf-8") as f:
            for raw_line in f:
                skip_rows += 1
                if raw_line.startswith("_"):
                    break  # next line is first deleted accession

        df = (
            pl.scan_csv(
                file_path,
                has_header=False,
                skip_rows=skip_rows,
                separator="\t",
                new_columns=["accession"],
                infer_schema_length=0,
                quote_char=None,
            )
            .with_columns(pl.col("accession").str.strip_chars())
            .filter(
                pl.col("accession").str.contains(
                    r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
                )
            )
            .with_columns((pl.lit("UniProtKB:") + pl.col("accession")).alias("object_id"))
            .select("object_id")
            .collect()
        )

        if df.is_empty():
            return []

        m_meta = self.get_mapping_metadata()
        fixed = {
            "subject_id": WITHDRAWN_ENTRY,
            "subject_label": WITHDRAWN_ENTRY_LABEL,
            "predicate_id": "oboInOwl:consider",
            "mapping_justification": m_meta["mapping_justification"],
            "subject_source": m_meta.get("subject_source"),
            "object_source": m_meta.get("object_source"),
            "mapping_tool": m_meta.get("mapping_tool"),
            "license": m_meta.get("license"),
            "comment": "Deleted accession with no replacement.",
        }
        rows = df.select("object_id").to_dicts()
        return self._build_mappings(rows, fixed, desc="Processing delac", total=len(rows))

    def _create_mapping_set(
        self, mappings: list[Mapping], mapping_type: str = "id"
    ) -> Sec2PriMappingSet:
        """Delegate to base class method."""
        return self.create_mapping_set(mappings, mapping_type)


__all__ = ["UniProtParser"]
