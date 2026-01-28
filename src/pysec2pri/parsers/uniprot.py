"""UniProt file parser for secondary-to-primary identifier mappings."""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path

import polars as pl
from tqdm import tqdm

from pysec2pri.models import (
    MappingCardinality,
    MappingSet,
    UniProtMapping,
    get_comment_for_cardinality,
)
from pysec2pri.parsers.base import BaseParser


class UniProtParser(BaseParser):
    """Parser for UniProt files using Polars where applicable.

    Extracts secondary-to-primary UniProt accession mappings from
    sec_ac.txt (secondary accessions) and delac_sp.txt (deleted accessions).
    Also extracts primary IDs from FASTA files.
    """

    datasource_name = "UniProt"
    sec_ac_url = (
        "https://ftp.ebi.ac.uk/pub/databases/uniprot/"
        "current_release/knowledgebase/complete/docs/sec_ac.txt"
    )
    delac_url = (
        "https://ftp.ebi.ac.uk/pub/databases/uniprot/"
        "current_release/knowledgebase/complete/docs/delac_sp.txt"
    )
    fasta_url = (
        "https://ftp.ebi.ac.uk/pub/databases/uniprot/"
        "current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    )

    def parse(
        self,
        input_path: Path | str,
        delac_path: Path | str | None = None,
        fasta_path: Path | str | None = None,
    ) -> MappingSet:
        """Parse UniProt mapping files.

        Args:
            input_path: Path to sec_ac.txt (secondary accessions file).
            delac_path: Path to delac_sp.txt (deleted accessions file).
            fasta_path: Path to FASTA file (for primary ID extraction only).

        Returns:
            MappingSet with UniProt identifier mappings.
        """
        input_path = Path(input_path)

        mappings: list[UniProtMapping] = []

        # Parse secondary accessions file
        sec_mappings = self._parse_sec_ac(input_path)
        mappings.extend(sec_mappings)

        # Parse deleted accessions if provided
        if delac_path:
            delac_path = Path(delac_path)
            del_mappings = self._parse_delac(delac_path)
            mappings.extend(del_mappings)

        mapping_set = MappingSet(
            mappings=mappings,
            datasource_name=self.datasource_name,
            version=self.version,
            mapping_set_id="omicsfixid_uniprot_01",
            mapping_set_description=(
                "Secondary to primary ID mappings for UniProt database, "
                "generated for the omicsFixID project."
            ),
            curie_map={
                "UniProt:": "http://identifiers.org/uniprot/",
                "IAO:": "http://purl.obolibrary.org/obo/IAO_",
                "oboInOwl:": "http://www.geneontology.org/formats/oboInOwl#",
            },
        )

        return mapping_set

    def _parse_sec_ac(self, file_path: Path) -> list[UniProtMapping]:
        """Parse the sec_ac.txt file using Polars.

        The sec_ac.txt format has header lines followed by whitespace-separated
        SECONDARY PRIMARY pairs.
        """
        # Read lines and find where data starts
        with file_path.open("r", encoding="utf-8") as f:
            lines = f.readlines()

        # Find data start (skip header lines)
        data_lines = []
        in_data = False
        pattern = re.compile(r"^[A-Z0-9]+\s+[A-Z0-9]+$")

        for line in lines:
            line = line.strip()
            if not in_data:
                if pattern.match(line):
                    in_data = True
                else:
                    continue
            if not line or line.startswith("_"):
                continue
            data_lines.append(line)

        # Parse with Polars - create DataFrame from parsed data
        rows = []
        for line in data_lines:
            parts = line.split()
            if len(parts) >= 2:
                rows.append({"secondary_id": parts[0], "primary_id": parts[1]})

        if not rows:
            return []

        df = pl.DataFrame(rows)

        # Compute cardinality using Polars aggregations
        primary_counts = (
            df.group_by("primary_id")
            .agg(pl.len().alias("count"))
            .to_dict(as_series=False)
        )
        primary_count_map = dict(
            zip(primary_counts["primary_id"], primary_counts["count"])
        )

        secondary_counts = (
            df.group_by("secondary_id")
            .agg(pl.len().alias("count"))
            .to_dict(as_series=False)
        )
        secondary_count_map = dict(
            zip(secondary_counts["secondary_id"], secondary_counts["count"])
        )

        # Create mappings
        mappings: list[UniProtMapping] = []

        rows_iter = df.iter_rows(named=True)
        if self.show_progress:
            rows_iter = tqdm(
                rows_iter, total=len(df), desc="Creating UniProt mappings"
            )

        for row in rows_iter:
            primary_id = row["primary_id"]
            secondary_id = row["secondary_id"]

            # Compute cardinality
            pri_count = primary_count_map.get(primary_id, 1)
            sec_count = secondary_count_map.get(secondary_id, 1)

            if pri_count > 1 and sec_count > 1:
                cardinality = MappingCardinality.MANY_TO_MANY
            elif pri_count > 1:
                cardinality = MappingCardinality.MANY_TO_ONE
            elif sec_count > 1:
                cardinality = MappingCardinality.ONE_TO_MANY
            else:
                cardinality = MappingCardinality.ONE_TO_ONE

            comment = get_comment_for_cardinality(cardinality)
            if self.version:
                comment += f" Release: {self.version}."

            mapping = UniProtMapping(
                primary_id=primary_id,
                secondary_id=secondary_id,
                mapping_cardinality=cardinality,
                comment=comment,
                source_url=self.sec_ac_url,
            )
            mappings.append(mapping)

        return mappings

    def _parse_delac(self, file_path: Path) -> list[UniProtMapping]:
        """Parse the delac_sp.txt file for deleted/withdrawn accessions."""
        # Read lines and find where data starts
        with file_path.open("r", encoding="utf-8") as f:
            lines = f.readlines()

        # Extract valid accession IDs (skip header/footer)
        accession_pattern = re.compile(r"^[A-Z0-9]{6,10}$")
        accessions = []

        in_data = False
        for i, line in enumerate(lines):
            line = line.strip()

            if not in_data:
                if accession_pattern.match(line):
                    in_data = True
                else:
                    continue

            if not line:
                continue

            # Skip footer lines
            if i >= len(lines) - 4:
                continue

            if accession_pattern.match(line):
                accessions.append(line)

        # Create mappings using Polars DataFrame
        if not accessions:
            return []

        df = pl.DataFrame({"secondary_id": accessions})

        mappings: list[UniProtMapping] = []

        rows_iter = df.iter_rows(named=True)
        if self.show_progress:
            rows_iter = tqdm(
                rows_iter, total=len(df), desc="Creating deleted mappings"
            )

        for row in rows_iter:
            secondary_id = row["secondary_id"]

            comment = "ID (subject) withdrawn/deprecated."
            if self.version:
                comment += f" Release: {self.version}."

            mapping = UniProtMapping(
                primary_id=self.normalize_primary_id(None),
                secondary_id=secondary_id,
                mapping_cardinality=MappingCardinality.ONE_TO_ZERO,
                comment=comment,
                source_url=self.delac_url,
            )
            mappings.append(mapping)

        return mappings

    def _parse_fasta_primary_ids(self, file_path: Path) -> list[str]:
        """Extract primary IDs from a UniProt FASTA file."""
        primary_ids: list[str] = []
        pattern = re.compile(r"^>.*?\|(.+?)\|")

        opener = gzip.open if file_path.suffix == ".gz" else open
        mode = "rt" if file_path.suffix == ".gz" else "r"

        with opener(file_path, mode, encoding="utf-8") as f:
            for line in f:
                if line.startswith(">"):
                    match = pattern.match(line)
                    if match:
                        primary_ids.append(match.group(1))

        return primary_ids


__all__ = ["UniProtParser"]
