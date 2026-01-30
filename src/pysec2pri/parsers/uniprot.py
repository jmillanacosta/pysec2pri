"""UniProt file parser for secondary-to-primary identifier mappings."""

from __future__ import annotations

import gzip
import re
from collections.abc import Iterable
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
        input_path: Path | str | None = None,
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
        mappings: list[UniProtMapping] = []

        # Parse secondary accessions file
        if input_path is not None:
            sec_mappings = self._parse_sec_ac(Path(input_path))
            mappings.extend(sec_mappings)

        # Parse deleted accessions if provided
        if delac_path is not None:
            del_mappings = self._parse_delac(Path(delac_path))
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

    def _iter_sec_ac_data_lines(self, file_path: Path) -> Iterable[str]:
        """Yield data lines from sec_ac.txt, skipping headers and invalid lines."""
        pattern = re.compile(r"^[A-Z0-9]+\s+[A-Z0-9]+$")

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

                yield line

    def _sec_ac_lines_to_df(self, lines: Iterable[str]) -> pl.DataFrame:
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

    def _compute_cardinality_maps(self, df: pl.DataFrame) -> tuple[dict[str, int], dict[str, int]]:
        """Compute primary and secondary occurrence counts."""
        primary_counts = (
            df.group_by("subject_id").agg(pl.len().alias("count")).to_dict(as_series=False)
        )
        primary_map = dict(zip(primary_counts["subject_id"], primary_counts["count"], strict=False))

        secondary_counts = (
            df.group_by("object_id").agg(pl.len().alias("count")).to_dict(as_series=False)
        )
        secondary_map = dict(
            zip(
                secondary_counts["object_id"],
                secondary_counts["count"],
                strict=False,
            )
        )

        return primary_map, secondary_map

    def _infer_cardinality(
        self,
        subject_id: str,
        object_id: str,
        primary_count_map: dict[str, int],
        secondary_count_map: dict[str, int],
    ) -> MappingCardinality:
        pri_count = primary_count_map.get(subject_id, 1)
        sec_count = secondary_count_map.get(object_id, 1)

        if pri_count > 1 and sec_count > 1:
            return MappingCardinality.MANY_TO_MANY
        if pri_count > 1:
            return MappingCardinality.MANY_TO_ONE
        if sec_count > 1:
            return MappingCardinality.ONE_TO_MANY
        return MappingCardinality.ONE_TO_ONE

    def _iter_uniprot_mappings(
        self,
        df: pl.DataFrame,
        primary_count_map: dict[str, int],
        secondary_count_map: dict[str, int],
    ) -> Iterable[UniProtMapping]:
        rows_iter = df.iter_rows(named=True)

        if self.show_progress:
            rows_iter = tqdm(
                rows_iter,
                total=len(df),
                desc="Creating UniProt mappings",
            )

        for row in rows_iter:
            subject_id = row["subject_id"]
            object_id = row["object_id"]

            cardinality = self._infer_cardinality(
                subject_id,
                object_id,
                primary_count_map,
                secondary_count_map,
            )

            # Use IAO:0100001 (term replaced by) for clear replacements
            predicate = "IAO:0100001"
            if cardinality in (
                MappingCardinality.ONE_TO_MANY,
                MappingCardinality.MANY_TO_MANY,
                MappingCardinality.ONE_TO_ZERO,
            ):
                predicate = "oboInOwl:consider"

            comment = get_comment_for_cardinality(cardinality)
            if self.version:
                comment += f" Release: {self.version}."

            yield UniProtMapping(
                subject_id=subject_id,
                predicate_id=predicate,
                object_id=object_id,
                mapping_cardinality=cardinality,
                comment=comment,
                source_url=self.sec_ac_url,
            )

    def _parse_sec_ac(self, file_path: Path) -> list[UniProtMapping]:
        """Parse the sec_ac.txt file and return UniProt mappings."""
        lines = self._iter_sec_ac_data_lines(file_path)
        df = self._sec_ac_lines_to_df(lines)

        if df.is_empty():
            return []

        primary_map, secondary_map = self._compute_cardinality_maps(df)

        return list(self._iter_uniprot_mappings(df, primary_map, secondary_map))

    def _iter_delac_accessions(self, file_path: Path) -> Iterable[str]:
        """Yield valid deleted/withdrawn accession IDs from delac_sp.txt."""
        accession_pattern = re.compile(r"^[A-Z0-9]{6,10}$")

        with file_path.open("r", encoding="utf-8") as f:
            lines = f.readlines()

        in_data = False
        footer_start = max(len(lines) - 4, 0)

        for i, raw_line in enumerate(lines):
            line = raw_line.strip()

            if not in_data:
                if accession_pattern.match(line):
                    in_data = True
                else:
                    continue

            if not line:
                continue

            if i >= footer_start:
                continue

            if accession_pattern.match(line):
                yield line

    def _delac_accessions_to_df(self, accessions: Iterable[str]) -> pl.DataFrame:
        """Convert deleted accession IDs into a Polars DataFrame."""
        acc_list = list(accessions)
        if not acc_list:
            return pl.DataFrame(schema={"object_id": pl.Utf8})

        return pl.DataFrame({"object_id": acc_list})

    def _iter_deleted_mappings(
        self,
        df: pl.DataFrame,
    ) -> Iterable[UniProtMapping]:
        rows_iter = df.iter_rows(named=True)

        if self.show_progress:
            rows_iter = tqdm(
                rows_iter,
                total=len(df),
                desc="Creating deleted mappings",
            )

        for row in rows_iter:
            object_id = row["object_id"]

            comment = "ID (subject) withdrawn/deprecated."
            if self.version:
                comment += f" Release: {self.version}."

            yield UniProtMapping(
                subject_id=self.normalize_subject_id(None),
                predicate_id="oboInOwl:consider",
                object_id=object_id,
                mapping_cardinality=MappingCardinality.ONE_TO_ZERO,
                comment=comment,
                source_url=self.delac_url,
            )

    def _parse_delac(self, file_path: Path) -> list[UniProtMapping]:
        """Parse the delac_sp.txt file for deleted/withdrawn accessions."""
        accessions = self._iter_delac_accessions(file_path)
        df = self._delac_accessions_to_df(accessions)

        if df.is_empty():
            return []

        return list(self._iter_deleted_mappings(df))

    def _parse_fasta_subject_ids(self, file_path: Path) -> list[str]:
        """Extract primary IDs from a UniProt FASTA file."""
        subject_ids: list[str] = []
        pattern = re.compile(r"^>.*?\|(.+?)\|")

        opener = gzip.open if file_path.suffix == ".gz" else open
        # Always open in text mode for str lines
        with opener(file_path, "rt", encoding="utf-8") as f:
            for line in f:
                if line.startswith(">"):
                    match = pattern.match(line)
                    if match:
                        subject_ids.append(match.group(1))
        return subject_ids


__all__ = ["UniProtParser"]
