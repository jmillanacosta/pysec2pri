"""NCBI Gene TSV file parser using Polars for memory efficiency."""

from __future__ import annotations

import gzip
from collections import defaultdict
from pathlib import Path

import polars as pl
from tqdm import tqdm

from pysec2pri.models import (
    MappingCardinality,
    MappingSet,
    NCBIGeneMapping,
    get_comment_for_cardinality,
)
from pysec2pri.parsers.base import BaseParser


class NCBIParser(BaseParser):
    """Parser for NCBI Gene TSV files using Polars.

    Extracts secondary-to-primary NCBI Gene identifier mappings including
    gene symbols from gene_history and gene_info files.
    """

    datasource_name = "NCBI"
    history_source_url = "https://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz"
    info_source_url = "https://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz"

    def parse(
        self,
        input_path: Path | str,
        gene_info_path: Path | str | None = None,
        tax_id: str = "9606",
    ) -> MappingSet:
        """Parse NCBI Gene files.

        Args:
            input_path: Path to gene_history file (can be .gz compressed).
            gene_info_path: Path to gene_info file for symbol information.
            tax_id: Taxonomy ID to filter by (default: "9606" for human).

        Returns:
            MappingSet with NCBI Gene identifier and symbol mappings.
        """
        input_path = Path(input_path)

        mappings: list[NCBIGeneMapping] = []

        # Parse gene_history file
        history_mappings, symbol_info = self._parse_gene_history(
            input_path, tax_id
        )
        mappings.extend(history_mappings)

        # Parse gene_info file if provided
        if gene_info_path:
            gene_info_path = Path(gene_info_path)
            symbol_mappings = self._parse_gene_info(
                gene_info_path, tax_id, symbol_info
            )
            mappings.extend(symbol_mappings)

        mapping_set = MappingSet(
            mappings=mappings,
            datasource_name=self.datasource_name,
            version=self.version,
            mapping_set_id="omicsfixid_ncbi_01",
            mapping_set_description=(
                "Secondary to primary ID mappings for NCBI database, "
                "generated for the omicsFixID project."
            ),
            comment=(
                "object_label represents an alternative symbol "
                "by which the subject_label is also known."
            ),
            curie_map={
                "NCBIGene:": "http://www.ncbi.nlm.nih.gov/gene/",
                "IAO:": "http://purl.obolibrary.org/obo/IAO_",
                "oboInOwl:": "http://www.geneontology.org/formats/oboInOwl#",
            },
        )

        return mapping_set

    def _decompress_if_needed(self, file_path: Path) -> Path:
        """Decompress gzip file if needed, return path to readable file."""
        if file_path.suffix == ".gz":
            # Polars can read gzip directly, but for very large files
            # it's sometimes faster to decompress first
            return file_path
        return file_path

    def _parse_gene_history(
        self,
        file_path: Path,
        tax_id: str,
    ) -> tuple[list[NCBIGeneMapping], dict[str, str]]:
        """Parse the gene_history file using Polars.

        Uses streaming to handle very large files efficiently.
        """
        # Read the file with Polars (handles gzip automatically)
        # Use lazy evaluation for memory efficiency
        df = pl.scan_csv(
            file_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=["-"],
        )

        # Filter by tax_id and collect
        df_filtered = df.filter(
            pl.col("#tax_id").cast(pl.Utf8) == tax_id
        ).collect()

        if self.show_progress:
            print(f"Filtered to {len(df_filtered)} rows for tax_id {tax_id}")

        # Build cardinality counts
        primary_counts: dict[str, int] = defaultdict(int)
        secondary_counts: dict[str, int] = defaultdict(int)
        all_secondary_ids: set[str] = set()
        symbol_info: dict[str, str] = {}

        # First pass: compute cardinality
        rows = df_filtered.iter_rows(named=True)
        if self.show_progress:
            rows = tqdm(
                rows,
                total=len(df_filtered),
                desc="Pass 1: Computing cardinality",
            )

        for row in rows:
            primary_id = str(row.get("GeneID") or "")
            secondary_id = str(row.get("Discontinued_GeneID") or "")
            sec_symbol = row.get("Discontinued_Symbol")

            primary_id = self.normalize_primary_id(primary_id)

            primary_counts[primary_id] += 1
            secondary_counts[secondary_id] += 1
            all_secondary_ids.add(secondary_id)
            if sec_symbol:
                symbol_info[secondary_id] = str(sec_symbol)

        # Second pass: create mappings
        mappings: list[NCBIGeneMapping] = []

        rows = df_filtered.iter_rows(named=True)
        if self.show_progress:
            rows = tqdm(
                rows,
                total=len(df_filtered),
                desc="Pass 2: Creating mappings",
            )

        for row in rows:
            primary_id = str(row.get("GeneID") or "")
            secondary_id = str(row.get("Discontinued_GeneID") or "")
            sec_symbol = row.get("Discontinued_Symbol")
            disc_date = row.get("Discontinue_Date")

            primary_id = self.normalize_primary_id(primary_id)

            # Compute cardinality
            pri_count = primary_counts.get(primary_id, 1)
            sec_count = secondary_counts.get(secondary_id, 1)

            if pri_count > 1 and sec_count > 1:
                cardinality = MappingCardinality.MANY_TO_MANY
            elif pri_count > 1:
                cardinality = MappingCardinality.MANY_TO_ONE
            elif sec_count > 1:
                cardinality = MappingCardinality.ONE_TO_MANY
            else:
                cardinality = MappingCardinality.ONE_TO_ONE

            if self.is_withdrawn_primary(primary_id):
                cardinality = MappingCardinality.ONE_TO_ZERO

            # Build comment
            comment_parts = []
            if disc_date:
                comment_parts.append(f"Withdrawn date: {disc_date}.")
            comment_parts.append(get_comment_for_cardinality(cardinality))
            if primary_id in all_secondary_ids:
                comment_parts.append("Object is also withdrawn.")
            if self.version:
                comment_parts.append(f"Release: {self.version}.")

            mapping = NCBIGeneMapping(
                primary_id=primary_id,
                secondary_id=secondary_id,
                primary_symbol=None,
                secondary_symbol=str(sec_symbol) if sec_symbol else None,
                mapping_cardinality=cardinality,
                comment=" ".join(comment_parts),
                source_url=self.history_source_url,
            )
            mappings.append(mapping)

        return mappings, symbol_info

    def _parse_gene_info(
        self,
        file_path: Path,
        tax_id: str,
        symbol_info: dict[str, str],
    ) -> list[NCBIGeneMapping]:
        """Parse the gene_info file for symbol mappings using Polars."""
        # Read with Polars lazy evaluation
        df = pl.scan_csv(
            file_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=["-"],
        )

        # Filter by tax_id
        df_filtered = df.filter(
            pl.col("#tax_id").cast(pl.Utf8) == tax_id
        ).collect()

        if self.show_progress:
            print(f"Filtered gene_info to {len(df_filtered)} rows")

        mappings: list[NCBIGeneMapping] = []

        rows = df_filtered.iter_rows(named=True)
        if self.show_progress:
            rows = tqdm(
                rows,
                total=len(df_filtered),
                desc="Processing gene_info",
            )

        for row in rows:
            primary_id = str(row.get("GeneID") or "")
            pri_symbol = row.get("Symbol")
            synonyms = row.get("Synonyms")
            auth_sym = row.get("Symbol_from_nomenclature_authority")

            if not primary_id:
                continue

            # Handle nomenclature authority symbol
            pri_symbol_str = str(pri_symbol) if pri_symbol else ""
            if auth_sym and str(auth_sym) != pri_symbol_str:
                pri_symbol_str = f"{pri_symbol_str}|{auth_sym}"

            # Process synonyms
            if synonyms:
                synonyms_str = str(synonyms)
                for syn in synonyms_str.split("|"):
                    syn = syn.strip()
                    if syn:
                        comment = "Unofficial symbol for the gene."
                        if self.version:
                            comment += f" Release: {self.version}."

                        mapping = NCBIGeneMapping(
                            primary_id=primary_id,
                            secondary_id=None,
                            primary_symbol=pri_symbol_str if pri_symbol_str else None,
                            secondary_symbol=syn,
                            comment=comment,
                            source_url=self.info_source_url,
                        )
                        mappings.append(mapping)

        return mappings


__all__ = ["NCBIParser"]
