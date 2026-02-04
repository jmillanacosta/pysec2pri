"""NCBI Gene TSV file parser for secondary-to-primary identifier mappings.

This parser extracts:
1. ID-to-ID mappings: discontinued Gene IDs -> current Gene IDs
2. Label-to-label mappings: gene symbol synonyms -> current symbols

Uses SSSOM-compliant MappingSet classes with cardinality computation.
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


class NCBIParser(BaseParser):
    """Parser for NCBI Gene TSV files using Polars.

    Extracts secondary-to-primary NCBI Gene identifier mappings including
    gene symbols from gene_history and gene_info files.

    Returns:
    - IdMappingSet for ID-to-ID mappings (discontinued Gene IDs)
    - LabelMappingSet for symbol mappings (gene synonyms)
    """

    datasource_name = "ncbi"

    @property
    def history_source_url(self) -> str:
        """Get the gene_history download URL from config."""
        return self.get_download_url("gene_history") or ""

    @property
    def info_source_url(self) -> str:
        """Get the gene_info download URL from config."""
        return self.get_download_url("gene_info") or ""

    def parse(
        self,
        input_path: Path | str | None = None,
        tax_id: str = "9606",
    ) -> Sec2PriMappingSet:
        """Parse NCBI gene_history file into an IdMappingSet.

        Args:
            input_path: Path to gene_history file (can be .gz compressed).
            tax_id: Taxonomy ID to filter by (default: "9606" for human).

        Returns:
            IdMappingSet with computed cardinalities based on IDs.
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        input_path = Path(input_path)

        # Parse gene_history for ID mappings
        mappings = self._parse_gene_history(input_path, tax_id)

        # Create IdMappingSet and compute cardinalities
        mapping_set = self._create_mapping_set(mappings, mapping_type="id")
        return mapping_set

    def parse_symbols(
        self,
        gene_info_path: Path | str | None,
        tax_id: str = "9606",
    ) -> Sec2PriMappingSet:
        """Parse NCBI gene_info file for symbol (label) mappings.

        Args:
            gene_info_path: Path to gene_info file.
            tax_id: Taxonomy ID to filter by (default: "9606" for human).

        Returns:
            LabelMappingSet with computed cardinalities based on labels.
        """
        if gene_info_path is None:
            raise ValueError("gene_info_path must not be None")
        gene_info_path = Path(gene_info_path)

        # Parse gene_info for symbol mappings
        mappings = self._parse_gene_info(gene_info_path, tax_id)

        # Create LabelMappingSet and compute cardinalities
        mapping_set = self._create_mapping_set(mappings, mapping_type="label")
        return mapping_set

    def parse_all(
        self,
        gene_history_path: Path | str | None,
        gene_info_path: Path | str | None,
        tax_id: str = "9606",
    ) -> tuple[Sec2PriMappingSet, Sec2PriMappingSet]:
        """Parse both gene_history and gene_info files.

        Args:
            gene_history_path: Path to gene_history file.
            gene_info_path: Path to gene_info file.
            tax_id: Taxonomy ID to filter by.

        Returns:
            Tuple of (IdMappingSet, LabelMappingSet).
        """
        id_mappings = self.parse(gene_history_path, tax_id)
        label_mappings = self.parse_symbols(gene_info_path, tax_id)
        return id_mappings, label_mappings

    def _parse_gene_history(
        self,
        file_path: Path,
        tax_id: str,
    ) -> list[Mapping]:
        """Parse the gene_history file for ID-to-ID mappings.

        Args:
            file_path: Path to gene_history file.
            tax_id: Taxonomy ID to filter by.

        Returns:
            List of SSSOM Mapping objects.
        """
        df = (
            pl.scan_csv(
                file_path,
                separator="\t",
                infer_schema_length=10000,
                null_values=["-"],
            )
            .filter(pl.col("#tax_id").cast(pl.Utf8) == tax_id)
            .collect()
        )

        if df.is_empty():
            return []

        m_meta = self.get_mapping_metadata()
        mappings: list[Mapping] = []

        rows = list(df.iter_rows(named=True))
        for row in self._progress(rows, desc="Processing gene_history"):
            subject_id = str(row.get("GeneID") or "")
            object_id = str(row.get("Discontinued_GeneID") or "")
            sec_symbol = row.get("Discontinued_Symbol")
            disc_date = row.get("Discontinue_Date")

            if not object_id:
                continue

            # Normalize withdrawn IDs
            subject_id = self.normalize_withdrawn_id(subject_id)

            # Determine if this is a withdrawn entry with no replacement
            if self.is_withdrawn_primary(subject_id):
                mapping = Mapping(
                    subject_id=WITHDRAWN_ENTRY,
                    object_id=f"NCBIGene:{object_id}",
                    subject_label=WITHDRAWN_ENTRY_LABEL,
                    object_label=str(sec_symbol) if sec_symbol else "",
                    predicate_id="oboInOwl:consider",
                    mapping_justification=m_meta["mapping_justification"],
                    subject_source=m_meta.get("subject_source"),
                    object_source=m_meta.get("object_source"),
                    mapping_tool=m_meta.get("mapping_tool"),
                    license=m_meta.get("license"),
                    comment=f"Withdrawn on {disc_date}." if disc_date else None,
                )
            else:
                # Normal replacement mapping
                mapping = Mapping(
                    subject_id=f"NCBIGene:{subject_id}",
                    object_id=f"NCBIGene:{object_id}",
                    object_label=str(sec_symbol) if sec_symbol else "",
                    predicate_id=m_meta["predicate_id"],
                    mapping_justification=m_meta["mapping_justification"],
                    subject_source=m_meta.get("subject_source"),
                    object_source=m_meta.get("object_source"),
                    mapping_tool=m_meta.get("mapping_tool"),
                    license=m_meta.get("license"),
                    comment=f"Discontinued on {disc_date}." if disc_date else None,
                )
            mappings.append(mapping)

        return mappings

    def _parse_gene_info(
        self,
        file_path: Path,
        tax_id: str,
    ) -> list[Mapping]:
        """Parse the gene_info file for symbol (label) mappings.

        Args:
            file_path: Path to gene_info file.
            tax_id: Taxonomy ID to filter by.

        Returns:
            List of SSSOM Mapping objects for symbol mappings.
        """
        df = (
            pl.scan_csv(
                file_path,
                separator="\t",
                infer_schema_length=10000,
                null_values=["-"],
            )
            .filter(pl.col("#tax_id").cast(pl.Utf8) == tax_id)
            .collect()
        )

        if df.is_empty():
            return []

        m_meta = self.get_mapping_metadata()
        mappings: list[Mapping] = []

        rows = list(df.iter_rows(named=True))
        for row in self._progress(rows, desc="Processing gene_info"):
            gene_id = str(row.get("GeneID") or "")
            pri_symbol = row.get("Symbol")
            synonyms = row.get("Synonyms")

            if not gene_id or not pri_symbol:
                continue

            pri_symbol_str = str(pri_symbol)
            curie_id = f"NCBIGene:{gene_id}"

            # Process synonyms
            if synonyms:
                synonyms_str = str(synonyms)
                for syn in synonyms_str.split("|"):
                    syn = syn.strip()
                    if syn:
                        mapping = Mapping(
                            subject_id=curie_id,
                            subject_label=pri_symbol_str,
                            object_id=curie_id,  # Same entity
                            object_label=syn,
                            predicate_id="oboInOwl:hasRelatedSynonym",
                            mapping_justification=m_meta["mapping_justification"],
                            subject_source=m_meta.get("subject_source"),
                            object_source=m_meta.get("object_source"),
                            mapping_tool=m_meta.get("mapping_tool"),
                            license=m_meta.get("license"),
                            comment="Gene symbol synonym.",
                        )
                        mappings.append(mapping)

        return mappings

    def _create_mapping_set(
        self, mappings: list[Mapping], mapping_type: str = "id"
    ) -> Sec2PriMappingSet:
        """Create an IdMappingSet or LabelMappingSet with config metadata.

        Delegates to BaseParser.create_mapping_set().
        """
        return self.create_mapping_set(mappings, mapping_type)


__all__ = ["NCBIParser"]
