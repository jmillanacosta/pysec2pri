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
        gene_info_path: Path | str | None = None,
    ) -> Sec2PriMappingSet:
        """Parse NCBI gene_history file into an IdMappingSet.

        Args:
            input_path: Path to gene_history file (can be .gz compressed).
            tax_id: Taxonomy ID to filter by (default: "9606" for human).
            gene_info_path: Optional path to the gene_info file.  When
                supplied, ``_primary_ids`` on the returned mapping set is
                populated with every current ``NCBIGene:<id>`` CURIE for the
                given taxonomy, not just those that appear as ``object_id`` in
                a discontinued-to-primary mapping.

        Returns:
            IdMappingSet with computed cardinalities based on IDs.
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        input_path = Path(input_path)
        self._resolve_version(input_path)

        # Parse gene_history for ID mappings
        mappings = self._parse_gene_history(input_path, tax_id)

        # Create IdMappingSet and compute cardinalities
        mapping_set = self._create_mapping_set(mappings, mapping_type="id")

        # Populate the full primary ID set when gene_info is available
        if gene_info_path is not None:
            object.__setattr__(
                mapping_set,
                "_primary_ids",
                self._extract_primary_ids(Path(gene_info_path), tax_id),
            )

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
        self._resolve_version(gene_info_path)

        # Parse gene_info for symbol mappings
        mappings = self._parse_gene_info(gene_info_path, tax_id)

        # Create LabelMappingSet and compute cardinalities
        mapping_set = self._create_mapping_set(mappings, mapping_type="label")
        # Populate full primary symbol set from the same file
        object.__setattr__(
            mapping_set, "_primary_symbols", self._extract_primary_symbols(gene_info_path, tax_id)
        )
        object.__setattr__(
            mapping_set, "_primary_ids", self._extract_primary_ids(gene_info_path, tax_id)
        )
        return mapping_set

    def parse_primary_ids(
        self,
        gene_info_path: Path | str | None,
        tax_id: str = "9606",
    ) -> Sec2PriMappingSet:
        """Return a mapping set containing the full list of current NCBI Gene primary IDs.

        Reads ``gene_info`` to extract every current Gene ID for the given
        taxonomy.  The returned mapping set has an empty ``mappings`` list;
        ``_primary_ids`` is populated with every current ``NCBIGene:<id>`` CURIE.

        Args:
            gene_info_path: Path to the gene_info file (can be .gz compressed).
            tax_id: Taxonomy ID to filter by (default: ``"9606"`` for human).

        Returns:
            :class:`~pysec2pri.parsers.base.IdMappingSet` with ``_primary_ids``
            populated.
        """
        if gene_info_path is None:
            raise ValueError("gene_info_path must not be None")
        gene_info_path = Path(gene_info_path)
        self._resolve_version(gene_info_path)
        ms = self._create_mapping_set([], mapping_type="id")
        object.__setattr__(ms, "_primary_ids", self._extract_primary_ids(gene_info_path, tax_id))
        return ms

    def parse_primary_symbols(
        self,
        gene_info_path: Path | str | None,
        tax_id: str = "9606",
    ) -> Sec2PriMappingSet:
        """Return a mapping set containing the full list of current NCBI Gene symbols.

        Reads ``gene_info`` to extract every current gene symbol for the given
        taxonomy.  The returned mapping set has an empty ``mappings`` list;
        ``_primary_symbols`` is populated.

        Args:
            gene_info_path: Path to the gene_info file (can be .gz compressed).
            tax_id: Taxonomy ID to filter by (default: ``"9606"`` for human).

        Returns:
            :class:`~pysec2pri.parsers.base.LabelMappingSet` with
            ``_primary_symbols`` populated.
        """
        if gene_info_path is None:
            raise ValueError("gene_info_path must not be None")
        gene_info_path = Path(gene_info_path)
        self._resolve_version(gene_info_path)
        ms = self._create_mapping_set([], mapping_type="label")
        object.__setattr__(
            ms, "_primary_symbols", self._extract_primary_symbols(gene_info_path, tax_id)
        )
        return ms

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
        id_mappings = self.parse(gene_history_path, tax_id, gene_info_path=gene_info_path)
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
        fixed_base = {
            "mapping_justification": m_meta["mapping_justification"],
            "subject_source": m_meta.get("subject_source"),
            "object_source": m_meta.get("object_source"),
            "mapping_tool": m_meta.get("mapping_tool"),
            "license": m_meta.get("license"),
        }

        rows_data: list[dict[str, str | None]] = []
        for row in df.iter_rows(named=True):
            subject_id = str(row.get("GeneID") or "")
            object_id = str(row.get("Discontinued_GeneID") or "")
            sec_symbol = row.get("Discontinued_Symbol")
            disc_date = row.get("Discontinue_Date")

            if not object_id:
                continue

            subject_id = self.normalize_withdrawn_id(subject_id)

            if self.is_withdrawn_primary(subject_id):
                rows_data.append(
                    {
                        "subject_id": f"NCBIGene:{object_id}",
                        "object_id": WITHDRAWN_ENTRY,
                        "subject_label": str(sec_symbol) if sec_symbol else "",
                        "object_label": WITHDRAWN_ENTRY_LABEL,
                        "predicate_id": "oboInOwl:consider",
                        "comment": f"Withdrawn on {disc_date}." if disc_date else None,
                    }
                )
            else:
                rows_data.append(
                    {
                        "subject_id": f"NCBIGene:{object_id}",
                        "object_id": f"NCBIGene:{subject_id}",
                        "subject_label": str(sec_symbol) if sec_symbol else "",
                        "predicate_id": m_meta["predicate_id"],
                        "predicate_label": m_meta.get("predicate_label"),
                        "comment": f"Discontinued on {disc_date}." if disc_date else None,
                    }
                )

        return self._build_mappings(
            rows_data, fixed_base, desc="Processing gene_history", total=len(rows_data)
        )

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
        fixed = {
            "mapping_justification": m_meta["mapping_justification"],
            "subject_source": m_meta.get("subject_source"),
            "object_source": m_meta.get("object_source"),
            "mapping_tool": m_meta.get("mapping_tool"),
            "license": m_meta.get("license"),
        }

        rows_data: list[dict[str, str | None]] = []
        for row in df.iter_rows(named=True):
            gene_id = str(row.get("GeneID") or "")
            pri_symbol = row.get("Symbol")
            synonyms = row.get("Synonyms")

            if not gene_id or not pri_symbol:
                continue

            pri_symbol_str = str(pri_symbol)
            curie_id = f"NCBIGene:{gene_id}"

            if synonyms:
                for syn in str(synonyms).split("|"):
                    syn = syn.strip()
                    if syn:
                        rows_data.append(
                            {
                                "subject_id": curie_id,
                                "subject_label": pri_symbol_str,
                                "object_id": curie_id,
                                "object_label": syn,
                                "_label_type": "alias",
                                "comment": "Gene symbol synonym.",
                            }
                        )

        return self._build_mappings(
            rows_data, fixed, desc="Processing gene_info", total=len(rows_data)
        )

    def _create_mapping_set(
        self, mappings: list[Mapping], mapping_type: str = "id"
    ) -> Sec2PriMappingSet:
        """Create an IdMappingSet or LabelMappingSet with config metadata.

        Delegates to BaseParser.create_mapping_set().
        """
        return self.create_mapping_set(mappings, mapping_type)

    def _extract_primary_ids(self, file_path: Path, tax_id: str) -> set[str]:
        """Extract all current NCBI Gene IDs from gene_info for a given taxonomy.

        Args:
            file_path: Path to the gene_info file.
            tax_id: Taxonomy ID string (e.g. ``"9606"``).

        Returns:
            Set of ``NCBIGene:<id>`` CURIEs.
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
        return {f"NCBIGene:{v}" for v in df["GeneID"].drop_nulls().cast(pl.Utf8).to_list()}

    def _extract_primary_symbols(self, file_path: Path, tax_id: str) -> dict[str, set[str]]:
        """Extract all current gene symbols from gene_info for a given taxonomy.

        Returns a ``dict`` mapping each symbol text to the set of primary
        ``NCBIGene:<GeneID>`` IDs that carry it.

        Args:
            file_path: Path to the gene_info file.
            tax_id: Taxonomy ID string (e.g. ``"9606"``).

        Returns:
            ``dict[symbol, set[NCBIGene:<id>]]``
        """
        df = (
            pl.scan_csv(
                file_path,
                separator="\t",
                infer_schema_length=10000,
                null_values=["-"],
            )
            .filter(pl.col("#tax_id").cast(pl.Utf8) == tax_id)
            .select(["GeneID", "Symbol"])
            .collect()
        )
        result: dict[str, set[str]] = {}
        for gene_id, symbol in df.drop_nulls().rows():
            result.setdefault(str(symbol), set()).add(f"NCBIGene:{gene_id}")
        return result


__all__ = ["NCBIParser"]
