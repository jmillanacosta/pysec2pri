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
    ALL_SPECIES,
    WITHDRAWN_ENTRY,
    WITHDRAWN_ENTRY_LABEL,
    BaseMappingSet,
    BaseParser,
)


def _filter_by_taxon(lazy: pl.LazyFrame, species: str) -> pl.LazyFrame:
    """Filter *lazy* to ``#tax_id == species``, or return it unfiltered for :data:`ALL_SPECIES`.

    NCBI's ``gene_info``/``gene_history`` are global files covering every
    organism NCBI tracks; ``species="all"`` processes every row instead of
    one taxon, which is a much larger (and slower) operation.
    """
    if species == ALL_SPECIES:
        return lazy
    return lazy.filter(pl.col("#tax_id").cast(pl.Utf8) == species)


def _ncbi_date_to_iso(value: object) -> str | None:
    """Convert an NCBI ``YYYYMMDD`` date (str or int) to ISO ``YYYY-MM-DD``.

    Args:
        value: Raw ``Discontinue_Date`` cell, e.g. ``20230115`` or ``"20230115"``.

    Returns:
        ISO-8601 date string, or ``None`` if *value* is missing/malformed.
    """
    if value is None:
        return None
    digits = str(value)
    if len(digits) != 8 or not digits.isdigit():
        return None
    return f"{digits[:4]}-{digits[4:6]}-{digits[6:]}"


class NCBIParser(BaseParser):
    """Parser for NCBI Gene TSV files using Polars.

    Extracts secondary-to-primary NCBI Gene identifier mappings including
    gene symbols from gene_history and gene_info files.

    Returns:
    - IdMappingSet for ID-to-ID mappings (discontinued Gene IDs)
    - LabelMappingSet for symbol mappings (gene synonyms)
    """

    datasource_name = "ncbi"

    def parse(
        self,
        input_path: Path | str | None = None,
        species: str = "9606",
        gene_info_path: Path | str | None = None,
    ) -> BaseMappingSet:
        """Parse NCBI gene_history file into an IdMappingSet.

        Args:
            input_path: Path to gene_history file (can be .gz compressed).
            species: NCBI taxon ID to filter by, or "all" to skip filtering
                entirely (default: "9606" for human).
            gene_info_path: Optional path to the gene_info file.  When
                supplied, ``_primary_ids`` on the returned mapping set is
                populated with every current ``NCBIGene:<id>`` CURIE for the
                given taxonomy.

        Returns:
            IdMappingSet with computed cardinalities based on IDs.
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        input_path = Path(input_path)
        self._resolve_version(input_path)
        self.species = species

        # Parse gene_history for ID mappings
        mappings = self._parse_gene_history(input_path, species)

        # Populate the full primary ID set when gene_info is available
        primary_ids = (
            self._extract_primary_ids(Path(gene_info_path), species)
            if gene_info_path is not None
            else None
        )
        return self.create_mapping_set(mappings, mapping_type="id", primary_ids=primary_ids)

    def parse_labels(
        self,
        gene_info_path: Path | str | None,
        species: str = "9606",
    ) -> BaseMappingSet:
        """Parse NCBI gene_info file for label (label) mappings.

        Args:
            gene_info_path: Path to gene_info file.
            species: NCBI taxon ID to filter by, or "all" to skip filtering
                entirely (default: "9606" for human).

        Returns:
            LabelMappingSet with computed cardinalities based on labels.
        """
        if gene_info_path is None:
            raise ValueError("gene_info_path must not be None")
        gene_info_path = Path(gene_info_path)
        self._resolve_version(gene_info_path)
        self.species = species

        # Parse gene_info for symbol mappings
        mappings = self._parse_gene_info(gene_info_path, species)

        # Populate full primary symbol and ID sets from the same file
        return self.create_mapping_set(
            mappings,
            mapping_type="label",
            primary_labels=self._extract_primary_labels(gene_info_path, species),
            primary_ids=self._extract_primary_ids(gene_info_path, species),
        )

    def parse_primary_ids(
        self,
        gene_info_path: Path | str | None,
        species: str = "9606",
    ) -> BaseMappingSet:
        """Return a mapping set containing the full list of current NCBI Gene primary IDs.

        Reads ``gene_info`` to extract every current Gene ID for the given
        taxonomy.  The returned mapping set has an empty ``mappings`` list;
        ``_primary_ids`` is populated with every current ``NCBIGene:<id>`` CURIE.

        Args:
            gene_info_path: Path to the gene_info file (can be .gz compressed).
            species: NCBI taxon ID to filter by, or ``"all"`` to skip
                filtering entirely (default: ``"9606"`` for human).

        Returns:
            :class:`~pysec2pri.parsers.base.IdMappingSet` with ``_primary_ids``
            populated.
        """
        if gene_info_path is None:
            raise ValueError("gene_info_path must not be None")
        gene_info_path = Path(gene_info_path)
        self._resolve_version(gene_info_path)
        self.species = species
        return self.create_mapping_set(
            [], mapping_type="id", primary_ids=self._extract_primary_ids(gene_info_path, species)
        )

    def parse_primary_labels(
        self,
        gene_info_path: Path | str | None,
        species: str = "9606",
    ) -> BaseMappingSet:
        """Return a mapping set containing the full list of current NCBI Gene labels.

        Reads ``gene_info`` to extract every current gene label for the given
        taxonomy.  The returned mapping set has an empty ``mappings`` list;
        ``_primary_labels`` is populated.

        Args:
            gene_info_path: Path to the gene_info file (can be .gz compressed).
            species: NCBI taxon ID to filter by, or ``"all"`` to skip
                filtering entirely (default: ``"9606"`` for human).

        Returns:
            :class:`~pysec2pri.parsers.base.LabelMappingSet` with
            ``_primary_labels`` populated.
        """
        if gene_info_path is None:
            raise ValueError("gene_info_path must not be None")
        gene_info_path = Path(gene_info_path)
        self._resolve_version(gene_info_path)
        self.species = species
        return self.create_mapping_set(
            [],
            mapping_type="label",
            primary_labels=self._extract_primary_labels(gene_info_path, species),
        )

    def parse_all(
        self,
        gene_history_path: Path | str | None,
        gene_info_path: Path | str | None,
        species: str = "9606",
    ) -> tuple[BaseMappingSet, BaseMappingSet]:
        """Parse both gene_history and gene_info files.

        Args:
            gene_history_path: Path to gene_history file.
            gene_info_path: Path to gene_info file.
            species: NCBI taxon ID to filter by, or ``"all"`` to process
                every organism in the file (see :data:`ALL_SPECIES`).

        Returns:
            Tuple of (IdMappingSet, LabelMappingSet).
        """
        id_mappings = self.parse(gene_history_path, species, gene_info_path=gene_info_path)
        label_mappings = self.parse_labels(gene_info_path, species)
        return id_mappings, label_mappings

    def _parse_gene_history(
        self,
        file_path: Path,
        species: str,
    ) -> list[Mapping]:
        """Parse the gene_history file for ID-to-ID mappings.

        Args:
            file_path: Path to gene_history file.
            species: NCBI taxon ID to filter by, or ``"all"`` to process
                every organism in the file (see :data:`ALL_SPECIES`).

        Returns:
            List of SSSOM Mapping objects.
        """
        df = _filter_by_taxon(
            pl.scan_csv(
                file_path,
                separator="\t",
                infer_schema_length=10000,
                null_values=["-"],
            ),
            species,
        ).collect()

        if df.is_empty():
            return []

        m_meta = self.get_mapping_metadata()
        fixed_base = self._fixed_mapping_fields()

        rows_data: list[dict[str, str | None]] = []
        for row in df.iter_rows(named=True):
            subject_id = str(row.get("GeneID") or "")
            object_id = str(row.get("Discontinued_GeneID") or "")
            sec_label = row.get("Discontinued_Symbol")
            disc_date = row.get("Discontinue_Date")
            mapping_date = _ncbi_date_to_iso(disc_date)
            row_species = str(row.get("#tax_id")) if row.get("#tax_id") is not None else None

            if not object_id:
                continue

            subject_id = self.normalize_withdrawn_id(subject_id)

            if self.is_withdrawn_primary(subject_id):
                rows_data.append(
                    {
                        "subject_id": f"NCBIGene:{object_id}",
                        "object_id": WITHDRAWN_ENTRY,
                        "subject_label": str(sec_label) if sec_label else "",
                        "object_label": WITHDRAWN_ENTRY_LABEL,
                        "predicate_id": "oboInOwl:consider",
                        "comment": f"Withdrawn on {disc_date}." if disc_date else None,
                        "mapping_date": mapping_date,
                        "species": row_species,
                    }
                )
            else:
                rows_data.append(
                    {
                        "subject_id": f"NCBIGene:{object_id}",
                        "object_id": f"NCBIGene:{subject_id}",
                        "subject_label": str(sec_label) if sec_label else "",
                        "predicate_id": m_meta["predicate_id"],
                        "predicate_label": m_meta.get("predicate_label"),
                        "comment": f"Discontinued on {disc_date}." if disc_date else None,
                        "mapping_date": mapping_date,
                        "species": row_species,
                    }
                )

        return self._build_mappings(
            rows_data, fixed_base, desc="Processing gene_history", total=len(rows_data)
        )

    def _parse_gene_info(
        self,
        file_path: Path,
        species: str,
    ) -> list[Mapping]:
        """Parse the gene_info file for symbol (label) mappings.

        Args:
            file_path: Path to gene_info file.
            species: NCBI taxon ID to filter by, or ``"all"`` to process
                every organism in the file (see :data:`ALL_SPECIES`).

        Returns:
            List of SSSOM Mapping objects for symbol mappings.
        """
        df = _filter_by_taxon(
            pl.scan_csv(
                file_path,
                separator="\t",
                infer_schema_length=10000,
                null_values=["-"],
            ),
            species,
        ).collect()

        if df.is_empty():
            return []

        fixed = self._fixed_mapping_fields()

        rows_data: list[dict[str, str | None]] = []
        for row in df.iter_rows(named=True):
            gene_id = str(row.get("GeneID") or "")
            pri_label = row.get("Symbol")
            synonyms = row.get("Synonyms")

            if not gene_id or not pri_label:
                continue

            pri_label_str = str(pri_label)
            curie_id = f"NCBIGene:{gene_id}"
            row_taxon = row.get("#tax_id")
            row_species = str(row_taxon) if row_taxon is not None else None

            if synonyms:
                for syn in str(synonyms).split("|"):
                    syn = syn.strip()
                    if syn:
                        rows_data.append(
                            {
                                "object_id": curie_id,
                                "subject_label": syn,  # synonym = secondary : subject
                                "subject_type": "rdfs literal",
                                "object_label": pri_label_str,  # current label = primary : object
                                "_label_type": "alias",
                                "comment": "Gene symbol synonym.",
                                "species": row_species,
                            }
                        )

        return self._build_mappings(
            rows_data, fixed, desc="Processing gene_info", total=len(rows_data)
        )

    def _extract_primary_ids(self, file_path: Path, species: str) -> set[str]:
        """Extract all current NCBI Gene IDs from gene_info for a given taxonomy.

        Args:
            file_path: Path to the gene_info file.
            species: NCBI taxon ID string (e.g. ``"9606"``), or ``"all"``.

        Returns:
            Set of ``NCBIGene:<id>`` CURIEs.
        """
        df = _filter_by_taxon(
            pl.scan_csv(
                file_path,
                separator="\t",
                infer_schema_length=10000,
                null_values=["-"],
            ),
            species,
        ).collect()
        return {f"NCBIGene:{v}" for v in df["GeneID"].drop_nulls().cast(pl.Utf8).to_list()}

    def _extract_primary_labels(self, file_path: Path, species: str) -> dict[str, set[str]]:
        """Extract all current gene symbols from gene_info for a given taxonomy.

        Returns a ``dict`` mapping each symbol text to the set of primary
        ``NCBIGene:<GeneID>`` IDs that carry it.

        Args:
            file_path: Path to the gene_info file.
            species: NCBI taxon ID string (e.g. ``"9606"``), or ``"all"``.

        Returns:
            ``dict[label, set[NCBIGene:<id>]]``
        """
        df = (
            _filter_by_taxon(
                pl.scan_csv(
                    file_path,
                    separator="\t",
                    infer_schema_length=10000,
                    null_values=["-"],
                ),
                species,
            )
            .select(["GeneID", "Symbol"])
            .collect()
        )
        result: dict[str, set[str]] = {}
        for gene_id, label in df.drop_nulls().rows():
            result.setdefault(str(label), set()).add(f"NCBIGene:{gene_id}")
        return result


__all__ = ["NCBIParser"]
