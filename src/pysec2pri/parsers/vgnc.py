"""VGNC TSV file parser for secondary-to-primary identifier mappings.

- ID-to-ID mappings: withdrawn/merged VGNC IDs -> current VGNC IDs.
- Label-to-label mappings: previous/alias gene symbols -> current symbols,
   from the gene-set file.
"""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

import polars as pl
from sssom_schema import Mapping

from pysec2pri.parsers._nomenclature import (
    ALIAS_SYMBOL,
    DATE_SYMBOL_CHANGED,
    PREV_SYMBOL,
    STATUS,
    SYMBOL,
    GeneNomenclatureParser,
)
from pysec2pri.parsers.base import (
    ALL_SPECIES,
    BaseMappingSet,
    LabelMappingSet,
)

# VGNC column names (case-insensitive matching used)
VGNC_ID = "vgnc_id"
WITHDRAWN_SYMBOL = "withdrawn_symbol"
TAXON_ID = "taxon_id"


class VGNCParser(GeneNomenclatureParser):
    """Parser for VGNC TSV files using Polars for memory efficiency.

    Extracts secondary-to-primary VGNC identifier mappings and symbol
    mappings from the VGNC withdrawn and gene-set files.

    Produces two mapping-set kinds:

    - IdMappingSet for ID-to-ID mappings (withdrawn/merged IDs).
    - LabelMappingSet for symbol mappings (alias/previous symbols), scoped to
      one species.
    """

    datasource_name = "vgnc"
    id_column: ClassVar[str] = VGNC_ID
    withdrawn_label_column: ClassVar[str] = WITHDRAWN_SYMBOL

    def __init__(
        self,
        version: str | None = None,
        show_progress: bool = True,
    ):
        """Initialize the VGNC parser.

        Args:
            version: Version/release identifier for the datasource.
            show_progress: Whether to show progress bars during parsing.
        """
        super().__init__(version=version, show_progress=show_progress)

    def parse(
        self,
        input_path: Path | str | None,
        complete_set_path: Path | str | None = None,
        species: str | None = None,
    ) -> BaseMappingSet:
        """Parse the VGNC withdrawn TSV file into an IdMappingSet.

        Always parses the *full*, unfiltered withdrawn file first (see
        module docstring); when *species* is given (and isn't
        :data:`ALL_SPECIES`), the result is then subset by resolving each
        mapping's primary VGNC ID against the gene-set file's ``taxon_id``
        column -- this requires *complete_set_path*.

        Args:
            input_path: Path to the VGNC withdrawn TSV file.
            complete_set_path: Path to the VGNC gene-set file. Required
                when *species* is given (used to resolve taxon IDs). When
                supplied, ``_primary_ids`` on the returned mapping set is
                populated with every current VGNC ID for *species* (or
                across all species, when *species* is ``None``/``"all"``).
            species: NCBI taxon ID to subset the output to, or
                :data:`ALL_SPECIES`. ``None`` (default) returns the full,
                unfiltered set across every species.

        Returns:
            IdMappingSet with computed cardinalities based on IDs.

        Raises:
            ValueError: If *species* is given without *complete_set_path*
                (there is no other way to resolve taxon IDs).
        """
        if input_path is None:
            raise ValueError("input_path must not be None")
        if species not in (None, ALL_SPECIES) and complete_set_path is None:
            raise ValueError("species subsetting requires complete_set_path to resolve taxon IDs.")
        input_path = Path(input_path)
        self._resolve_version(input_path)
        if species not in (None, ALL_SPECIES):
            self.species = species

        taxon_by_id: dict[str, str] | None = None
        if complete_set_path is not None:
            taxon_by_id = self._taxon_by_vgnc_id(Path(complete_set_path))

        mappings = self._parse_withdrawn(input_path, taxon_by_id=taxon_by_id)

        if species not in (None, ALL_SPECIES) and taxon_by_id is not None:
            mappings = [
                m
                for m in mappings
                if taxon_by_id.get(str(getattr(m, "object_id", "") or "")) == str(species)
            ]

        primary_ids: set[str] | None = None
        if complete_set_path is not None:
            primary_ids = self._extract_primary_ids(Path(complete_set_path))
            if species not in (None, ALL_SPECIES) and taxon_by_id is not None:
                primary_ids = {vid for vid in primary_ids if taxon_by_id.get(vid) == str(species)}
        return self.create_mapping_set(mappings, mapping_type="id", primary_ids=primary_ids)

    def parse_primary_ids(
        self,
        complete_set_path: Path | str | None,
        species: str | None = None,
    ) -> BaseMappingSet:
        """Return a mapping set whose only content is the full primary ID list.

        Reads the VGNC gene-set file to extract every current VGNC ID,
        optionally subset to *species*.

        Args:
            complete_set_path: Path to the VGNC gene-set TSV file.
            species: NCBI taxon ID to subset the result to, or
                :data:`ALL_SPECIES`. ``None`` (default) returns the full,
                unfiltered set across every species.

        Returns:
            :class:`~pysec2pri.parsers.base.IdMappingSet` with no mappings
            and ``_primary_ids`` populated.
        """
        if complete_set_path is None:
            raise ValueError("complete_set_path must not be None")
        complete_set_path = Path(complete_set_path)
        self._resolve_version(complete_set_path)
        if species not in (None, ALL_SPECIES):
            self.species = species

        primary_ids = self._extract_primary_ids(complete_set_path)
        if species not in (None, ALL_SPECIES):
            taxon_by_id = self._taxon_by_vgnc_id(complete_set_path)
            primary_ids = {vid for vid in primary_ids if taxon_by_id.get(vid) == str(species)}

        return self.create_mapping_set([], mapping_type="id", primary_ids=primary_ids)

    def parse_labels(
        self,
        complete_set_path: Path | str | None,
        species: str,
        statuses: list[str] | None = None,
    ) -> LabelMappingSet:
        """Parse the VGNC gene-set file for symbol (label) mappings, scoped to one species.

        Args:
            complete_set_path: Path to the VGNC gene-set TSV file.
            species: NCBI taxon ID to filter by, or :data:`ALL_SPECIES` to
                process every species together (see module docstring for
                why that changes ambiguity detection). Required at this
                layer -- callers needing config's ``species.default``
                fallback (see ``config/vgnc.yaml``) should resolve it
                themselves, as :mod:`pysec2pri.api`'s
                ``generate_vgnc_labels`` does.
            statuses: Entry statuses to include (e.g. ``["Approved"]``).
                If ``None`` (default), all entries are included.

        Returns:
            LabelMappingSet with computed cardinalities based on labels.
        """
        if complete_set_path is None:
            raise ValueError("complete_set_path must not be None")
        complete_set_path = Path(complete_set_path)
        self._resolve_version(complete_set_path)
        self.species = species

        mappings = self._parse_gene_set(complete_set_path, species, statuses=statuses)
        return self.create_mapping_set(
            mappings,
            mapping_type="label",
            primary_labels=self._extract_primary_labels(complete_set_path, species),
        )

    def parse_primary_labels(
        self,
        complete_set_path: Path | str | None,
        species: str,
    ) -> BaseMappingSet:
        """Return a mapping set whose only content is the full primary Symbol list.

        Reads the VGNC gene-set file to extract every current approved
        symbol for *species*, storing it in ``_primary_labels``. The
        ``mappings`` list is intentionally left empty; this mapping set
        exists only to drive ``to_pri_labels()``.

        Args:
            complete_set_path: Path to the VGNC gene-set TSV file.
            species: NCBI taxon ID to filter by, or :data:`ALL_SPECIES`.

        Returns:
            :class:`~pysec2pri.parsers.base.LabelMappingSet` with no
            mappings and ``_primary_labels`` populated.
        """
        if complete_set_path is None:
            raise ValueError("complete_set_path must not be None")
        complete_set_path = Path(complete_set_path)
        self._resolve_version(complete_set_path)
        self.species = species

        return self.create_mapping_set(
            [],
            mapping_type="label",
            primary_labels=self._extract_primary_labels(complete_set_path, species),
        )

    def parse_all(
        self,
        withdrawn_path: Path | str | None,
        complete_set_path: Path | str | None,
        species: str,
    ) -> tuple[BaseMappingSet, BaseMappingSet]:
        """Parse both the withdrawn and gene-set files.

        Args:
            withdrawn_path: Path to the VGNC withdrawn TSV file.
            complete_set_path: Path to the VGNC gene-set TSV file.
            species: NCBI taxon ID to filter the label mappings by.

        Returns:
            Tuple of (IdMappingSet, LabelMappingSet).
        """
        id_mappings = self.parse(withdrawn_path, complete_set_path=complete_set_path)
        label_mappings = self.parse_labels(complete_set_path, species)
        return id_mappings, label_mappings

    def _taxon_by_vgnc_id(self, file_path: Path) -> dict[str, str]:
        """Return ``{vgnc_id: taxon_id}`` for every current gene in the gene-set file.

        Used to resolve which species a withdrawn entry's *replacement*
        gene belongs to, since the withdrawn file ``taxon_id``
        column is not populated upstream (see module docstring).

        Args:
            file_path: Path to the VGNC gene-set TSV file.

        Returns:
            ``dict[vgnc_id, taxon_id]``.
        """
        df = self._read_tsv(file_path)
        id_col = self._find_column(df.columns, self.id_column)
        taxon_col = self._find_column(df.columns, TAXON_ID)
        if id_col is None or taxon_col is None:
            raise ValueError(f"Could not find vgnc_id/taxon_id columns in {file_path}")
        return {
            str(vid): str(taxon)
            for vid, taxon in df.select([id_col, taxon_col]).drop_nulls().rows()
        }

    def _extract_primary_labels(self, file_path: Path, species: str) -> dict[str, set[str]]:
        """Extract all current approved VGNC symbols for *species*.

        Returns a ``dict`` mapping each symbol text to the set of primary
        VGNC IDs (within *species*) that carry that symbol.

        Args:
            file_path: Path to the VGNC gene-set TSV file.
            species: NCBI taxon ID to filter by, or :data:`ALL_SPECIES` to
                process every species together.

        Returns:
            ``dict[label, set[vgnc_id]]``
        """
        df = self._read_tsv(file_path)
        symbol_col = self._find_column(df.columns, SYMBOL)
        id_col = self._find_column(df.columns, self.id_column)
        status_col = self._find_column(df.columns, STATUS)
        taxon_col = self._find_column(df.columns, TAXON_ID)
        if symbol_col is None or id_col is None or taxon_col is None:
            raise ValueError(f"Could not find required columns in {file_path}")

        filtered = (
            df
            if species == ALL_SPECIES
            else df.filter(pl.col(taxon_col).cast(pl.Utf8) == str(species))
        )
        if status_col:
            filtered = filtered.filter(pl.col(status_col) == "Approved")

        return self._labels_by_id(filtered, id_col, symbol_col)

    def _parse_gene_set(
        self,
        file_path: Path,
        species: str,
        statuses: list[str] | None = None,
    ) -> list[Mapping]:
        """Parse the VGNC gene-set file for symbol (label) mappings, scoped to *species*.

        Args:
            file_path: Path to the VGNC gene-set TSV file.
            species: NCBI taxon ID to filter by, or :data:`ALL_SPECIES` to
                process every species together.
            statuses: Entry statuses to include (e.g. ``["Approved"]``).
                If ``None`` (default), all entries are included.

        Returns:
            List of SSSOM Mapping objects for label mappings.
        """
        df = self._read_tsv(file_path)

        taxon_col = self._find_column(df.columns, TAXON_ID)
        status_col = self._find_column(df.columns, STATUS)
        id_col = self._find_column(df.columns, self.id_column)
        label_col = self._find_column(df.columns, SYMBOL)
        alias_col = self._find_column(df.columns, ALIAS_SYMBOL)
        prev_col = self._find_column(df.columns, PREV_SYMBOL)
        date_changed_col = self._find_column(df.columns, DATE_SYMBOL_CHANGED)

        if not all([taxon_col, id_col, label_col]):
            raise ValueError(f"Missing required columns in {file_path}")
        assert id_col is not None
        assert label_col is not None
        assert taxon_col is not None

        df_species = (
            df
            if species == ALL_SPECIES
            else df.filter(pl.col(taxon_col).cast(pl.Utf8) == str(species))
        )
        if statuses is not None and status_col:
            df_species = df_species.filter(pl.col(status_col).is_in(statuses))

        return self._build_label_mappings(
            df_species,
            id_col=id_col,
            label_col=label_col,
            alias_col=alias_col,
            prev_col=prev_col,
            date_changed_col=date_changed_col,
            taxon_col=taxon_col,
        )


__all__ = ["VGNCParser"]
