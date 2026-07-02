"""VGNC TSV file parser for secondary-to-primary identifier mappings.

VGNC (Vertebrate Gene Nomenclature Committee) names genes for vertebrate
species that lack their own nomenclature committee. Structurally it
publishes the same two-file shape as HGNC (a withdrawn file and a complete
gene-set file), but as a *single global* TSV per file covering every VGNC
species at once, with a ``taxon_id`` column -- closer to how NCBI's
``gene_info``/``gene_history`` are filtered by ``#tax_id`` than to Ensembl's
one-URL-per-species layout.

This parser extracts:
1. ID-to-ID mappings: withdrawn/merged VGNC IDs -> current VGNC IDs. The
   withdrawn file's own ``taxon_id`` column is not populated upstream, so
   :meth:`VGNCParser.parse` always parses the *full*, unfiltered withdrawn
   file first; an optional ``species`` then *subsets the output* by
   resolving each mapping's primary (replacement) VGNC ID against the
   gene-set file's ``taxon_id`` column. Withdrawn entries with no
   resolvable replacement can't be attributed to a species and are
   dropped when subsetting (they could belong to any one).
2. Label-to-label mappings: previous/alias gene symbols -> current symbols,
   from the gene-set file. Symbols are *not* unique across species (the same
   approved symbol can legitimately name orthologous genes in different
   species), so every label-mapping method requires an explicit ``species``
   (NCBI taxon ID) filter -- this is what lets cardinality/ambiguity
   detection and ``to_pri_labels()`` stay scoped to one species' namespace
   instead of flagging cross-species homonyms as ambiguous. Passing
   :data:`ALL_SPECIES` processes every species together instead, in which
   case a shared symbol is really ambiguous (there's no other context
   to tell the species apart) and gets flagged accordingly.

Uses SSSOM-compliant MappingSet classes with cardinality computation.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
from sssom_schema import Mapping

from pysec2pri.parsers.base import (
    WITHDRAWN_ENTRY,
    WITHDRAWN_ENTRY_LABEL,
    BaseMappingSet,
    BaseParser,
    LabelMappingSet,
)

# VGNC column names (case-insensitive matching used)
VGNC_ID = "vgnc_id"
SYMBOL = "symbol"
WITHDRAWN_SYMBOL = "withdrawn_symbol"
ALIAS_SYMBOL = "alias_symbol"
PREV_SYMBOL = "prev_symbol"
DATE_SYMBOL_CHANGED = "date_symbol_changed"
STATUS = "status"
TAXON_ID = "taxon_id"

#: Sentinel accepted by every species-aware method: skip taxon filtering/
#: subsetting entirely and process every species together.
ALL_SPECIES = "all"


class VGNCParser(BaseParser):
    """Parser for VGNC TSV files using Polars for memory efficiency.

    Extracts secondary-to-primary VGNC identifier mappings and symbol
    mappings from the VGNC withdrawn and gene-set files.

    Returns:
    - IdMappingSet for ID-to-ID mappings (withdrawn/merged IDs)
    - LabelMappingSet for symbol mappings (alias/previous symbols),
      scoped to one species
    """

    datasource_name = "vgnc"

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

    def _product_slug(self) -> str | None:
        """Species (NCBI taxon ID) a label-mapping run was filtered to.

        ``None`` for ID mappings (the withdrawn file is never species
        filtered, see module docstring), so only label-mapping runs fold a
        taxon ID into ``mapping_set_id``/``record_id`` (see
        :meth:`~pysec2pri.parsers.base.BaseParser._product_slug`).
        """
        return getattr(self, "species", None)

    @property
    def withdrawn_source_url(self) -> str:
        """Get the withdrawn file download URL from config."""
        return self.get_download_url("withdrawn") or ""

    @property
    def gene_set_source_url(self) -> str:
        """Get the gene-set download URL from config."""
        return self.get_download_url("complete") or ""

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

        mappings = self._parse_withdrawn(input_path)

        taxon_by_id: dict[str, str] | None = None
        if complete_set_path is not None:
            taxon_by_id = self._taxon_by_vgnc_id(Path(complete_set_path))

        if species not in (None, ALL_SPECIES) and taxon_by_id is not None:
            mappings = [
                m
                for m in mappings
                if taxon_by_id.get(str(getattr(m, "object_id", "") or "")) == str(species)
            ]

        mapping_set = self._create_mapping_set(mappings, mapping_type="id")

        if complete_set_path is not None:
            primary_ids = self._extract_primary_ids(Path(complete_set_path))
            if species not in (None, ALL_SPECIES) and taxon_by_id is not None:
                primary_ids = {vid for vid in primary_ids if taxon_by_id.get(vid) == str(species)}
            object.__setattr__(mapping_set, "_primary_ids", primary_ids)
        return mapping_set

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

        mapping_set = self._create_mapping_set([], mapping_type="id")
        object.__setattr__(mapping_set, "_primary_ids", primary_ids)
        return mapping_set

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
        mapping_set = self._create_mapping_set(mappings, mapping_type="label")
        object.__setattr__(
            mapping_set,
            "_primary_labels",
            self._extract_primary_labels(complete_set_path, species),
        )
        return mapping_set

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

        mapping_set = self._create_mapping_set([], mapping_type="label")
        object.__setattr__(
            mapping_set,
            "_primary_labels",
            self._extract_primary_labels(complete_set_path, species),
        )
        return mapping_set

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

    def _extract_primary_ids(self, file_path: Path) -> set[str]:
        """Extract all current VGNC IDs from the gene-set file (every species).

        Args:
            file_path: Path to the VGNC gene-set TSV file.

        Returns:
            Set of all VGNC IDs present in the gene-set file.
        """
        df = pl.read_csv(
            file_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=[""],
        )
        vgnc_id_col = self._find_column(df.columns, VGNC_ID)
        if vgnc_id_col is None:
            raise ValueError(f"Could not find vgnc_id column in {file_path}")
        return {str(val) for val in df[vgnc_id_col].drop_nulls().to_list()}

    def _taxon_by_vgnc_id(self, file_path: Path) -> dict[str, str]:
        """Return ``{vgnc_id: taxon_id}`` for every current gene in the gene-set file.

        Used to resolve which species a withdrawn entry's *replacement*
        gene belongs to, since the withdrawn file's own ``taxon_id``
        column is not populated upstream (see module docstring).

        Args:
            file_path: Path to the VGNC gene-set TSV file.

        Returns:
            ``dict[vgnc_id, taxon_id]``.
        """
        df = pl.read_csv(
            file_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=[""],
        )
        vgnc_id_col = self._find_column(df.columns, VGNC_ID)
        taxon_col = self._find_column(df.columns, TAXON_ID)
        if vgnc_id_col is None or taxon_col is None:
            raise ValueError(f"Could not find vgnc_id/taxon_id columns in {file_path}")
        return {
            str(vid): str(taxon)
            for vid, taxon in df.select([vgnc_id_col, taxon_col]).drop_nulls().rows()
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
        df = pl.read_csv(
            file_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=[""],
        )
        symbol_col = self._find_column(df.columns, SYMBOL)
        vgnc_id_col = self._find_column(df.columns, VGNC_ID)
        status_col = self._find_column(df.columns, STATUS)
        taxon_col = self._find_column(df.columns, TAXON_ID)
        if symbol_col is None or vgnc_id_col is None or taxon_col is None:
            raise ValueError(f"Could not find required columns in {file_path}")

        filtered = (
            df
            if species == ALL_SPECIES
            else df.filter(pl.col(taxon_col).cast(pl.Utf8) == str(species))
        )
        if status_col:
            filtered = filtered.filter(pl.col(status_col) == "Approved")

        result: dict[str, set[str]] = {}
        for id_, label in filtered.select([vgnc_id_col, symbol_col]).drop_nulls().rows():
            result.setdefault(str(label), set()).add(str(id_))
        return result

    def _parse_withdrawn(self, file_path: Path) -> list[Mapping]:
        """Parse the VGNC withdrawn file for ID-to-ID mappings.

        Args:
            file_path: Path to the VGNC withdrawn TSV file.

        Returns:
            List of SSSOM Mapping objects.
        """
        df = pl.read_csv(
            file_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=[""],
        )

        merged_col = self._find_merged_column(df.columns, [])
        if merged_col is None:
            raise ValueError(f"Could not find merged_into_report column in {file_path}")

        vgnc_id_col = self._find_column(df.columns, VGNC_ID)
        if vgnc_id_col is None:
            raise ValueError(f"Could not find vgnc_id column in {file_path}")

        status_col = self._find_column(df.columns, STATUS)
        label_col = self._find_column(df.columns, WITHDRAWN_SYMBOL)

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
            vgnc_id = row.get(vgnc_id_col)
            if not vgnc_id:
                continue

            merged_info = row.get(merged_col)
            status = row.get(status_col) if status_col else None
            label = row.get(label_col) if label_col else None

            # Case 1: Withdrawn with no replacement
            if not merged_info and status and "Entry Withdrawn" in str(status):
                rows_data.append(
                    {
                        "subject_id": vgnc_id,
                        "object_id": WITHDRAWN_ENTRY,
                        "subject_label": label or "",
                        "object_label": WITHDRAWN_ENTRY_LABEL,
                        "predicate_id": "oboInOwl:consider",
                        "comment": "Withdrawn entry with no replacement.",
                        "record_id": self._record_id(
                            self._record_namespace(),
                            WITHDRAWN_ENTRY,
                            vgnc_id,
                        ),
                    }
                )
                continue

            # Case 2: Merged into another entry
            if merged_info:
                parsed = self._parse_merged_info(merged_info)
                if parsed:
                    target_id, target_label = parsed
                    rows_data.append(
                        {
                            "subject_id": vgnc_id,
                            "object_id": target_id,
                            "subject_label": label or "",
                            "object_label": target_label or "",
                            "predicate_id": m_meta["predicate_id"],
                            "predicate_label": m_meta.get("predicate_label"),
                            "record_id": self._record_id(
                                self._record_namespace(),
                                target_id,
                                vgnc_id,
                            ),
                        }
                    )

        return self._build_mappings(
            rows_data, fixed, desc="Processing withdrawn", total=len(rows_data)
        )

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
        df = pl.read_csv(
            file_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=[""],
        )

        taxon_col = self._find_column(df.columns, TAXON_ID)
        status_col = self._find_column(df.columns, STATUS)
        vgnc_id_col = self._find_column(df.columns, VGNC_ID)
        label_col = self._find_column(df.columns, SYMBOL)
        alias_col = self._find_column(df.columns, ALIAS_SYMBOL)
        prev_col = self._find_column(df.columns, PREV_SYMBOL)
        date_changed_col = self._find_column(df.columns, DATE_SYMBOL_CHANGED)

        if not all([taxon_col, vgnc_id_col, label_col]):
            raise ValueError(f"Missing required columns in {file_path}")
        assert vgnc_id_col is not None
        assert label_col is not None
        assert taxon_col is not None

        df_species = (
            df
            if species == ALL_SPECIES
            else df.filter(pl.col(taxon_col).cast(pl.Utf8) == str(species))
        )
        if statuses is not None and status_col:
            df_species = df_species.filter(pl.col(status_col).is_in(statuses))

        m_meta = self.get_mapping_metadata()
        fixed = {
            "mapping_justification": m_meta["mapping_justification"],
            "subject_source": m_meta.get("subject_source"),
            "object_source": m_meta.get("object_source"),
            "mapping_tool": m_meta.get("mapping_tool"),
            "license": m_meta.get("license"),
        }

        rows_data: list[dict[str, str | None]] = []
        for row in df_species.iter_rows(named=True):
            vgnc_id = row.get(vgnc_id_col)
            label = row.get(label_col)
            if not vgnc_id or not label:
                continue

            alias_str = row.get(alias_col) if alias_col else None
            prev_str = row.get(prev_col) if prev_col else None
            aliases = self._split_labels(labels_str=alias_str) if alias_str else []
            prev_labels = self._split_labels(labels_str=prev_str) if prev_str else []
            # The date the current symbol was set, i.e. when its previous
            # symbol(s) became secondary. VGNC records only the most recent
            # change, so with multiple prev_symbol entries this date applies
            # exactly to the latest rename and approximately to earlier ones.
            symbol_changed_date = row.get(date_changed_col) if date_changed_col else None

            for alias in aliases:
                rows_data.append(
                    {
                        "object_id": vgnc_id,
                        "subject_label": alias,
                        "subject_type": "rdfs literal",
                        "object_label": label,
                        "_label_type": "alias",
                        "comment": "Alias symbol mapping.",
                        "record_id": self._record_id(
                            self._record_namespace(),
                            vgnc_id,
                            alias,
                        ),
                    }
                )

            for prev in prev_labels:
                rows_data.append(
                    {
                        "object_id": vgnc_id,
                        "subject_label": prev,
                        "subject_type": "rdfs literal",
                        "object_label": label,
                        "_label_type": "previous",
                        "comment": "Previous symbol mapping.",
                        "mapping_date": symbol_changed_date,
                        "record_id": self._record_id(
                            self._record_namespace(),
                            vgnc_id,
                            prev,
                        ),
                    }
                )

        return self._build_mappings(
            rows_data, fixed, desc="Processing symbols", total=len(rows_data)
        )

    def _create_mapping_set(
        self, mappings: list[Mapping], mapping_type: str = "id"
    ) -> BaseMappingSet:
        """Create an IdMappingSet or LabelMappingSet with config metadata.

        Delegates to BaseParser.create_mapping_set().
        """
        return self.create_mapping_set(mappings, mapping_type)


__all__ = ["VGNCParser"]
