"""Ensembl core flat-file parser for secondary-to-primary gene identifier mappings.

Reads the MySQL flat-file dumps Ensembl publishes per release/species:

- ``stable_id_event.txt.gz``: cumulative old (secondary) -> new (primary)
  stable-ID history. Each release's table already holds the *entire* history
  up to that release, so one release's mappings describe the full state of
  Ensembl at that point in time (no cross-release chain resolution needed
  here -- :mod:`pysec2pri.consolidate` handles cross-release temporality).
- ``mapping_session.txt.gz``: maps each ``mapping_session_id`` to the
  ``created`` timestamp used for the per-row ``mapping_date``.
- ``gene.txt.gz`` / ``xref.txt.gz`` / ``external_synonym.txt.gz``: the
  current gene set, its display labels, and symbol synonyms.

This parser extracts:
1. ID-to-ID mappings: retired/old Ensembl gene IDs -> current Ensembl gene IDs
2. Label-to-label mappings: gene symbol synonyms -> current display label

Uses SSSOM-compliant MappingSet classes with cardinality computation.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import httpx
import polars as pl
from sssom_schema import Mapping

from pysec2pri.parsers.base import (
    WITHDRAWN_ENTRY,
    WITHDRAWN_ENTRY_LABEL,
    BaseDownloader,
    BaseParser,
    Sec2PriMappingSet,
)

# Column layouts verified against the Ensembl release 115 core schema

_STABLE_ID_EVENT_COLUMNS = [
    "old_stable_id",
    "old_version",
    "new_stable_id",
    "new_version",
    "mapping_session_id",
    "type",
    "score",
]
_MAPPING_SESSION_COLUMNS = [
    "mapping_session_id",
    "old_db_name",
    "new_db_name",
    "old_release",
    "new_release",
    "old_assembly",
    "new_assembly",
    "created",
]
_GENE_COLUMNS = [
    "gene_id",
    "biotype",
    "analysis_id",
    "seq_region_id",
    "seq_region_start",
    "seq_region_end",
    "seq_region_strand",
    "display_xref_id",
    "source",
    "description",
    "is_current",
    "canonical_transcript_id",
    "stable_id",
    "version",
    "created_date",
    "modified_date",
]
_XREF_COLUMNS = [
    "xref_id",
    "external_db_id",
    "dbprimary_acc",
    "display_label",
    "version",
    "description",
    "info_type",
    "info_text",
]
_EXTERNAL_SYNONYM_COLUMNS = ["xref_id", "synonym"]

# Explicit dtypes where the natural column type (e.g. stable_id_event.score,
# mostly 0 with rare floats) would otherwise be mis-inferred from a sample.
_DType = type[pl.DataType]
_STABLE_ID_EVENT_DTYPES: dict[str, _DType] = {
    "old_version": pl.Int64,
    "new_version": pl.Int64,
    "mapping_session_id": pl.Int64,
    "score": pl.Float64,
}
_GENE_DTYPES: dict[str, _DType] = {
    "display_xref_id": pl.Int64,
    "is_current": pl.Int64,
}
_XREF_DTYPES: dict[str, _DType] = {"xref_id": pl.Int64}
_EXTERNAL_SYNONYM_DTYPES: dict[str, _DType] = {"xref_id": pl.Int64}


def _scan_ensembl_tsv(
    path: Path, columns: list[str], schema_overrides: dict[str, _DType] | None = None
) -> pl.LazyFrame:
    """Lazily scan a header-less Ensembl flat-file dump."""
    return pl.scan_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=columns,
        null_values=["\\N"],
        quote_char=None,
        infer_schema_length=10000,
        schema_overrides=schema_overrides,
    )


_ENSEMBL_REST_SPECIES_URL = "https://rest.ensembl.org/info/species?content-type=application/json"


def _lookup_ensembl_species_token(taxon_id: str | int) -> str | None:
    """Resolve a taxon ID to its Ensembl species token via Ensembl's own REST API.

    ``config/ensembl.yaml``'s ``species.available`` only curates a handful
    of species for CLI defaults/help; Ensembl itself publishes 300+ (e.g.
    dog, pig, chicken). This is the live fallback for anything not in that
    static list, queried on demand rather than hardcoded, since the list
    changes every release.

    Several taxon IDs map to more than one entry (per-breed/per-strain
    assemblies, e.g. multiple dog breeds): the *shortest* matching name is
    returned, which is reliably the canonical species-level entry (breed/
    strain variants always have additional characters appended).

    Args:
        taxon_id: Canonical NCBI taxon ID.

    Returns:
        The Ensembl species token (e.g. ``"canis_lupus_familiaris"``), or
        ``None`` if the REST API is unreachable or has no matching entry.
    """
    try:
        with httpx.Client(follow_redirects=True, timeout=30.0) as client:
            response = client.get(_ENSEMBL_REST_SPECIES_URL)
            response.raise_for_status()
            data = response.json()
    except httpx.HTTPError:
        return None

    matches = [
        str(entry["name"])
        for entry in data.get("species", [])
        if entry.get("name") and str(entry.get("taxon_id")) == str(taxon_id)
    ]
    return min(matches, key=len) if matches else None


def _ensembl_date_to_iso(value: object) -> str | None:
    """Convert a ``mapping_session.created`` timestamp to ISO ``YYYY-MM-DD``.

    Args:
        value: Raw ``created`` cell, e.g. ``"2002-09-03 11:19:48"``.

    Returns:
        ISO-8601 date string, or ``None`` if *value* is missing/malformed.
    """
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    return text.split(" ", 1)[0]


class EnsemblParser(BaseParser):
    """Parser for Ensembl core flat-file dumps using Polars.

    Each release's ``stable_id_event`` table is cumulative, so a single
    parse describes the *whole state* of Ensembl gene IDs at that release
    (no chain-walking across releases).

    Returns:
    - IdMappingSet for ID-to-ID mappings (retired Ensembl gene IDs)
    - LabelMappingSet for symbol mappings (external gene synonyms)
    """

    datasource_name = "ensembl"

    def __init__(
        self,
        version: str | None = None,
        show_progress: bool = True,
        species: str | int = 9606,
    ) -> None:
        """Initialize the Ensembl parser.

        Args:
            version: Ensembl release number (e.g. ``"115"``).
            show_progress: Whether to show progress bars during parsing.
            species: Canonical NCBI taxon ID of the species being parsed.
                Unlike NCBI, this never filters rows (the input files are
                already species-specific downloads) -- it only labels the
                output (see :meth:`_product_slug`).
        """
        super().__init__(version=version, show_progress=show_progress)
        self.species = species

    def _product_slug(self) -> str | None:
        """Species (NCBI taxon ID) the parsed files belong to.

        Different species are disjoint datasets at the same release, so the
        taxon ID is folded into ``mapping_set_id``/``record_id`` (see
        :meth:`~pysec2pri.parsers.base.BaseParser._product_slug`). When this
        instance built its mappings via :meth:`parse_label_history`, a
        trailing ``consolidate`` segment marks the IRI as a derived,
        cross-release product rather than a single release's parse --
        ``version`` is then the *last* release that walk consolidated.
        """
        species_slug = str(self.species) if self.species is not None else None
        if getattr(self, "_is_label_history", False):
            return f"{species_slug}/consolidate" if species_slug else "consolidate"
        return species_slug

    def parse(
        self,
        stable_id_event_path: Path | str | None = None,
        mapping_session_path: Path | str | None = None,
        gene_path: Path | str | None = None,
    ) -> Sec2PriMappingSet:
        """Parse ``stable_id_event`` (+ ``mapping_session``) into an IdMappingSet.

        Args:
            stable_id_event_path: Path to ``stable_id_event.txt`` (can be
                ``.gz`` compressed).
            mapping_session_path: Optional path to ``mapping_session.txt``,
                used to resolve each row's ``mapping_date``. When omitted,
                rows fall back to the set-level release date.
            gene_path: Optional path to ``gene.txt``. When supplied,
                ``_primary_ids`` is populated with every current gene ID.

        Returns:
            IdMappingSet with computed cardinalities based on IDs.
        """
        if stable_id_event_path is None:
            raise ValueError("stable_id_event_path must not be None")
        stable_id_event_path = Path(stable_id_event_path)
        self._resolve_version(stable_id_event_path)

        mappings = self._parse_stable_id_event(stable_id_event_path, mapping_session_path)
        mapping_set = self.create_mapping_set(mappings, mapping_type="id")

        if gene_path is not None:
            object.__setattr__(
                mapping_set, "_primary_ids", self._extract_primary_ids(Path(gene_path))
            )

        return mapping_set

    def parse_labels(
        self,
        gene_path: Path | str | None = None,
        xref_path: Path | str | None = None,
        external_synonym_path: Path | str | None = None,
    ) -> Sec2PriMappingSet:
        """Parse external gene synonyms into a LabelMappingSet.

        Args:
            gene_path: Path to ``gene.txt``.
            xref_path: Path to ``xref.txt``.
            external_synonym_path: Optional path to ``external_synonym.txt``.
                When omitted, the returned set carries only the full primary
                label set (no synonym mappings).

        Returns:
            LabelMappingSet with computed cardinalities based on labels.
        """
        if gene_path is None or xref_path is None:
            raise ValueError("gene_path and xref_path must not be None")
        gene_path = Path(gene_path)
        xref_path = Path(xref_path)
        self._resolve_version(gene_path)

        mappings: list[Mapping] = []
        if external_synonym_path is not None:
            mappings = self._parse_external_synonyms(
                gene_path, xref_path, Path(external_synonym_path)
            )

        mapping_set = self.create_mapping_set(mappings, mapping_type="label")
        object.__setattr__(
            mapping_set, "_primary_labels", self._extract_primary_labels(gene_path, xref_path)
        )
        object.__setattr__(mapping_set, "_primary_ids", self._extract_primary_ids(gene_path))
        return mapping_set

    def parse_primary_ids(self, gene_path: Path | str | None = None) -> Sec2PriMappingSet:
        """Return a mapping set containing the full list of current Ensembl gene IDs.

        Args:
            gene_path: Path to ``gene.txt`` (can be ``.gz`` compressed).

        Returns:
            :class:`~pysec2pri.parsers.base.IdMappingSet` with ``_primary_ids``
            populated.
        """
        if gene_path is None:
            raise ValueError("gene_path must not be None")
        gene_path = Path(gene_path)
        self._resolve_version(gene_path)
        ms = self.create_mapping_set([], mapping_type="id")
        object.__setattr__(ms, "_primary_ids", self._extract_primary_ids(gene_path))
        return ms

    def parse_primary_labels(
        self,
        gene_path: Path | str | None = None,
        xref_path: Path | str | None = None,
    ) -> Sec2PriMappingSet:
        """Return a mapping set containing the full list of current Ensembl gene labels.

        Args:
            gene_path: Path to ``gene.txt``.
            xref_path: Path to ``xref.txt``.

        Returns:
            :class:`~pysec2pri.parsers.base.LabelMappingSet` with
            ``_primary_labels`` populated.
        """
        if gene_path is None or xref_path is None:
            raise ValueError("gene_path and xref_path must not be None")
        gene_path = Path(gene_path)
        xref_path = Path(xref_path)
        self._resolve_version(gene_path)
        ms = self.create_mapping_set([], mapping_type="label")
        object.__setattr__(
            ms, "_primary_labels", self._extract_primary_labels(gene_path, xref_path)
        )
        return ms

    def parse_all(
        self,
        stable_id_event_path: Path | str | None,
        mapping_session_path: Path | str | None,
        gene_path: Path | str | None,
        xref_path: Path | str | None,
        external_synonym_path: Path | str | None,
    ) -> tuple[Sec2PriMappingSet, Sec2PriMappingSet]:
        """Parse the full set of Ensembl core flat files.

        Args:
            stable_id_event_path: Path to ``stable_id_event.txt``.
            mapping_session_path: Path to ``mapping_session.txt``.
            gene_path: Path to ``gene.txt``.
            xref_path: Path to ``xref.txt``.
            external_synonym_path: Path to ``external_synonym.txt``.

        Returns:
            Tuple of (IdMappingSet, LabelMappingSet).
        """
        id_mappings = self.parse(stable_id_event_path, mapping_session_path, gene_path=gene_path)
        label_mappings = self.parse_labels(gene_path, xref_path, external_synonym_path)
        return id_mappings, label_mappings

    def _parse_stable_id_event(
        self,
        file_path: Path,
        mapping_session_path: Path | str | None,
    ) -> list[Mapping]:
        """Parse ``stable_id_event`` for gene ID-to-ID mappings.

        Args:
            file_path: Path to ``stable_id_event.txt``.
            mapping_session_path: Optional path to ``mapping_session.txt``.

        Returns:
            List of SSSOM Mapping objects.
        """
        try:
            df = (
                _scan_ensembl_tsv(file_path, _STABLE_ID_EVENT_COLUMNS, _STABLE_ID_EVENT_DTYPES)
                .filter(pl.col("type") == "gene")
                .collect()
            )
        except pl.exceptions.NoDataError:
            return []
        if df.is_empty():
            return []

        date_by_session: dict[str, str | None] = {}
        if mapping_session_path is not None:
            sessions = _scan_ensembl_tsv(
                Path(mapping_session_path), _MAPPING_SESSION_COLUMNS
            ).collect()
            date_by_session = {
                str(row["mapping_session_id"]): _ensembl_date_to_iso(row["created"])
                for row in sessions.iter_rows(named=True)
            }

        m_meta = self.get_mapping_metadata()
        fixed_base = {
            "mapping_justification": m_meta["mapping_justification"],
            "subject_source": m_meta.get("subject_source"),
            "object_source": m_meta.get("object_source"),
            "mapping_tool": m_meta.get("mapping_tool"),
            "license": m_meta.get("license"),
        }

        record_ns = self._record_namespace()
        rows_data: list[dict[str, Any]] = []
        for row in df.iter_rows(named=True):
            old_id = row.get("old_stable_id")
            new_id = row.get("new_stable_id")

            # No old (secondary) ID: a brand-new gene, not a sec->pri mapping.
            # old == new: an annotation-only version bump, not an ID change.
            if not old_id or old_id == new_id:
                continue

            session_id = str(row.get("mapping_session_id") or "")
            mapping_date = date_by_session.get(session_id)
            score = row.get("score")

            if not new_id:
                rows_data.append(
                    {
                        "subject_id": f"ENSEMBL:{old_id}",
                        "object_id": WITHDRAWN_ENTRY,
                        "object_label": WITHDRAWN_ENTRY_LABEL,
                        "predicate_id": "oboInOwl:consider",
                        "mapping_date": mapping_date,
                        "record_id": self._record_id(record_ns, WITHDRAWN_ENTRY, old_id),
                    }
                )
            else:
                row_fields: dict[str, Any] = {
                    "subject_id": f"ENSEMBL:{old_id}",
                    "object_id": f"ENSEMBL:{new_id}",
                    "predicate_id": m_meta["predicate_id"],
                    "predicate_label": m_meta.get("predicate_label"),
                    "mapping_date": mapping_date,
                    "record_id": self._record_id(record_ns, new_id, old_id),
                }
                if score and float(score) > 0:
                    row_fields["confidence"] = float(score)
                rows_data.append(row_fields)

        return self._build_mappings(
            rows_data, fixed_base, desc="Processing stable_id_event", total=len(rows_data)
        )

    def _parse_external_synonyms(
        self,
        gene_path: Path,
        xref_path: Path,
        external_synonym_path: Path,
    ) -> list[Mapping]:
        """Parse ``external_synonym`` joined through genes' display xref.

        Only synonyms attached to the *same* xref a gene displays as its
        canonical label are picked up (there is no ``object_xref`` table in
        the required input list to resolve arbitrary gene<->xref links).

        Args:
            gene_path: Path to ``gene.txt``.
            xref_path: Path to ``xref.txt``.
            external_synonym_path: Path to ``external_synonym.txt``.

        Returns:
            List of SSSOM Mapping objects for synonym mappings.
        """
        genes = (
            _scan_ensembl_tsv(gene_path, _GENE_COLUMNS, _GENE_DTYPES)
            .filter(pl.col("is_current") == 1)
            .select(["display_xref_id", "stable_id"])
        )
        xrefs = _scan_ensembl_tsv(xref_path, _XREF_COLUMNS, _XREF_DTYPES).select(
            ["xref_id", "display_label"]
        )
        synonyms = _scan_ensembl_tsv(
            external_synonym_path, _EXTERNAL_SYNONYM_COLUMNS, _EXTERNAL_SYNONYM_DTYPES
        )

        df = (
            genes.join(xrefs, left_on="display_xref_id", right_on="xref_id", how="inner")
            .join(synonyms, left_on="display_xref_id", right_on="xref_id", how="inner")
            .select(["stable_id", "display_label", "synonym"])
            .drop_nulls()
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

        record_ns = self._record_namespace()
        rows_data: list[dict[str, Any]] = []
        for stable_id, display_label, synonym in df.rows():
            curie_id = f"ENSEMBL:{stable_id}"
            rows_data.append(
                {
                    "object_id": curie_id,
                    "subject_label": synonym,
                    "subject_type": "rdfs literal",
                    "object_label": display_label,
                    "_label_type": "alias",
                    "comment": "Ensembl gene external synonym.",
                    "record_id": self._record_id(record_ns, curie_id, synonym),
                }
            )

        return self._build_mappings(
            rows_data, fixed, desc="Processing external_synonym", total=len(rows_data)
        )

    def _extract_primary_ids(self, file_path: Path) -> set[str]:
        """Extract all current Ensembl gene IDs from ``gene.txt``.

        Args:
            file_path: Path to ``gene.txt``.

        Returns:
            Set of ``ENSEMBL:<stable_id>`` CURIEs.
        """
        df = (
            _scan_ensembl_tsv(file_path, _GENE_COLUMNS, _GENE_DTYPES)
            .filter(pl.col("is_current") == 1)
            .select("stable_id")
            .collect()
        )
        return {f"ENSEMBL:{v}" for v in df["stable_id"].drop_nulls().to_list()}

    def _extract_primary_labels(self, gene_path: Path, xref_path: Path) -> dict[str, set[str]]:
        """Extract all current gene display labels from ``gene.txt`` + ``xref.txt``.

        Returns a ``dict`` mapping each label text to the set of primary
        ``ENSEMBL:<stable_id>`` IDs that carry it.

        Args:
            gene_path: Path to ``gene.txt``.
            xref_path: Path to ``xref.txt``.

        Returns:
            ``dict[label, set[ENSEMBL:<stable_id>]]``
        """
        genes = (
            _scan_ensembl_tsv(gene_path, _GENE_COLUMNS, _GENE_DTYPES)
            .filter(pl.col("is_current") == 1)
            .select(["display_xref_id", "stable_id"])
        )
        xrefs = _scan_ensembl_tsv(xref_path, _XREF_COLUMNS, _XREF_DTYPES).select(
            ["xref_id", "display_label"]
        )
        df = (
            genes.join(xrefs, left_on="display_xref_id", right_on="xref_id", how="inner")
            .select(["stable_id", "display_label"])
            .drop_nulls()
            .collect()
        )
        result: dict[str, set[str]] = {}
        for stable_id, label in df.rows():
            result.setdefault(str(label), set()).add(f"ENSEMBL:{stable_id}")
        return result

    def current_label_snapshot(
        self, gene_path: Path | str, xref_path: Path | str
    ) -> dict[str, str]:
        """Return ``{bare stable_id -> current display label}`` for every current gene.

        Unlike :meth:`_extract_primary_labels` (keyed by label, for ambiguity
        detection and ``to_pri_labels()``), this is keyed by gene -- the
        natural shape for diffing one release's label snapshot against the
        next (see :func:`pysec2pri.consolidate.build_label_history`).

        Args:
            gene_path: Path to ``gene.txt``.
            xref_path: Path to ``xref.txt``.

        Returns:
            ``dict[stable_id, label]`` (bare stable IDs, no ``ENSEMBL:`` prefix).
        """
        genes = (
            _scan_ensembl_tsv(Path(gene_path), _GENE_COLUMNS, _GENE_DTYPES)
            .filter(pl.col("is_current") == 1)
            .select(["display_xref_id", "stable_id"])
        )
        xrefs = _scan_ensembl_tsv(Path(xref_path), _XREF_COLUMNS, _XREF_DTYPES).select(
            ["xref_id", "display_label"]
        )
        df = (
            genes.join(xrefs, left_on="display_xref_id", right_on="xref_id", how="inner")
            .select(["stable_id", "display_label"])
            .drop_nulls()
            .collect()
        )
        return {str(stable_id): str(label) for stable_id, label in df.rows()}

    def parse_label_history(
        self, transitions: Iterable[tuple[str, str, str, str | None]]
    ) -> Sec2PriMappingSet:
        """Build a LabelMappingSet from precomputed previous->current label transitions.

        Args:
            transitions: Iterable of ``(stable_id, prev_label, curr_label,
                mapping_date)`` tuples, one per gene whose display label
                changed between two releases (see
                :func:`pysec2pri.consolidate.build_label_history`).

        Returns:
            LabelMappingSet with ``IAO:0100001`` ("term replaced by") mappings.
        """
        self._is_label_history = True
        m_meta = self.get_mapping_metadata()
        fixed = {
            "mapping_justification": m_meta["mapping_justification"],
            "subject_source": m_meta.get("subject_source"),
            "object_source": m_meta.get("object_source"),
            "mapping_tool": m_meta.get("mapping_tool"),
            "license": m_meta.get("license"),
        }

        record_ns = self._record_namespace()
        rows_data: list[dict[str, Any]] = []
        for stable_id, prev_label, curr_label, mapping_date in transitions:
            curie_id = f"ENSEMBL:{stable_id}"
            rows_data.append(
                {
                    "object_id": curie_id,
                    "subject_label": prev_label,
                    "subject_type": "rdfs literal",
                    "object_label": curr_label,
                    "_label_type": "previous",
                    "mapping_date": mapping_date,
                    "comment": "Derived from snapshots across releases.",
                    "record_id": self._record_id(record_ns, curie_id, prev_label),
                }
            )

        mappings = self._build_mappings(
            rows_data, fixed, desc="Building label history", total=len(rows_data)
        )
        return self.create_mapping_set(mappings, mapping_type="label")


class EnsemblDownloader(BaseDownloader):
    """Downloader for Ensembl core flat-file dumps.

    Download paths are templated by release, species, and assembly suffix
    (e.g. ``homo_sapiens_core_115_38``); the assembly suffix is
    auto-discovered per release via :meth:`_resolve_core_dir`.
    """

    datasource_name = "ensembl"

    def __init__(
        self,
        version: str | None = None,
        species: str | int = 9606,
        assembly: str | None = None,
        show_progress: bool = True,
    ) -> None:
        """Initialize the Ensembl downloader.

        Args:
            version: Ensembl release number (e.g. ``"115"``).
            species: Canonical NCBI taxon ID, resolved to Ensembl's own
                species token via ``config.species_token()``.
            assembly: Force a specific assembly suffix (e.g. ``"38"``),
                skipping auto-discovery.
            show_progress: Whether to show progress bars.
        """
        super().__init__(version=version, show_progress=show_progress)
        self.species = species
        self.assembly = assembly

    def _resolve_core_dir(self, version: str, token: str) -> str:
        r"""Discover the per-species assembly suffix for *version*.

        Lists ``.../release-{version}/mysql/`` and matches
        ``{token}_core_{version}_(\\d+)``; falls back to
        ``parse_options.default_assembly`` when the listing can't be
        fetched or no match is found.

        Args:
            version: Ensembl release number.
            token: Species token (e.g. ``"homo_sapiens"``).

        Returns:
            Assembly suffix string (e.g. ``"38"``).
        """
        default_assembly = "38"
        if self._config:
            default_assembly = str(
                self._config.parse_options.get("default_assembly", default_assembly)
            )

        index_url = f"https://ftp.ensembl.org/pub/release-{version}/mysql/"
        try:
            with httpx.Client(follow_redirects=True, timeout=30.0) as client:
                response = client.get(index_url)
                response.raise_for_status()
        except httpx.HTTPError:
            return default_assembly

        match = re.search(
            rf"{re.escape(token)}_core_{re.escape(str(version))}_(\d+)", response.text
        )
        return match.group(1) if match else default_assembly

    def _resolve_species_token(self, taxon_id: str | int) -> str:
        """Resolve *taxon_id* to its Ensembl species token.

        Tries ``config/ensembl.yaml``'s static ``species.available`` map
        first (no network); falls back to a live lookup against Ensembl's
        REST species list (see :func:`_lookup_ensembl_species_token`) for
        taxon IDs not curated there -- Ensembl publishes 300+ species, of
        which ``available`` only lists a handful for CLI defaults/help.

        Args:
            taxon_id: Canonical NCBI taxon ID.

        Returns:
            The Ensembl species token.

        Raises:
            ValueError: If *taxon_id* is unknown both statically and live.
        """
        if self._config:
            try:
                return self._config.species_token(taxon_id)
            except ValueError:
                pass

        token = _lookup_ensembl_species_token(taxon_id)
        if token is not None:
            return token

        known = ""
        if self._config:
            available = (self._config.species or {}).get("available") or {}
            known = ", ".join(sorted(str(k) for k in available))
        raise ValueError(
            f"Unknown species taxon ID {taxon_id!r} for Ensembl: not in config/ensembl.yaml's "
            f"static list ({known or '(none configured)'}) and not found via Ensembl's live "
            "species list either."
        )

    def get_download_urls(
        self,
        version: str | None = None,
        **kwargs: Any,
    ) -> dict[str, str]:
        """Get Ensembl download URLs for a release/species.

        Args:
            version: Ensembl release number.
            **kwargs: Optional ``species``/``assembly`` overrides.

        Returns:
            Dictionary with file URLs keyed by table name.
        """
        v = version or self.version
        if v is None:
            raise ValueError("Ensembl requires an explicit release version.")
        if not self._config:
            raise ValueError("Ensembl config not loaded")

        species = kwargs.get("species", self.species)
        token = self._resolve_species_token(species)
        assembly = kwargs.get("assembly") or self.assembly or self._resolve_core_dir(str(v), token)

        return {
            key: url.format(version=v, species=token, assembly=assembly)
            for key, url in self._config.download_urls.items()
        }

    def download(
        self,
        output_dir: Path,
        version: str | None = None,
        decompress: bool = True,
        keys: list[str] | None = None,
        **kwargs: Any,
    ) -> dict[str, Path]:
        """Download Ensembl core flat files.

        Args:
            output_dir: Directory to save files.
            version: Ensembl release number.
            decompress: Whether to decompress ``.gz`` files.
            keys: Optional list of file-key names to download (defaults to
                all keys in ``ensembl.yaml``'s ``download_urls``).
            **kwargs: Forwarded to :meth:`get_download_urls`.

        Returns:
            Dictionary mapping file keys to downloaded paths.
        """
        urls = self.get_download_urls(version, **kwargs)
        if keys is not None:
            urls = {k: v for k, v in urls.items() if k in keys}
        return self._download_urls(urls, output_dir, decompress)

    def list_versions(self) -> list[str]:
        """List all Ensembl release numbers published on the FTP server.

        Returns:
            Sorted list of release-number strings, e.g. ``["100", ..., "115"]``.
        """
        with httpx.Client(follow_redirects=True, timeout=30.0) as client:
            response = client.get("https://ftp.ensembl.org/pub/")
            response.raise_for_status()
        versions = sorted({int(m) for m in re.findall(r"release-(\d+)/", response.text)})
        return [str(v) for v in versions]


__all__ = ["EnsemblDownloader", "EnsemblParser"]
