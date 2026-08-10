"""UniProt file parser for secondary-to-primary identifier mappings.

This parser extracts ID-to-ID mappings:
- Secondary accessions -> primary accessions (from sec_ac.txt), optionally
  filtered by whether the primary is reviewed (Swiss-Prot) or not (TrEMBL)
- Deleted accessions -> sssom:NoTermFound (from delac_sp.txt; Swiss-Prot
  only.

Uses SSSOM-compliant IdMappingSet with cardinality computation.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
from sssom_schema import Mapping

from pysec2pri.logging import logger
from pysec2pri.parsers.base import (
    WITHDRAWN_ENTRY,
    WITHDRAWN_ENTRY_LABEL,
    BaseMappingSet,
    BaseParser,
)

#: UniProtKB accession syntax (https://www.uniprot.org/help/accession_numbers).
_ACCESSION_PATTERN = r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"


class UniProtParser(BaseParser):
    """Parser for UniProt files using Polars.

    Extracts secondary-to-primary UniProt accession mappings from
    sec_ac.txt (secondary accessions) and delac_sp.txt (deleted accessions).

    Returns IdMappingSet for all mappings (UniProt only has ID mappings).
    """

    datasource_name = "uniprot"

    def __init__(
        self,
        version: str | None = None,
        show_progress: bool = True,
        subset: str | None = None,
    ) -> None:
        """Initialize the UniProt parser.

        Args:
            version: Version/release identifier for the datasource.
            show_progress: Whether to show progress bars during parsing.
            subset: ``"swissprot"`` (reviewed primaries only), ``"trembl"``
                (unreviewed primaries only, no deletions), or ``"all"``
                (default; unfiltered, but still no TrEMBL deletions).
        """
        super().__init__(version=version, show_progress=show_progress)
        if subset is None and self._config is not None:
            subset = self._config.default_subset()
        self.subset = subset

    def parse(
        self,
        input_path: Path | str | None = None,
        delac_path: Path | str | None = None,
        acindex_path: Path | str | None = None,
    ) -> BaseMappingSet:
        """Parse UniProt mapping files into an IdMappingSet.

        Args:
            input_path: Path to sec_ac.txt (secondary accessions file).
            delac_path: Path to delac_sp.txt (deleted Swiss-Prot accessions).
            acindex_path: Path to acindex.txt (current Swiss-Prot
                accessions), used to split by reviewed/unreviewed primary.

        Returns:
            IdMappingSet with computed cardinalities based on IDs.
        """
        self._resolve_version(
            Path(input_path)
            if input_path is not None
            else (Path(delac_path) if delac_path is not None else None)
        )

        keep_reviewed = {"swissprot": True, "trembl": False}.get(self.subset or "")
        reviewed_acs = (
            self._extract_primary_ids_from_acindex(Path(acindex_path))
            if acindex_path is not None and keep_reviewed is not None
            else None
        )

        mappings: list[Mapping] = []

        if input_path is not None:
            mappings.extend(self._parse_sec_ac(Path(input_path), reviewed_acs, keep_reviewed))

        if self.subset != "trembl" and delac_path is not None:
            mappings.extend(self._parse_delac(Path(delac_path)))

        return self.create_mapping_set(mappings)

    def _parse_sec_ac(
        self,
        file_path: Path,
        reviewed_acs: set[str] | None = None,
        keep_reviewed: bool | None = None,
    ) -> list[Mapping]:
        """Parse sec_ac.txt for secondary -> primary accession mappings.

        Args:
            file_path: Path to sec_ac.txt file.
            reviewed_acs: Reviewed (Swiss-Prot) primary accessions, from
                ``acindex.txt``.
            keep_reviewed: Keep rows reviewed if ``True``, unreviewed if ``False``.

        Returns:
            List of SSSOM Mapping objects.
        """
        # Count lines until (and including) the separator row.
        skip_rows = 0
        found_separator = False
        with file_path.open("r", encoding="utf-8") as f:
            for raw_line in f:
                skip_rows += 1
                if raw_line.startswith("_"):
                    found_separator = True
                    break  # next line is first data row
        if not found_separator:
            logger.warning("No '_' separator line in %s; parsing from the first line.", file_path)
            skip_rows = 0

        meta_ns = self._record_namespace()

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
                pl.col("subject_id").is_not_null(),
                pl.col("object_id").is_not_null(),
                pl.col("subject_id").str.len_chars() > 0,
                pl.col("object_id").str.len_chars() > 0,
                pl.col("subject_id") != pl.col("object_id"),
            )
            .with_columns(
                pl.concat_str([pl.lit("UniProtKB:"), pl.col("subject_id")]).alias("subject_id"),
                pl.concat_str([pl.lit("UniProtKB:"), pl.col("object_id")]).alias("object_id"),
            )
            .collect()
        )

        if reviewed_acs is not None:
            is_reviewed = pl.col("object_id").is_in(reviewed_acs)
            df = df.filter(is_reviewed if keep_reviewed else ~is_reviewed)

        if df.is_empty():
            return []

        subj = df["subject_id"].to_list()
        obj = df["object_id"].to_list()
        df = df.with_columns(
            pl.Series(
                "record_id",
                [self._record_id(meta_ns, o, s) for o, s in zip(obj, subj, strict=True)],
                dtype=pl.Utf8,
            ),
            pl.Series(
                "pair_key",
                [self._pair_hash(o, s) for o, s in zip(obj, subj, strict=True)],
                dtype=pl.Utf8,
            ),
        )

        from pysec2pri.consolidate import load_mapping_dates

        consolidated = load_mapping_dates("uniprot", mapping_sets="ids")

        m_meta = self.get_mapping_metadata()
        fixed = {
            **self._fixed_mapping_fields(),
            "predicate_id": m_meta["predicate_id"],
            "predicate_label": m_meta.get("predicate_label"),
        }
        rows = df.select(["subject_id", "object_id", "record_id", "pair_key"]).to_dicts()
        for row in rows:
            # The consolidated index is keyed by the pair hash
            row["mapping_date"] = consolidated.get(row.pop("pair_key"))
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
        found_separator = False
        with file_path.open("r", encoding="utf-8") as f:
            for raw_line in f:
                skip_rows += 1
                if raw_line.startswith("_"):
                    found_separator = True
                    break  # next line is first deleted accession
        if not found_separator:
            logger.warning("No '_' separator line in %s; parsing from the first line.", file_path)
            skip_rows = 0

        meta_ns = self._record_namespace()

        accessions = (
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
            .filter(pl.col("accession").str.len_chars() > 0)
            .collect()
        )

        df = accessions.filter(pl.col("accession").str.contains(_ACCESSION_PATTERN))
        dropped = accessions.height - df.height
        if dropped:
            logger.warning(
                "%s: dropped %d of %d entries not matching the UniProt accession pattern.",
                file_path.name,
                dropped,
                accessions.height,
            )

        if df.is_empty():
            return []

        df = df.select(
            pl.concat_str([pl.lit("UniProtKB:"), pl.col("accession")]).alias("object_id")
        )

        # subject_id == object_id in this dataset
        obj = df["object_id"].to_list()
        df = df.with_columns(
            pl.Series("record_id", [self._record_id(meta_ns, o, o) for o in obj], dtype=pl.Utf8),
            pl.Series("pair_key", [self._pair_hash(o, o) for o in obj], dtype=pl.Utf8),
        )

        from pysec2pri.consolidate import load_mapping_dates

        consolidated = load_mapping_dates("uniprot", mapping_sets="ids")

        fixed = {
            **self._fixed_mapping_fields(),
            "subject_id": WITHDRAWN_ENTRY,
            "subject_label": WITHDRAWN_ENTRY_LABEL,
            "predicate_id": "oboInOwl:consider",
            "comment": "Deleted accession with no replacement.",
        }
        rows = df.select(["object_id", "record_id", "pair_key"]).to_dicts()
        for row in rows:
            # The consolidated index is keyed by the pair hash.
            row["mapping_date"] = consolidated.get(row.pop("pair_key"))
        return self._build_mappings(rows, fixed, desc="Processing delac", total=len(rows))

    def parse_primary_ids(
        self,
        acindex_path: Path | str | None = None,
    ) -> BaseMappingSet:
        """Return a mapping set containing the full list of current UniProt primary ACs.

        Parses ``acindex.txt`` (or a gzip-compressed variant) to extract every
        accession number that currently appears in UniProtKB/Swiss-Prot.

        For versioned (archived) releases, the file is extracted from
        ``knowledgebase-docs-only{version}.tar.gz`` (see :mod:`pysec2pri.downloads.uniprot`).

        Args:
            acindex_path: Local path to ``acindex.txt`` (plain or ``.gz``).
                Auto-downloaded from the current release when ``None``.

        Returns:
            :class:`~pysec2pri.parsers.base.IdMappingSet` with no mappings and
            ``_primary_ids`` populated with all current ``UniProtKB:<AC>`` CURIEs.
        """
        if acindex_path is None:
            from pysec2pri.api import _auto_download

            files, version, release_date = _auto_download("uniprot", None, keys=["acindex"])
            acindex_path = files["acindex"]
            self.version = version
            self.release_date = release_date
        acindex_path = Path(str(acindex_path))
        self._resolve_version(acindex_path)

        primary_ids = self._extract_primary_ids_from_acindex(acindex_path)
        return self.create_mapping_set([], mapping_type="id", primary_ids=primary_ids)

    def _extract_primary_ids_from_acindex(self, file_path: Path) -> set[str]:
        """Parse ``acindex.txt`` and return the set of all AC numbers.

        Skips the header block (everything up to and including the ``__________``
        separator line) and extracts the first whitespace-delimited token of
        each subsequent non-empty line.

        Args:
            file_path: Path to ``acindex.txt`` (plain or ``.gz``).

        Returns:
            Set of ``UniProtKB:<AC>`` CURIEs.
        """
        import gzip

        opener = gzip.open if file_path.suffix == ".gz" else open
        primary_ids: set[str] = set()
        in_data = False
        with opener(file_path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                stripped = line.strip()
                if not in_data:
                    if stripped.startswith("__"):
                        in_data = True
                    continue
                if not stripped:
                    continue
                token = stripped.split()[0]
                if token:
                    primary_ids.add(f"UniProtKB:{token}")
        return primary_ids


__all__ = ["UniProtParser"]
