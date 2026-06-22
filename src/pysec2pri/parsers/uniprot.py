"""UniProt file parser for secondary-to-primary identifier mappings.

This parser extracts ID-to-ID mappings:
- Secondary accessions -> primary accessions (from sec_ac.txt)
- Deleted accessions -> sssom:NoTermFound (from delac_sp.txt)

Uses SSSOM-compliant IdMappingSet with cardinality computation.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
from sssom_schema import Mapping

from pysec2pri.parsers.base import (
    WITHDRAWN_ENTRY,
    WITHDRAWN_ENTRY_LABEL,
    BaseDownloader,
    BaseParser,
    Sec2PriMappingSet,
)


class UniProtParser(BaseParser):
    """Parser for UniProt files using Polars.

    Extracts secondary-to-primary UniProt accession mappings from
    sec_ac.txt (secondary accessions) and delac_sp.txt (deleted accessions).

    Returns IdMappingSet for all mappings (UniProt only has ID mappings).
    """

    datasource_name = "uniprot"

    @property
    def sec_ac_url(self) -> str:
        """Get the sec_ac.txt download URL from config."""
        return self.get_download_url("sec_ac") or ""

    @property
    def delac_url(self) -> str:
        """Get the delac_sp.txt download URL from config."""
        return self.get_download_url("delac_sp") or ""

    def parse(
        self,
        input_path: Path | str | None = None,
        delac_path: Path | str | None = None,
    ) -> Sec2PriMappingSet:
        """Parse UniProt mapping files into an IdMappingSet.

        Args:
            input_path: Path to sec_ac.txt (secondary accessions file).
            delac_path: Path to delac_sp.txt (deleted accessions file).

        Returns:
            IdMappingSet with computed cardinalities based on IDs.
        """
        # Resolve version
        self._resolve_version(
            Path(input_path)
            if input_path is not None
            else (Path(delac_path) if delac_path is not None else None)
        )

        mappings: list[Mapping] = []

        if input_path is not None:
            mappings.extend(self._parse_sec_ac(Path(input_path)))

        if delac_path is not None:
            mappings.extend(self._parse_delac(Path(delac_path)))

        return self._create_mapping_set(mappings)

    def _parse_sec_ac(self, file_path: Path) -> list[Mapping]:
        """Parse sec_ac.txt for secondary -> primary accession mappings.

        Args:
            file_path: Path to sec_ac.txt file.

        Returns:
            List of SSSOM Mapping objects.
        """
        # Count lines until (and including) the separator row.
        skip_rows = 0
        with file_path.open("r", encoding="utf-8") as f:
            for raw_line in f:
                skip_rows += 1
                if raw_line.startswith("_"):
                    break  # next line is first data row

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

        if df.is_empty():
            return []

        df = df.with_columns(
            pl.struct(["subject_id", "object_id"])
            .map_elements(
                lambda x: self._record_id(
                    meta_ns,
                    x["object_id"],
                    x["subject_id"],
                ),
                return_dtype=pl.Utf8,
                strategy="thread_local",
            )
            .alias("record_id"),
            pl.struct(["subject_id", "object_id"])
            .map_elements(
                lambda x: self._pair_hash(x["object_id"], x["subject_id"]),
                return_dtype=pl.Utf8,
                strategy="thread_local",
            )
            .alias("pair_key"),
        )

        from pysec2pri.consolidate import load_mapping_dates

        consolidated = load_mapping_dates("uniprot", mapping_sets="ids")

        m_meta = self.get_mapping_metadata()
        fixed = {
            "predicate_id": m_meta["predicate_id"],
            "predicate_label": m_meta.get("predicate_label"),
            "mapping_justification": m_meta["mapping_justification"],
            "subject_source": m_meta.get("subject_source"),
            "object_source": m_meta.get("object_source"),
            "mapping_tool": m_meta.get("mapping_tool"),
            "license": m_meta.get("license"),
        }
        rows = df.select(["subject_id", "object_id", "record_id", "pair_key"]).to_dicts()
        for row in rows:
            # The consolidated index is keyed by the version-independent pair
            # hash, not the whole record_id.
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
        with file_path.open("r", encoding="utf-8") as f:
            for raw_line in f:
                skip_rows += 1
                if raw_line.startswith("_"):
                    break  # next line is first deleted accession

        meta_ns = self._record_namespace()

        df = (
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
            .filter(
                pl.col("accession").str.contains(
                    r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
                )
            )
            .with_columns(
                pl.concat_str([pl.lit("UniProtKB:"), pl.col("accession")]).alias("object_id")
            )
            .select("object_id")
            .collect()
        )

        if df.is_empty():
            return []

        df = df.with_columns(
            pl.struct(["object_id"])
            .map_elements(
                lambda x: self._record_id(
                    meta_ns,
                    x["object_id"],
                    x["object_id"],  # subject_id == object_id in this dataset
                ),
                return_dtype=pl.Utf8,
                strategy="thread_local",
            )
            .alias("record_id"),
            pl.struct(["object_id"])
            .map_elements(
                lambda x: self._pair_hash(x["object_id"], x["object_id"]),
                return_dtype=pl.Utf8,
                strategy="thread_local",
            )
            .alias("pair_key"),
        )

        from pysec2pri.consolidate import load_mapping_dates

        consolidated = load_mapping_dates("uniprot", mapping_sets="ids")

        m_meta = self.get_mapping_metadata()
        fixed = {
            "subject_id": WITHDRAWN_ENTRY,
            "subject_label": WITHDRAWN_ENTRY_LABEL,
            "predicate_id": "oboInOwl:consider",
            "mapping_justification": m_meta["mapping_justification"],
            "subject_source": m_meta.get("subject_source"),
            "object_source": m_meta.get("object_source"),
            "mapping_tool": m_meta.get("mapping_tool"),
            "license": m_meta.get("license"),
            "comment": "Deleted accession with no replacement.",
        }
        rows = df.select(["object_id", "record_id", "pair_key"]).to_dicts()
        for row in rows:
            # The consolidated index is keyed by the version-independent pair
            # hash, not the whole record_id.
            row["mapping_date"] = consolidated.get(row.pop("pair_key"))
        return self._build_mappings(rows, fixed, desc="Processing delac", total=len(rows))

    def _create_mapping_set(
        self, mappings: list[Mapping], mapping_type: str = "id"
    ) -> Sec2PriMappingSet:
        """Delegate to base class method."""
        return self.create_mapping_set(mappings, mapping_type)

    def parse_primary_ids(
        self,
        acindex_path: Path | str | None = None,
    ) -> Sec2PriMappingSet:
        """Return a mapping set containing the full list of current UniProt primary ACs.

        Parses ``acindex.txt`` (or a gzip-compressed variant) to extract every
        accession number that currently appears in UniProtKB/Swiss-Prot.  The
        file lists one AC per row (after the ``__________`` separator line);
        only the first whitespace-delimited token of each data line is taken.

        For versioned (legacy) releases the file can be found at::

            https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/
            release-{version}/knowledgebase/docs/acindex.txt.gz

        Args:
            acindex_path: Local path to ``acindex.txt`` (plain or ``.gz``).
                Auto-downloaded from the current release when ``None``.

        Returns:
            :class:`~pysec2pri.parsers.base.IdMappingSet` with no mappings and
            ``_primary_ids`` populated with all current ``UniProtKB:<AC>`` CURIEs.
        """
        if acindex_path is None:
            from pysec2pri.api import _auto_download

            files, release_date = _auto_download("uniprot", None, keys=["acindex"])
            acindex_path = files["acindex"]
            self.release_date = release_date
        acindex_path = Path(str(acindex_path))
        self._resolve_version(acindex_path)

        primary_ids = self._extract_primary_ids_from_acindex(acindex_path)
        ms = self._create_mapping_set([], mapping_type="id")
        object.__setattr__(ms, "_primary_ids", primary_ids)
        return ms

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


class UniProtDownloader(BaseDownloader):
    """Downloader for UniProt data files."""

    datasource_name = "uniprot"

    def get_download_urls(
        self,
        version: str | None = None,
        **kwargs: object,
    ) -> dict[str, str]:
        """Get UniProt download URLs for *version*, or latest."""
        from pysec2pri.download import _get_uniprot_urls_for_version

        if version:
            return _get_uniprot_urls_for_version(version)
        if self._config:
            return dict(self._config.download_urls)
        raise ValueError("UniProt config not loaded")

    def download(
        self,
        output_dir: Path,
        version: str | None = None,
        decompress: bool = True,
        **kwargs: object,
    ) -> dict[str, Path]:
        """Download UniProt files into *output_dir*."""
        urls = self.get_download_urls(version)
        return self._download_urls(urls, output_dir, decompress)

    def list_versions(self) -> list[str]:
        """List all available UniProt previous-release versions.

        Scrapes the UniProt FTP previous_releases directory for version
        strings.

        Returns:
            Sorted list of version strings
            (e.g. ``["2024_01", "2024_02", ...]``).

        Raises:
            ValueError: If the archive URL is not configured.
        """
        import re

        import httpx

        if not self._config or not self._config.archive_url:
            raise ValueError("UniProt archive URL not configured")
        with httpx.Client(follow_redirects=True, timeout=30.0) as client:
            response = client.get(self._config.archive_url)
            response.raise_for_status()
        # FTP HTML index: links like "release-2024_01/"
        matches = re.findall(r'href="release-(\d{4}_\d{2})/', response.text)
        return sorted(set(matches))


__all__ = ["UniProtDownloader", "UniProtParser"]
