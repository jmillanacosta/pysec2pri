"""ChEBI parser for secondary-to-primary identifier mappings.

Supports two formats:
1. TSV flat files (releases >= 245): secondary_ids.tsv, names.tsv, compounds.tsv
2. SDF files (releases < 245): chebi_3_stars.sdf or ChEBI_complete.sdf

Extracts:
1. ID-to-ID mappings: secondary ChEBI IDs -> primary ChEBI IDs
2. Label-to-label mappings: synonyms -> primary names

Uses SSSOM-compliant MappingSet classes with cardinality computation.
"""

from __future__ import annotations

from functools import cache
from pathlib import Path
from typing import TYPE_CHECKING, Any

import polars as pl
from sssom_schema import Mapping

from pysec2pri.parsers.base import (
    BaseDownloader,
    BaseParser,
    Sec2PriMappingSet,
)

if TYPE_CHECKING:
    pass

# Threshold version where TSV format was introduced
NEW_FORMAT_VERSION = 245


class ChEBIParser(BaseParser):
    """Parser for ChEBI data files.

    Supports both TSV flat files (>= release 245) and (legacy) SDF files.
    Extracts secondary-to-primary ChEBI identifier mappings and
    name-to-synonym relationships.

    Returns an IdMappingSet for ID mappings (cardinality computed on IDs)
    and can optionally include synonym mappings via LabelMappingSet.
    """

    datasource_name = "chebi"

    def __init__(
        self,
        version: str | None = None,
        show_progress: bool = True,
        subset: str = "3star",
    ) -> None:
        """Initialize the ChEBI parser.

        Args:
            version: Version/release identifier for the datasource.
            show_progress: Whether to show progress bars during parsing.
            subset: "3star" or "complete" - which compound subset to use.
                For TSV format, filters by stars in compounds.tsv.
                For SDF format, determines which file to download.
        """
        super().__init__(version=version, show_progress=show_progress)
        self.subset = subset

    @property
    def source_url(self) -> str:
        """Get the default download URL from config."""
        return self.get_download_url("sdf") or ""

    def _is_new_format(self) -> bool:
        """Check if we should use new TSV format based on version."""
        if self.version is None:
            return True  # Default to new format for latest
        try:
            return int(self.version) >= NEW_FORMAT_VERSION
        except ValueError:
            return True  # Default to new if version is not numeric

    def parse(
        self,
        input_path: Path | str | None = None,
        *,
        secondary_ids_path: Path | str | None = None,
        compounds_path: Path | str | None = None,
    ) -> Sec2PriMappingSet:
        """Parse ChEBI data into an IdMappingSet.

        Accepts three calling conventions:

        - ``input_path`` is a **directory**: expects ``secondary_ids.tsv``
          (and optionally ``compounds.tsv``) inside it (TSV format >= 245).
        - ``input_path`` is an **SDF file**: legacy format (< 245).
        - Keyword args ``secondary_ids_path`` / ``compounds_path``: explicit
          TSV paths (kept for backwards compatibility).

        Args:
            input_path: Path to an SDF file, or a directory of TSV files.
            secondary_ids_path: Explicit path to secondary_ids.tsv (TSV format).
            compounds_path: Explicit path to compounds.tsv for 3-star filtering.

        Returns:
            IdMappingSet with computed cardinalities.
        """
        sid_path, cpd_path = self._resolve_tsv_paths(input_path, secondary_ids_path, compounds_path)
        if sid_path is not None:
            self._resolve_version(sid_path)
            raw_mappings = _parse_secondary_ids_tsv(
                sid_path,
                compounds_path=cpd_path,
                subset=self.subset,
                show_progress=self.show_progress,
            )
        elif input_path is not None:
            self._resolve_version(Path(input_path))
            raw_mappings, _ = _parse_chebi_sdf_fast(
                Path(input_path),
                show_progress=self.show_progress,
            )
        else:
            raise ValueError("Must provide input_path (SDF or TSV dir) or secondary_ids_path")

        mappings = self._build_id_mappings(raw_mappings)
        return self._create_mapping_set(mappings, mapping_type="id")

    def parse_synonyms(
        self,
        input_path: Path | str | None = None,
        *,
        names_path: Path | str | None = None,
        compounds_path: Path | str | None = None,
    ) -> Sec2PriMappingSet:
        """Parse ChEBI data into a LabelMappingSet for synonyms.

        Accepts three calling conventions:

        - ``input_path`` is a **directory**: expects ``names.tsv``
          (and optionally ``compounds.tsv``) inside it (TSV format >= 245).
        - ``input_path`` is an **SDF file**: legacy format (< 245).
        - Keyword args ``names_path`` / ``compounds_path``: explicit
          TSV paths (kept for backwards compatibility).

        Args:
            input_path: Path to an SDF file, or a directory of TSV files.
            names_path: Explicit path to names.tsv (TSV format).
            compounds_path: Explicit path to compounds.tsv for 3-star filtering.

        Returns:
            LabelMappingSet with computed cardinalities based on labels.
        """
        n_path, cpd_path = self._resolve_tsv_paths(
            input_path, names_path, compounds_path, tsv_key="names.tsv"
        )
        if n_path is not None:
            self._resolve_version(n_path)
            raw_mappings = _parse_names_tsv(
                n_path,
                compounds_path=cpd_path,
                subset=self.subset,
                show_progress=self.show_progress,
            )
        elif input_path is not None:
            self._resolve_version(Path(input_path))
            _, raw_mappings = _parse_chebi_sdf_fast(
                Path(input_path),
                show_progress=self.show_progress,
            )
        else:
            raise ValueError("Must provide input_path (SDF or TSV dir) or names_path")

        mappings = self._build_label_mappings(raw_mappings)
        return self._create_mapping_set(mappings, mapping_type="label")

    # Internal helpers

    def _resolve_tsv_paths(
        self,
        input_path: Path | str | None,
        primary_tsv: Path | str | None,
        compounds: Path | str | None,
        tsv_key: str = "secondary_ids.tsv",
    ) -> tuple[Path | None, Path | None]:
        """Resolve TSV file paths from either a directory or explicit paths.

        Returns (primary_tsv_path, compounds_path) or (None, None) if not TSV.
        """
        if primary_tsv is not None:
            cpd = Path(compounds) if compounds else None
            return Path(primary_tsv), cpd

        if input_path is not None:
            p = Path(input_path)
            if p.is_dir():
                tsv_file = p / tsv_key
                cpd_file = p / "compounds.tsv"
                return (
                    tsv_file if tsv_file.exists() else None,
                    cpd_file if cpd_file.exists() else None,
                )
        return None, None

    def _build_id_mappings(self, raw_id_mappings: list[tuple[str, str]]) -> list[Mapping]:
        """Build Mapping objects for secondary->primary ID mappings."""
        m_meta = self.get_mapping_metadata()
        fixed = {
            "predicate_id": m_meta["predicate_id"],
            "predicate_label": m_meta.get("predicate_label"),
            "mapping_justification": m_meta["mapping_justification"],
            "subject_source": m_meta.get("subject_source"),
            "object_source": m_meta.get("object_source"),
            "mapping_tool": m_meta.get("mapping_tool"),
            "confidence": m_meta.get("confidence"),
            "license": m_meta.get("license"),
        }
        rows = [{"subject_id": sec, "object_id": pri} for pri, sec in raw_id_mappings]
        return self._build_mappings(rows, fixed, desc="Creating ID mappings", total=len(rows))

    def _build_label_mappings(self, raw_name_mappings: list[tuple[str, str, str]]) -> list[Mapping]:
        """Build Mapping objects for label/synonym mappings."""
        m_meta = self.get_mapping_metadata()
        fixed = {
            "mapping_justification": m_meta["mapping_justification"],
            "subject_source": m_meta.get("subject_source"),
            "object_source": m_meta.get("object_source"),
            "mapping_tool": m_meta.get("mapping_tool"),
            "license": m_meta.get("license"),
        }
        rows = [
            {
                "subject_id": sid,
                "subject_label": pname,
                "object_id": sid,
                "object_label": syn,
                "_label_type": "alias",
            }
            for sid, pname, syn in raw_name_mappings
        ]
        return self._build_mappings(rows, fixed, desc="Creating synonym mappings", total=len(rows))

    def _create_mapping_set(
        self, mappings: list[Mapping], mapping_type: str = "id"
    ) -> Sec2PriMappingSet:
        """Create an IdMappingSet or LabelMappingSet with metadata from config.

        Delegates to BaseParser.create_mapping_set().
        """
        return self.create_mapping_set(mappings, mapping_type)


# TSV parsing functions (new format >= 245)


@cache
def _get_3star_compound_ids(
    compounds_path: Path,
    show_progress: bool = True,
) -> set[int]:
    """Get set of compound IDs with 3 stars from compounds.tsv.

    Args:
        compounds_path: Path to compounds.tsv file.
        show_progress: Whether to show progress.

    Returns:
        Set of compound IDs (as integers) with stars == 3.
    """
    df = pl.read_csv(
        compounds_path,
        separator="\t",
        columns=["id", "stars"],
        schema_overrides={"id": pl.Int64, "stars": pl.Int64},
    )

    three_star_ids = set(df.filter(pl.col("stars") == 3)["id"].to_list())

    return three_star_ids


def _parse_secondary_ids_tsv(
    secondary_ids_path: Path,
    compounds_path: Path | None = None,
    subset: str = "3star",
    show_progress: bool = True,
) -> list[tuple[str, str]]:
    """Parse secondary_ids.tsv into (primary_id, secondary_id) tuples.

    Format: compound_id | secondary_id

    Args:
        secondary_ids_path: Path to secondary_ids.tsv file.
        compounds_path: Path to compounds.tsv for 3-star filtering.
        subset: "3star" or "complete".
        show_progress: Whether to show progress.

    Returns:
        List of (primary_curie, secondary_curie) tuples.
    """
    df = pl.read_csv(
        secondary_ids_path,
        separator="\t",
        schema_overrides={"compound_id": pl.Int64, "secondary_id": pl.Int64},
    )

    # Filter to 3-star compounds if requested
    if subset == "3star" and compounds_path is not None:
        three_star_ids = _get_3star_compound_ids(compounds_path, show_progress)
        df = df.filter(pl.col("compound_id").is_in(three_star_ids))

    # Build mapping tuples with CHEBI: prefix
    mappings: list[tuple[str, str]] = df.select(
        [
            (pl.lit("CHEBI:") + pl.col("compound_id").cast(pl.Utf8)).alias("primary_id"),
            (pl.lit("CHEBI:") + pl.col("secondary_id").cast(pl.Utf8)).alias("secondary_id"),
        ]
    ).rows()

    return mappings


def _parse_names_tsv(
    names_path: Path,
    compounds_path: Path | None = None,
    subset: str = "3star",
    show_progress: bool = True,
) -> list[tuple[str, str, str]]:
    """Parse names.tsv into (subject_id, primary_name, synonym) tuples.

    Format: id|compound_id|name|type|status_id|adapted|language_code|ascii_name

    For each compound, the name with the smallest id is considered the primary
    name, and all other names are synonyms.

    Args:
        names_path: Path to names.tsv file.
        compounds_path: Path to compounds.tsv for 3-star filtering.
        subset: "3star" or "complete".
        show_progress: Whether to show progress.

    Returns:
        List of (subject_curie, primary_name, synonym) tuples.
    """
    df = pl.read_csv(
        names_path,
        separator="\t",
        columns=["id", "compound_id", "name"],
        schema_overrides={"id": pl.Int64, "compound_id": pl.Int64, "name": pl.Utf8},
    )

    # Filter to 3-star compounds if requested
    if subset == "3star" and compounds_path is not None:
        three_star_ids = _get_3star_compound_ids(compounds_path, show_progress)
        df = df.filter(pl.col("compound_id").is_in(three_star_ids))

    # For each compound, find the primary name (smallest id)
    # and all other names as synonyms
    primary_names = df.group_by("compound_id").agg(
        [
            pl.col("name").sort_by("id").first().alias("primary_name"),
        ]
    )

    # Join to get primary name for each row
    df_with_primary = df.join(primary_names, on="compound_id")

    # Filter to only synonym rows (where name != primary_name)
    synonyms_df = df_with_primary.filter(pl.col("name") != pl.col("primary_name"))

    # Build mapping tuples
    mappings: list[tuple[str, str, str]] = synonyms_df.select(
        [
            (pl.lit("CHEBI:") + pl.col("compound_id").cast(pl.Utf8)).alias("subject_id"),
            pl.col("primary_name"),
            pl.col("name"),
        ]
    ).rows()

    return mappings


# SDF parsing functions (legacy format < 245)


def _read_sdf_lines(input_path: Path) -> list[str]:
    """Read an SDF file into a list of lines."""
    with input_path.open("r", encoding="utf-8") as f:
        return f.readlines()


def _extract_value(lines: list[str], i: int, total_lines: int) -> tuple[str | None, int]:
    """Extract a value from the next line after a tag."""
    i += 1
    if i < total_lines:
        value = lines[i].strip()
        return value, i
    return None, i


def _extract_secondary_ids(
    lines: list[str],
    i: int,
    total_lines: int,
    current_subject_id: str | None,
    raw_id_mappings: list[tuple[str, str]],
) -> int:
    """Extract secondary IDs from the lines."""
    i += 1
    while i < total_lines:
        sec_line = lines[i].strip()
        if sec_line.startswith("CHEBI:"):
            if current_subject_id:
                # Handle semicolon-separated IDs on same line
                for sec_id in sec_line.split(";"):
                    sec_id = sec_id.strip()
                    if sec_id.startswith("CHEBI:"):
                        raw_id_mappings.append((current_subject_id, sec_id))
            i += 1
        else:
            break
    return i


def _extract_synonyms(
    lines: list[str],
    i: int,
    total_lines: int,
    current_subject_id: str | None,
    current_name: str | None,
    raw_name_mappings: list[tuple[str, str, str]],
) -> int:
    """Extract synonyms from the lines."""
    i += 1
    while i < total_lines:
        syn_line = lines[i].strip()
        if syn_line and not syn_line.startswith(">"):
            if current_subject_id and current_name:
                raw_name_mappings.append((current_subject_id, current_name, syn_line))
            i += 1
        else:
            break
    return i


def _process_sdf_line(
    lines: list[str],
    i: int,
    total_lines: int,
    current_subject_id: str | None,
    current_name: str | None,
    raw_id_mappings: list[tuple[str, str]],
    raw_name_mappings: list[tuple[str, str, str]],
) -> tuple[int, str | None, str | None, bool]:
    """Process a single SDF line and return updated state."""
    line = lines[i].strip()

    if "<chebi id>" in line.lower():
        value, i = _extract_value(lines, i, total_lines)
        if value:
            # Handle both "CHEBI:123" and "123" formats
            if value.startswith("CHEBI:"):
                value = f"CHEBI:{value[6:]}"  # Normalize
            else:
                value = f"CHEBI:{value}"
        return i, value, current_name, False

    if "<chebi name>" in line.lower():
        value, i = _extract_value(lines, i, total_lines)
        return i, current_subject_id, value, False

    if "<secondary" in line.lower():
        i = _extract_secondary_ids(lines, i, total_lines, current_subject_id, raw_id_mappings)
        return i, current_subject_id, current_name, True

    if "<synonym" in line.lower():
        i = _extract_synonyms(
            lines,
            i,
            total_lines,
            current_subject_id,
            current_name,
            raw_name_mappings,
        )
        return i, current_subject_id, current_name, True

    if line.startswith("$$$$"):
        return i, None, None, False

    return i, current_subject_id, current_name, False


def _parse_chebi_sdf_fast(
    input_path: Path,
    show_progress: bool = True,
) -> tuple[list[tuple[str, str]], list[tuple[str, str, str]]]:
    """Parse ChEBI SDF file into raw ID and synonym mappings.

    Args:
        input_path: Path to the SDF file.
        show_progress: Whether to display a progress bar.

    Returns:
        Tuple of:
            - list of (primary_id, secondary_id) for ID mappings
            - list of (subject_id, primary_name, synonym) for label mappings
    """
    from tqdm import tqdm

    lines = _read_sdf_lines(input_path)

    raw_id_mappings: list[tuple[str, str]] = []
    raw_name_mappings: list[tuple[str, str, str]] = []

    current_subject_id: str | None = None
    current_name: str | None = None

    i = 0
    total_lines = len(lines)

    pbar = tqdm(total=total_lines, desc="Parsing ChEBI SDF") if show_progress else None

    while i < total_lines:
        i, current_subject_id, current_name, skip = _process_sdf_line(
            lines,
            i,
            total_lines,
            current_subject_id,
            current_name,
            raw_id_mappings,
            raw_name_mappings,
        )
        if not skip:
            i += 1
        if pbar:
            pbar.update(1)

    if pbar:
        pbar.close()

    return raw_id_mappings, raw_name_mappings


class ChEBIDownloader(BaseDownloader):
    """Downloader for ChEBI data files.

    Supports both TSV flat files (>= release 245) and legacy SDF files.
    """

    datasource_name = "chebi"

    def __init__(
        self,
        version: str | None = None,
        show_progress: bool = True,
        subset: str = "3star",
        use_sdf: bool = False,
    ) -> None:
        """Initialize the ChEBI downloader.

        Args:
            version: Version/release identifier.
            show_progress: Whether to show progress bars.
            subset: "3star" or "complete" - which compound subset to use.
            use_sdf: Force SDF format even for releases >= 245.
        """
        super().__init__(version=version, show_progress=show_progress)
        self.subset = subset
        self.use_sdf = use_sdf

    def get_download_urls(
        self,
        version: str | None = None,
        **kwargs: Any,
    ) -> dict[str, str]:
        """Get ChEBI download URLs based on version.

        For version >= 245 (and use_sdf=False): returns TSV flat file URLs
        For version < 245 (or use_sdf=True): returns SDF file URL

        URLs are loaded from chebi.yaml config file.

        Args:
            version: Specific release version.
            **kwargs: Optional 'subset' and 'use_sdf' overrides.

        Returns:
            Dictionary with file URLs keyed by type.
        """
        v = version or self.version
        subset = kwargs.get("subset", self.subset)
        use_sdf = kwargs.get("use_sdf", self.use_sdf)

        if not self._config:
            raise ValueError("ChEBI config not loaded")

        download_urls = self._config.download_urls

        # Determine format: use new TSV if >= threshold AND not forcing SDF
        use_tsv = self.is_new_format(v) and not use_sdf

        if use_tsv:
            # New TSV format (>= 245)
            new_urls: dict[str, str] = download_urls.get("new", {})
            return {
                "secondary_ids": new_urls["secondary_ids"].format(version=v),
                "names": new_urls["names"].format(version=v),
                "compounds": new_urls["compounds"].format(version=v),
            }
        else:
            # SDF format (< 245 legacy OR >= 245 with --use-sdf)
            sdf_key = "sdf_3star" if subset == "3star" else "sdf_complete"

            if not self.is_new_format(v):
                # Legacy releases (< 245): use legacy URLs
                legacy_urls: dict[str, str] = download_urls.get("legacy", {})
                url = legacy_urls[sdf_key].format(version=v)
            else:
                # New releases (>= 245) with --use-sdf: use new SDF URLs
                new_urls = download_urls.get("new", {})
                url = new_urls[sdf_key].format(version=v)

            return {"sdf": url}

    def download(
        self,
        output_dir: Path,
        version: str | None = None,
        decompress: bool = True,
        **kwargs: Any,
    ) -> dict[str, Path]:
        """Download ChEBI files.

        Args:
            output_dir: Directory to save files.
            version: Specific version to download.
            decompress: Whether to decompress .gz files.
            **kwargs: Optional 'subset' and 'use_sdf' overrides.

        Returns:
            Dictionary mapping file keys to downloaded paths.
        """
        v = version or self.version
        urls = self.get_download_urls(v, **kwargs)

        return self._download_urls(urls, output_dir, decompress)

    def get_format(self, version: str | None = None) -> str:
        """Get the format that will be used for a given version.

        Args:
            version: Version to check. If None, uses self.version.

        Returns:
            "tsv" or "sdf"
        """
        v = version or self.version
        use_tsv = self.is_new_format(v) and not self.use_sdf
        return "tsv" if use_tsv else "sdf"


__all__ = ["ChEBIDownloader", "ChEBIParser"]
