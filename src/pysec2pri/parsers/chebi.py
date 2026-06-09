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
          TSV paths.

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
        ms = self._create_mapping_set(mappings, mapping_type="id")
        # Populate primary IDs when compound data is available
        if cpd_path is not None and cpd_path.exists():
            object.__setattr__(ms, "_primary_ids", self._extract_primary_ids(cpd_path))
        return ms

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
        ms = self._create_mapping_set(mappings, mapping_type="label")
        # Populate primary symbols when compound data is available
        if cpd_path is not None and cpd_path.exists():
            object.__setattr__(ms, "_primary_symbols", self._extract_primary_symbols(cpd_path))
        return ms

    def parse_primary_ids(
        self,
        input_path: Path | str | None = None,
        *,
        compounds_path: Path | str | None = None,
    ) -> Sec2PriMappingSet:
        """Return a mapping set containing the full list of current ChEBI primary IDs.

        Reads ``compounds.tsv`` (TSV releases >= 245) to extract every current
        ChEBI compound ID.  The returned mapping set has an empty ``mappings``
        list; its ``_primary_ids`` store is populated with every current ChEBI
        ID (CHEBI: prefixed) so that ``to_pri_ids()`` produces the authoritative
        complete list.

        Args:
            input_path: Path to a directory containing ``compounds.tsv``, or
                directly to ``compounds.tsv`` itself.
            compounds_path: Explicit path to ``compounds.tsv`` (overrides
                ``input_path``).

        Returns:
            :class:`~pysec2pri.parsers.base.IdMappingSet` with no mappings and
            ``_primary_ids`` populated with all current ChEBI IDs.
        """
        cpd_path = self._resolve_compounds_path(input_path, compounds_path)
        if cpd_path is None or not cpd_path.exists():
            raise ValueError(
                "Could not locate compounds.tsv. Pass input_path (directory or file) "
                "or compounds_path."
            )
        self._resolve_version(cpd_path)
        ms = self._create_mapping_set([], mapping_type="id")
        object.__setattr__(ms, "_primary_ids", self._extract_primary_ids(cpd_path))
        return ms

    def parse_primary_symbols(
        self,
        input_path: Path | str | None = None,
        *,
        compounds_path: Path | str | None = None,
    ) -> Sec2PriMappingSet:
        """Return a mapping set containing the full list of current ChEBI compound names.

        Reads ``compounds.tsv`` to extract every current compound's canonical
        name.  The returned mapping set has an empty ``mappings`` list; its
        ``_primary_symbols`` store is populated.

        Args:
            input_path: Path to a directory containing ``compounds.tsv``, or
                directly to ``compounds.tsv`` itself.
            compounds_path: Explicit path to ``compounds.tsv``.

        Returns:
            :class:`~pysec2pri.parsers.base.LabelMappingSet` with no mappings
            and ``_primary_symbols`` populated.
        """
        cpd_path = self._resolve_compounds_path(input_path, compounds_path)
        if cpd_path is None or not cpd_path.exists():
            raise ValueError(
                "Could not locate compounds.tsv. Pass input_path (directory or file) "
                "or compounds_path."
            )
        self._resolve_version(cpd_path)
        ms = self._create_mapping_set([], mapping_type="label")
        object.__setattr__(ms, "_primary_symbols", self._extract_primary_symbols(cpd_path))
        return ms

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
                "subject_label": syn,  # synonym = secondary : subject
                "object_id": sid,
                "object_label": pname,  # primary name : object
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

    def _resolve_compounds_path(
        self,
        input_path: Path | str | None,
        compounds_path: Path | str | None,
    ) -> Path | None:
        """Resolve compounds.tsv path from either a directory or explicit path."""
        if compounds_path is not None:
            return Path(compounds_path)
        if input_path is not None:
            p = Path(input_path)
            if p.is_dir():
                candidate = p / "compounds.tsv"
                return candidate if candidate.exists() else None
            # Treat as direct path to compounds.tsv
            if p.name.startswith("compounds"):
                return p
        return None

    def _extract_primary_ids(self, compounds_path: Path) -> set[str]:
        """Extract all current ChEBI IDs from compounds.tsv.

        Reads the ``id`` column and returns a set of ``CHEBI:<n>`` CURIEs.
        When ``subset == "3star"`` only 3-star compounds are included.

        Args:
            compounds_path: Path to ``compounds.tsv``.

        Returns:
            Set of ``CHEBI:<id>`` strings.
        """
        cols = ["id", "stars"] if self.subset == "3star" else ["id"]
        df = pl.read_csv(
            compounds_path,
            separator="\t",
            columns=cols,
            schema_overrides={"id": pl.Int64},
        )
        if self.subset == "3star" and "stars" in df.columns:
            df = df.filter(pl.col("stars") == 3)
        return {f"CHEBI:{v}" for v in df["id"].drop_nulls().to_list()}

    def _extract_primary_symbols(self, compounds_path: Path) -> dict[str, set[str]]:
        """Extract all current ChEBI compound names from compounds.tsv.

        Returns a ``dict`` mapping each name to the set of primary
        ``CHEBI:<id>`` IDs that carry it.  When ``subset == "3star"`` only
        3-star compounds are included.

        Args:
            compounds_path: Path to ``compounds.tsv``.

        Returns:
            ``dict[name, set[CHEBI:<id>]]``
        """
        cols = ["id", "name", "stars"] if self.subset == "3star" else ["id", "name"]
        df = pl.read_csv(
            compounds_path,
            separator="\t",
            columns=cols,
            schema_overrides={"id": pl.Int64},
        )
        if self.subset == "3star" and "stars" in df.columns:
            df = df.filter(pl.col("stars") == 3)
        result: dict[str, set[str]] = {}
        for chebi_id, name in df.select(["id", "name"]).drop_nulls().rows():
            result.setdefault(str(name), set()).add(f"CHEBI:{chebi_id}")
        return result


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

    The authoritative primary (canonical) name for each compound is taken from
    the ``name`` column of ``compounds.tsv``: the same value displayed on the
    ChEBI website and stored as the top-level label in the ChEBI ontology.
    Every entry in ``names.tsv`` that differs from that canonical name is
    treated as a synonym.

    ``compounds_path`` is required because there is no source-authorised way
    to identify the primary name from ``names.tsv`` alone.

    Args:
        names_path: Path to names.tsv file.
        compounds_path: Path to compounds.tsv.  Required for primary-name
            resolution and, when subset=="3star", for star-count filtering.
        subset: "3star" or "complete".
        show_progress: Unused; kept for API compatibility.

    Returns:
        List of ``(subject_curie, primary_name, synonym)`` tuples where
        ``primary_name`` is the canonical name from ``compounds.tsv`` and
        ``synonym`` is the alternative name from ``names.tsv``.

    Raises:
        ValueError: When ``compounds_path`` is not provided.
    """
    if compounds_path is None:
        raise ValueError(
            "compounds_path is required by _parse_names_tsv to resolve primary names. "
            "Pass the path to compounds.tsv."
        )

    # Load canonical names from compounds.tsv (ChEBI-authoritative primary label)
    cpd_cols = ["id", "name", "stars"] if subset == "3star" else ["id", "name"]
    cpd_df = pl.read_csv(
        compounds_path,
        separator="\t",
        columns=cpd_cols,
        schema_overrides={"id": pl.Int64, "name": pl.Utf8},
        null_values=[""],
    )
    if subset == "3star" and "stars" in cpd_df.columns:
        cpd_df = cpd_df.filter(pl.col("stars") == 3)

    # Rename for the join: compounds.id : compound_id, compounds.name : primary_name
    cpd_df = cpd_df.select(
        [
            pl.col("id").alias("compound_id"),
            pl.col("name").alias("primary_name"),
        ]
    ).drop_nulls()

    # Load all alternative name entries from names.tsv
    df = pl.read_csv(
        names_path,
        separator="\t",
        columns=["compound_id", "name"],
        schema_overrides={"compound_id": pl.Int64, "name": pl.Utf8},
    )

    # Restrict to 3-star compounds (those present in cpd_df after star filtering)
    if subset == "3star":
        three_star_ids = set(cpd_df["compound_id"].to_list())
        df = df.filter(pl.col("compound_id").is_in(three_star_ids))

    # Join to attach the canonical primary name from compounds.tsv
    df_with_primary = df.join(cpd_df, on="compound_id")

    # Exclude trivial self-mappings (name in names.tsv == canonical name)
    synonyms_df = df_with_primary.filter(pl.col("name") != pl.col("primary_name"))

    # Build mapping tuples: (subject_curie, primary_name, synonym)
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
            # New TSV format (>= 245): flat keys in chebi.yaml
            return {
                "secondary_ids": download_urls["secondary_ids"].format(version=v),
                "names": download_urls["names"].format(version=v),
                "compounds": download_urls["compounds"].format(version=v),
            }
        else:
            # SDF format: keys live in chebi_sdf.yaml
            from pysec2pri.parsers.base import get_datasource_config

            sdf_config = get_datasource_config("chebi_sdf")
            sdf_urls = sdf_config.download_urls
            sdf_key = "sdf_3star" if subset == "3star" else "sdf_complete"
            url = sdf_urls[sdf_key].format(version=v)
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

    def list_versions(self) -> list[str]:
        """List all available ChEBI archive release numbers.

        Queries both the current archive (TSV releases >= 245) and the legacy
        archive (SDF-only releases < 245) and returns all release numbers
        sorted in ascending order.

        Returns:
            Sorted list of version strings (e.g. ``["100", "101", ..., "251"]``).

        Raises:
            ValueError: If the archive URL is not configured.
        """
        import re

        import httpx

        if not self._config or not self._config.archive_url:
            raise ValueError("ChEBI archive URL not configured")

        versions: set[str] = set()

        # Derive legacy archive base URL from the legacy SDF URL template in config:
        # e.g. ".../chebi_legacy/archive/rel{version}/SDF/..." -> ".../chebi_legacy/archive/"
        legacy_base: str | None = None
        legacy_urls: dict[str, str] = self._config.download_urls.get("legacy", {})
        legacy_template = next(iter(legacy_urls.values()), None)
        if legacy_template:
            # Strip from "rel{version}" onwards to get the directory listing URL
            idx = legacy_template.find("rel{version}")
            if idx != -1:
                legacy_base = legacy_template[:idx]

        with httpx.Client(follow_redirects=True, timeout=30.0) as client:
            # New archive (>= 245)
            response = client.get(self._config.archive_url)
            response.raise_for_status()
            versions.update(re.findall(r"rel(\d+)", response.text))

            # Legacy archive (< 245)
            if legacy_base:
                try:
                    legacy_response = client.get(legacy_base)
                    legacy_response.raise_for_status()
                    versions.update(re.findall(r"rel(\d+)", legacy_response.text))
                except httpx.HTTPError:
                    pass  # Legacy archive unavailable, still return what we have

        return sorted(versions, key=int)


__all__ = ["ChEBIDownloader", "ChEBIParser"]
