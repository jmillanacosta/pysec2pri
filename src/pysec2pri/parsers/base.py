"""Base parser class for all datasource parsers."""

from __future__ import annotations

import re
from abc import ABC, abstractmethod
from collections.abc import Iterable
from dataclasses import dataclass, field
from datetime import date
from importlib import resources as _importlib_resources
from pathlib import Path
from typing import TYPE_CHECKING, Any, TypeVar, cast

import yaml
from sssom_schema import Mapping, MappingCardinalityEnum, MappingSet
from tqdm import tqdm

from pysec2pri.logging import logger
from pysec2pri.version import VERSION

if TYPE_CHECKING:
    import pandas as pd
    import rdflib
    from sssom import sssom_document

_T = TypeVar("_T")

# Values for withdrawn entries

WITHDRAWN_ENTRY = "sssom:NoTermFound"
WITHDRAWN_ENTRY_LABEL = "Withdrawn Entry"

# Config directory path

CONFIG_DIR = Path(_importlib_resources.files("pysec2pri.config"))  # type: ignore[arg-type]


@dataclass
class DatasourceConfig:
    """Configuration for a biological database datasource loaded from YAML."""

    name: str
    prefix: str
    curie_base_url: str
    default_output_filename: str = ""
    available_outputs: list[str] = field(default_factory=list)
    download_urls: dict[str, Any] = field(default_factory=dict)
    primary_file_key: str = ""
    id_pattern: str = ""
    archive_url: str = ""
    input_file_types: list[str] = field(default_factory=list)
    source: str = ""
    homepage: str = ""
    data_license: str = ""
    # SPARQL-based datasources (e.g., Wikidata)
    sparql_endpoint: str = ""
    queries: dict[str, str] = field(default_factory=dict)
    # For now, only ChEBI: version threshold for new TSV format.
    # Use if release files change location or serialization.
    new_format_version: int | None = None
    # Full metadata from YAML
    mappingset_metadata: dict[str, Any] = field(default_factory=dict)
    mapping_metadata: dict[str, Any] = field(default_factory=dict)


def load_config(datasource_name: str) -> dict[str, Any]:
    """Load configuration from a YAML file for a datasource.

    Args:
        datasource_name: Name of the datasource (e.g., 'chebi', 'hgnc').

    Returns:
        Dictionary with the full YAML configuration.

    Raises:
        FileNotFoundError: If the config file does not exist.
    """
    config_path = CONFIG_DIR / f"{datasource_name.lower()}.yaml"
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with config_path.open("r", encoding="utf-8") as f:
        result: dict[str, Any] = yaml.safe_load(f)
        return result


def get_datasource_config(datasource_name: str) -> DatasourceConfig:
    """Load and parse a DatasourceConfig from YAML.

    Args:
        datasource_name: Name of the datasource (e.g., 'chebi', 'hgnc').

    Returns:
        DatasourceConfig object populated from YAML.
    """
    raw = load_config(datasource_name)

    # Extract curie_map to get prefix and curie_base_url
    curie_map = raw.get("mappingset", {}).get("curie_map", {})
    # The first entry that's not a standard prefix is the datasource prefix
    standard_prefixes = {"IAO", "oboInOwl", "semapv", "skos", "sssom", "sec2pri"}
    prefix = ""
    curie_base_url = ""
    for k, v in curie_map.items():
        if k not in standard_prefixes:
            prefix = k
            curie_base_url = v
            break

    # Determine name - special case for wikidata
    if datasource_name.lower() == "wikidata":
        name = "Wikidata"
    else:
        name = datasource_name.upper()

    return DatasourceConfig(
        name=name,
        prefix=prefix,
        curie_base_url=curie_base_url,
        default_output_filename=raw.get("default_output_filename", ""),
        available_outputs=raw.get("available_outputs", []),
        download_urls=raw.get("download_urls", {}),
        primary_file_key=raw.get("primary_file_key", ""),
        id_pattern=raw.get("id_pattern", ""),
        archive_url=raw.get("archive_url", ""),
        input_file_types=raw.get("input_file_types", []),
        source=raw.get("source", ""),
        homepage=raw.get("homepage", ""),
        data_license=raw.get("data_license", ""),
        sparql_endpoint=raw.get("sparql_endpoint", ""),
        queries=raw.get("queries", {}),
        new_format_version=raw.get("new_format_version"),
        mappingset_metadata=raw.get("mappingset", {}),
        mapping_metadata=raw.get("mapping", {}),
    )


# Base Downloader Class


class BaseDownloader(ABC):
    """Abstract base class for datasource downloaders.

    Provides shared download logic that can be inherited by datasource-specific
    downloaders. Handles file downloads, URL construction, and version detection.
    """

    datasource_name: str = ""
    _config: DatasourceConfig | None = None

    def __init__(
        self,
        version: str | None = None,
        show_progress: bool = True,
    ) -> None:
        """Initialize the downloader.

        Args:
            version: Version/release identifier for the datasource.
            show_progress: Whether to show progress bars during downloads.
        """
        self.version = version
        self.show_progress = show_progress

        # Load config from YAML
        if self.datasource_name:
            try:
                self._config = get_datasource_config(self.datasource_name.lower())
            except FileNotFoundError:
                self._config = None

    @property
    def config(self) -> DatasourceConfig | None:
        """Get the loaded configuration."""
        return self._config

    @property
    def new_format_version(self) -> int | None:
        """Get the version threshold for new format (if any)."""
        if self._config:
            return self._config.new_format_version
        return None

    def is_new_format(self, version: str | None = None) -> bool:
        """Check if a version uses the new format.

        Args:
            version: Version to check. If None, uses self.version.

        Returns:
            True if version >= new_format_version threshold.
        """
        v = version or self.version
        threshold = self.new_format_version

        if threshold is None:
            return True  # No threshold means always "new" format

        if v is None:
            return True  # Default to new format for latest

        try:
            return int(v) >= threshold
        except ValueError:
            return True  # Default to new if version is not numeric

    @abstractmethod
    def get_download_urls(
        self,
        version: str | None = None,
        **kwargs: Any,
    ) -> dict[str, str]:
        """Get download URLs for the datasource.

        Args:
            version: Specific version to get URLs for.
            **kwargs: Additional options (e.g., subset, force_format).

        Returns:
            Dictionary mapping file keys to URLs.
        """

    @abstractmethod
    def download(
        self,
        output_dir: Path,
        version: str | None = None,
        decompress: bool = True,
        **kwargs: Any,
    ) -> dict[str, Path]:
        """Download files for the datasource.

        Args:
            output_dir: Directory to save downloaded files.
            version: Specific version to download.
            decompress: Whether to decompress .gz files.
            **kwargs: Additional options.

        Returns:
            Dictionary mapping file keys to downloaded paths.
        """

    def _download_file(
        self,
        url: str,
        output_path: Path,
        decompress_gz: bool = True,
        timeout: float | None = None,
        description: str | None = None,
    ) -> Path:
        """Download a file from URL to the specified path.

        Args:
            url: URL to download from.
            output_path: Where to save the file.
            decompress_gz: Whether to decompress .gz files automatically.
            timeout: Request timeout in seconds.
            description: Description for the progress bar.

        Returns:
            Path to the downloaded file.
        """
        # Import here to avoid circular imports
        from pysec2pri.download import download_file

        return download_file(
            url,
            output_path,
            decompress_gz=decompress_gz,
            timeout=timeout,
            show_progress=self.show_progress,
            description=description,
        )

    def _download_urls(
        self,
        urls: dict[str, str],
        output_dir: Path,
        decompress: bool = True,
    ) -> dict[str, Path]:
        """Download files from URLs to output directory.

        Args:
            urls: Dictionary mapping file keys to URLs.
            output_dir: Directory to save files.
            decompress: Whether to decompress .gz files.

        Returns:
            Dictionary mapping file keys to downloaded paths.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        downloaded: dict[str, Path] = {}

        for key, url in urls.items():
            filename = url.split("/")[-1]

            if decompress and filename.endswith(".gz"):
                filename = filename[:-3]

            output_path = output_dir / filename
            logger.info("Downloading %s: %s", key, url)
            self._download_file(url, output_path, decompress_gz=decompress)
            downloaded[key] = output_path
            logger.info("Saved to: %s", output_path)

        return downloaded


class Sec2PriMappingSet(MappingSet):  # type: ignore[misc]
    """A MappingSet for Sec2Pri, with helpers for cardinality and export.

    Attributes:
        _primary_ids: Private store for the full authoritative primary ID set.
            Kept private so sssom serialisers never include it in any output.
            Access it only through :meth:`to_pri_ids`.
    """

    # Primaries are private to sssom's schema
    # Populated by parsers that have access to the full primary ID list
    # (e.g. HGNCParser when the complete set file is provided).
    _primary_ids: set[str]

    def __init__(self, *args: object, **kwargs: object) -> None:
        """Initialise the mapping set and the private primary-IDs store."""
        super().__init__(*args, **kwargs)
        object.__setattr__(self, "_primary_ids", set())

    # Export helpers

    def _default_stem(self) -> str:
        """Derive a base filename stem from mapping set metadata."""
        ms_id: str = str(getattr(self, "mapping_set_id", None) or "")
        if ms_id:
            stem = ms_id.rstrip("/").rsplit("/", 1)[-1]
        else:
            stem = str(getattr(self, "mapping_set_title", None) or "mapping_set")
            stem = stem.lower().replace(" ", "_")
        version = getattr(self, "mapping_set_version", None)
        if version:
            stem = f"{stem}_{version}"
        return stem

    def _resolve_path(self, output_path: Path | str | None, suffix: str) -> Path:
        """Return *output_path* if given, else auto-generate one."""
        if output_path is not None:
            return Path(output_path)
        return Path(f"{self._default_stem()}{suffix}")

    def to_sssom(self, output_path: Path | str | None = None) -> sssom_document.MappingSetDocument:
        """Return an SSSOM ``MappingSetDocument``, optionally writing to TSV.

        Args:
            output_path: If given, the document is also serialised to an SSSOM
                TSV file at this path

        Returns:
            :class:`sssom.sssom_document.MappingSetDocument` for the mapping set.
        """
        import curies
        from sssom.sssom_document import MappingSetDocument

        raw_curie_map: object = self.curie_map or {}
        records: list[curies.Record] = []
        if isinstance(raw_curie_map, dict):
            for k, v in raw_curie_map.items():
                if isinstance(v, str):
                    uri_prefix: str = v
                elif hasattr(v, "prefix_url"):
                    uri_prefix = cast(str, v.prefix_url)
                else:
                    continue
                records.append(curies.Record(prefix=k, uri_prefix=uri_prefix))
        converter = curies.Converter(records=records)
        doc = MappingSetDocument(mapping_set=self, converter=converter)

        if output_path is not None:
            from pysec2pri.exports import write_sssom

            write_sssom(self, self._resolve_path(output_path, "_sssom.tsv"))

        return doc

    def to_rdf(
        self,
        output_path: Path | str | None = None,
        serialisation: str = "turtle",
    ) -> rdflib.Graph:
        """Return an RDFLib graph, optionally writing it to a file.

        When *output_path* is given (or auto-generated via the ``save``
        dispatcher), the graph is also serialised to disk.  Either way
        the :class:`rdflib.Graph` is returned so callers can query or
        manipulate it directly.

        Args:
            output_path: Destination path. Pass a path (or ``None`` to
                auto-generate one) to persist the graph.  If you only want
                the in-memory graph without touching the file-system, call
                ``to_rdf()`` with no arguments and ignore the path attribute.
            serialisation: RDFLib serialisation format (default: ``"turtle"``).

        Returns:
            :class:`rdflib.Graph` containing all mappings as RDF triples.
        """
        import io

        import rdflib
        from sssom.writers import write_rdf as _sssom_write_rdf

        from pysec2pri.exports import _to_msdf_via_sssom_parser, write_rdf

        msdf = _to_msdf_via_sssom_parser(self)
        if msdf is None:
            raise ValueError("Failed to convert mapping set to RDF.")

        buf = io.StringIO()
        _sssom_write_rdf(msdf, buf, serialisation=serialisation)
        g = rdflib.Graph()
        g.parse(data=buf.getvalue(), format=serialisation)

        if output_path is not None:
            write_rdf(self, self._resolve_path(output_path, ".ttl"), serialisation=serialisation)

        return g

    def to_json(self, output_path: Path | str | None = None) -> dict[str, Any]:
        """Return the mapping set as a JSON-compatible ``dict``, optionally writing to file.

        Args:
            output_path: If given, the JSON is also written to this path.

        Returns:
            ``dict`` representation of the mapping set in SSSOM JSON format.
        """
        import io
        import json

        from sssom.writers import write_json as _sssom_write_json

        from pysec2pri.exports import _to_msdf_via_sssom_parser

        msdf = _to_msdf_via_sssom_parser(self)
        if msdf is None:
            raise ValueError("Failed to convert mapping set to JSON.")
        buf = io.StringIO()
        _sssom_write_json(msdf, buf)
        data: dict[str, Any] = json.loads(buf.getvalue())

        if output_path is not None:
            path = self._resolve_path(output_path, ".json")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(buf.getvalue(), encoding="utf-8")

        return data

    def to_owl(
        self, output_path: Path | str | None = None, serialisation: str = "turtle"
    ) -> rdflib.Graph:
        """Return an OWL ``rdflib.Graph``, optionally writing to file.

        Args:
            output_path: If given, the graph is also serialised to this path.
            serialisation: RDFLib serialisation format (default: ``"turtle"``).

        Returns:
            :class:`rdflib.Graph` containing OWL axioms for the mapping set.
        """
        import io

        import rdflib
        from sssom.writers import write_owl as _sssom_write_owl

        from pysec2pri.exports import _to_msdf_via_sssom_parser

        msdf = _to_msdf_via_sssom_parser(self)
        if msdf is None:
            raise ValueError("Failed to convert mapping set to OWL.")
        buf = io.StringIO()
        _sssom_write_owl(msdf, buf, serialisation=serialisation)
        g = rdflib.Graph()
        g.parse(data=buf.getvalue(), format=serialisation)

        if output_path is not None:
            from pysec2pri.exports import write_owl

            write_owl(
                self, self._resolve_path(output_path, "_owl.ttl"), serialisation=serialisation
            )

        return g

    def _save_shared(
        self,
        fmt: str,
        output_path: Path | str | None,
        **kwargs: object,
    ) -> Path | None:
        """Write one of the shared formats (sssom/rdf/json/owl).

        Returns the written :class:`Path`, or ``None`` if *fmt* is not a
        shared format (caller should handle it).
        """
        if fmt in ("rdf", "owl", "json"):
            from collections.abc import Callable as _Callable

            from pysec2pri.exports import write_json, write_owl, write_rdf

            _write_fns: dict[str, tuple[_Callable[..., Path], str]] = {
                "rdf": (write_rdf, ".ttl"),
                "owl": (write_owl, "_owl.ttl"),
                "json": (write_json, ".json"),
            }
            fn, suffix = _write_fns[fmt]
            return fn(self, self._resolve_path(output_path, suffix), **kwargs)

        if fmt == "sssom":
            from pysec2pri.exports import write_sssom

            return write_sssom(self, self._resolve_path(output_path, "_sssom.tsv"))

        return None

    def save(
        self,
        fmt: str,
        output_path: Path | str | None = None,
        **kwargs: object,
    ) -> Path:
        """Write to any supported format by name.

        Shared formats: ``"sssom"``, ``"rdf"``, ``"json"``, ``"owl"``.
        Subclasses override this to add type-specific formats.

        Args:
            fmt: Format key (see above).
            output_path: Destination path. Auto-generated if ``None``.
            **kwargs: Forwarded to the format-specific writer.

        Returns:
            Path to the written file.

        Raises:
            ValueError: For unknown format keys.
        """
        shared = self._save_shared(fmt, output_path, **kwargs)
        if shared is not None:
            return shared
        raise ValueError(f"Unknown format {fmt!r}. Choose from: json, owl, rdf, sssom")

    # Cardinality helpers

    def _compute_cardinalities(self, on: str = "id") -> None:
        """Compute and set mapping_cardinality on all mappings.

        'on' can be 'id' (uses subject_id/object_id) or 'label'.
        """
        if not self.mappings:  # type: ignore[has-type]
            return

        mappings = self._normalize_mappings()

        if on == "label":
            sec_field, pri_field = "subject_label", "object_label"
        else:
            sec_field, pri_field = "subject_id", "object_id"

        import polars as pl

        sec_vals = [str(getattr(m, sec_field, None) or "") for m in mappings]
        pri_vals = [str(getattr(m, pri_field, None) or "") for m in mappings]

        cardinalities: list[str] = (
            pl.DataFrame({"sec": sec_vals, "pri": pri_vals})
            .with_columns(
                pl.col("sec").count().over("sec").alias("s_count"),
                pl.col("pri").count().over("pri").alias("p_count"),
            )
            .select(
                pl.when(pl.col("s_count") == 1, pl.col("p_count") == 1)
                .then(pl.lit("1:1"))
                .when(pl.col("s_count") > 1, pl.col("p_count") == 1)
                .then(pl.lit("n:1"))
                .when(pl.col("s_count") == 1, pl.col("p_count") > 1)
                .then(pl.lit("1:n"))
                .when(pl.col("s_count") > 1, pl.col("p_count") > 1)
                .then(pl.lit("n:n"))
                .otherwise(pl.lit("1:0"))
                .alias("cardinality")
            )
            .get_column("cardinality")
            .to_list()
        )

        for m, card in zip(mappings, cardinalities, strict=False):
            m.mapping_cardinality = MappingCardinalityEnum(card)

        self.mappings = mappings

    def _normalize_mappings(self) -> list[Mapping]:
        """Normalize mappings to a list of Mapping objects.

        Returns:
            List of Mapping objects.
        """
        mappings = self.mappings
        if not isinstance(mappings, list):
            mappings = [mappings]
        for i, m in enumerate(mappings):
            if isinstance(m, dict):
                mappings[i] = Mapping(**m)
        return mappings


class IdMappingSet(Sec2PriMappingSet):
    """Mapping set for ID-based (secondary to primary identifier) mappings."""

    def compute_cardinalities(self) -> None:
        """Compute cardinalities using subject_id and object_id fields."""
        self._compute_cardinalities(on="id")

    def to_sec2pri(self, output_path: Path | str | None = None) -> pd.DataFrame:
        """Return a ``DataFrame`` of secondary to primary ID mappings.

        Columns: ``subject_id`` (secondary), ``object_id`` (primary),
        ``predicate_id``, ``mapping_cardinality``.

        Args:
            output_path: If given, the DataFrame is also written as a TSV file.

        Returns:
            :class:`pandas.DataFrame` with one row per mapping.
        """
        import pandas as pd

        rows = [
            {
                "subject_id": str(getattr(m, "subject_id", "") or ""),
                "object_id": str(getattr(m, "object_id", "") or ""),
                "predicate_id": str(getattr(m, "predicate_id", "") or ""),
                "mapping_cardinality": str(getattr(m, "mapping_cardinality", "") or ""),
            }
            for m in (self.mappings or [])
        ]
        df = pd.DataFrame(
            rows, columns=["subject_id", "object_id", "predicate_id", "mapping_cardinality"]
        )

        if output_path is not None:
            path = self._resolve_path(output_path, "_sec2pri.tsv")
            path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(path, sep="\t", index=False)

        return df

    def to_pri_ids(self, output_path: Path | str | None = None) -> list[str]:
        """Return a sorted list of unique primary IDs, optionally writing to TXT.

        When ``_primary_ids`` is populated (e.g. from the HGNC complete set)
        that set is used.  Otherwise primary IDs are derived from the unique
        ``object_id`` values in the mappings.

        Args:
            output_path: If given, the IDs are also written one-per-line to a
                text file.

        Returns:
            Sorted list of unique primary ID strings.
        """
        private: set[str] = (
            object.__getattribute__(self, "_primary_ids")
            if hasattr(self, "_primary_ids")
            else set()
        )

        if private:
            ids = sorted(private)
        else:
            ids = sorted(
                {str(getattr(m, "object_id", None) or "") for m in (self.mappings or [])} - {""}
            )

        if output_path is not None:
            path = self._resolve_path(output_path, "_pri_ids.txt")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("\n".join(ids) + "\n", encoding="utf-8")

        return ids

    def save(
        self,
        fmt: str,
        output_path: Path | str | None = None,
        **kwargs: object,
    ) -> Path:
        """Write to any supported format by name.

        Formats: ``"sssom"``, ``"rdf"``, ``"json"``, ``"owl"``,
        ``"sec2pri"``, ``"pri_ids"``.

        Args:
            fmt: Format key (see above).
            output_path: Destination path. Auto-generated if ``None``.
            **kwargs: Forwarded to the format-specific writer.

        Returns:
            Path to the written file.

        Raises:
            ValueError: For unknown format keys.
        """
        shared = self._save_shared(fmt, output_path, **kwargs)
        if shared is not None:
            return shared

        if fmt == "sec2pri":
            self.to_sec2pri(output_path)
            return self._resolve_path(output_path, "_sec2pri.tsv")

        if fmt == "pri_ids":
            self.to_pri_ids(output_path)
            return self._resolve_path(output_path, "_pri_ids.txt")

        raise ValueError(
            f"Unknown format {fmt!r}. Choose from: json, owl, pri_ids, rdf, sec2pri, sssom"
        )


class LabelMappingSet(Sec2PriMappingSet):
    """Mapping set for label-based (previous/alias symbol to current symbol) mappings."""

    def compute_cardinalities(self) -> None:
        """Compute cardinalities using subject_label and object_label."""
        self._compute_cardinalities(on="label")

    def to_symbol_sec2pri(self, output_path: Path | str | None = None) -> pd.DataFrame:
        """Return a ``DataFrame`` of previous/alias symbol to current symbol mappings.

        Columns: ``subject_id``, ``subject_label`` (secondary/previous symbol),
        ``object_id``, ``object_label`` (primary/current symbol),
        ``predicate_id``, ``mapping_cardinality``.

        Args:
            output_path: If given, the DataFrame is also written as a TSV file.

        Returns:
            :class:`pandas.DataFrame` with one row per symbol mapping.
        """
        import pandas as pd

        rows = [
            {
                "subject_id": str(getattr(m, "subject_id", "") or ""),
                "subject_label": str(getattr(m, "subject_label", "") or ""),
                "object_id": str(getattr(m, "object_id", "") or ""),
                "object_label": str(getattr(m, "object_label", "") or ""),
                "predicate_id": str(getattr(m, "predicate_id", "") or ""),
                "mapping_cardinality": str(getattr(m, "mapping_cardinality", "") or ""),
            }
            for m in (self.mappings or [])
        ]
        df = pd.DataFrame(
            rows,
            columns=[
                "subject_id",
                "subject_label",
                "object_id",
                "object_label",
                "predicate_id",
                "mapping_cardinality",
            ],
        )

        if output_path is not None:
            path = self._resolve_path(output_path, "_symbol_sec2pri.tsv")
            path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(path, sep="\t", index=False)

        return df

    def to_pri_symbols(self, output_path: Path | str | None = None) -> list[str]:
        """Return a sorted list of unique current/primary symbols, optionally writing to TXT.

        Derived from the unique ``object_label`` values in the mappings.

        Args:
            output_path: If given, the symbols are also written one-per-line to
                a text file.

        Returns:
            Sorted list of unique primary symbol strings.
        """
        symbols = sorted(
            {str(getattr(m, "object_label", None) or "") for m in (self.mappings or [])} - {""}
        )

        if output_path is not None:
            path = self._resolve_path(output_path, "_pri_symbols.txt")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("\n".join(symbols) + "\n", encoding="utf-8")

        return symbols

    def to_name2synonym(self, output_path: Path | str | None = None) -> pd.DataFrame:
        """Return a name to synonym ``DataFrame``, optionally writing to TSV.

        Columns: ``subject_id``, ``subject_label`` (primary name),
        ``object_label`` (synonym/previous name).

        Args:
            output_path: If given, the DataFrame is also written as a TSV file.

        Returns:
            :class:`pandas.DataFrame` with label mapping rows.
        """
        import pandas as pd

        rows = [
            {
                "subject_id": str(getattr(m, "subject_id", "") or ""),
                "subject_label": str(getattr(m, "subject_label", "") or ""),
                "object_label": str(getattr(m, "object_label", "") or ""),
            }
            for m in (self.mappings or [])
            if getattr(m, "subject_label", None) or getattr(m, "object_label", None)
        ]
        df = pd.DataFrame(rows, columns=["subject_id", "subject_label", "object_label"])

        if output_path is not None:
            path = self._resolve_path(output_path, "_name2synonym.tsv")
            path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(path, sep="\t", index=False)

        return df

    def save(
        self,
        fmt: str,
        output_path: Path | str | None = None,
        **kwargs: object,
    ) -> Path:
        """Write to any supported format by name.

        Formats: ``"sssom"``, ``"rdf"``, ``"json"``, ``"owl"``,
        ``"symbol_sec2pri"`` (``"symbol2prev"`` is a deprecated alias),
        ``"pri_symbols"``, ``"name2synonym"``.

        Args:
            fmt: Format key (see above).
            output_path: Destination path. Auto-generated if ``None``.
            **kwargs: Forwarded to the format-specific writer.

        Returns:
            Path to the written file.

        Raises:
            ValueError: For unknown format keys.
        """
        shared = self._save_shared(fmt, output_path, **kwargs)
        if shared is not None:
            return shared

        if fmt in ("symbol_sec2pri", "symbol2prev"):
            self.to_symbol_sec2pri(output_path)
            return self._resolve_path(output_path, "_symbol_sec2pri.tsv")

        if fmt == "pri_symbols":
            self.to_pri_symbols(output_path)
            return self._resolve_path(output_path, "_pri_symbols.txt")

        if fmt == "name2synonym":
            self.to_name2synonym(output_path)
            return self._resolve_path(output_path, "_name2synonym.tsv")

        raise ValueError(
            f"Unknown format {fmt!r}. Choose from: "
            "json, name2synonym, owl, pri_symbols, rdf, sssom, symbol_sec2pri"
        )


class BaseParser(ABC):
    """Abstract base class for all datasource parsers.

    Each parser is responsible for reading files from a specific datasource
    and extracting secondary-to-primary identifier Mapping Sets.
    """

    # To be overridden by subclasses
    datasource_name: str = ""
    default_source_url: str = ""
    _config: DatasourceConfig | None = None

    def __init__(
        self,
        version: str | None = None,
        show_progress: bool = True,
        config_name: str | None = None,
    ):
        """Initialize the parser.

        Args:
            version: Version/release identifier for the datasource.
            show_progress: Whether to show progress bars during parsing.
            config_name: Name of config file to load (defaults to class name).
        """
        self.version = version
        self.show_progress = show_progress

        # Load config from YAML if available
        if config_name:
            self._config = get_datasource_config(config_name)
        elif self.datasource_name:
            try:
                self._config = get_datasource_config(self.datasource_name.lower())
            except FileNotFoundError:
                self._config = None

    @property
    def config(self) -> DatasourceConfig | None:
        """Get the loaded configuration."""
        return self._config

    def get_download_url(self, key: str) -> str | None:
        """Get a download URL from config by key."""
        if self._config:
            return self._config.download_urls.get(key)
        return None

    def get_curie_map(self) -> dict[str, str]:
        """Get the CURIE map from config."""
        if self._config and self._config.mappingset_metadata:
            result: dict[str, str] = self._config.mappingset_metadata.get("curie_map", {})
            return result
        return {}

    def get_mappingset_metadata(self) -> dict[str, Any]:
        """Get mapping set metadata from config."""
        if self._config:
            result: dict[str, Any] = self._config.mappingset_metadata
            return result
        return {}

    def get_mapping_metadata(self) -> dict[str, Any]:
        """Get mapping metadata from config."""
        if self._config:
            result: dict[str, Any] = self._config.mapping_metadata
            return result
        return {}

    def load_metadata(self, yaml_path: str) -> dict[str, Any]:
        """Load metadata from a YAML config file."""
        with open(yaml_path, encoding="utf-8") as f:
            result: dict[str, Any] = yaml.safe_load(f)
            return result

    def apply_metadata_to_mappingset(
        self,
        mappingset: MappingSet,
        metadata: dict[str, Any],
    ) -> None:
        """Apply metadata to a MappingSet and its Mappings."""
        # Set MappingSet fields
        for key, value in metadata.get("mappingset", {}).items():
            if hasattr(mappingset, key) and value is not None:
                setattr(mappingset, key, value)
        # Set Mapping fields
        if hasattr(mappingset, "mappings") and mappingset.mappings:
            for mapping in mappingset.mappings:
                for key, value in metadata.get("mapping", {}).items():
                    if hasattr(mapping, key) and value is not None:
                        setattr(mapping, key, value)

    def _extract_version_from_file(self, file_path: Path) -> str | None:
        """Extract a version string embedded in a data file's header.

        Override in subclasses where the source file contains release
        metadata (e.g. ``Release: 2026_01`` in UniProt flat files).

        Args:
            file_path: Path to the data file to inspect.

        Returns:
            Version string, or ``None`` if not found.
        """
        return None

    def _resolve_version(self, file_path: Path | None = None) -> str:
        """Resolve the dataset version to use for source version fields.

        Resolution order:
        1. ``self.version`` if already set explicitly.
        2. Version extracted from file header via ``_extract_version_from_file``.
        3. ISO date or release token found in the filename stem
           (e.g. ``withdrawn_2026-04-07.txt`` -> ``2026-04-07``,
           ``chebi_245.sdf`` -> ``245``).
        4. File modification date (ISO-8601) when a path is provided.
        5. Today's date as a last resort.

        Sets ``self.version`` to the resolved value so that
        ``create_mapping_set`` picks it up for ``subject_source_version`` /
        ``object_source_version`` automatically.

        Args:
            file_path: Optional path to the primary input file.

        Returns:
            Resolved version string.
        """
        if self.version:
            return self.version

        if file_path is not None:
            file_path = Path(file_path)
            # 1. Try header-embedded version (parser-specific override)
            extracted = self._extract_version_from_file(file_path)
            if extracted:
                self.version = extracted
                return self.version
            # 2. Try ISO date (YYYY-MM-DD) in the filename stem
            iso_match = re.search(r"\d{4}-\d{2}-\d{2}", file_path.stem)
            if iso_match:
                self.version = iso_match.group(0)
                return self.version
            # 3. Try a plain numeric/semver token in the filename stem
            #    e.g. "chebi_245" -> "245", "gene_history_v2" -> "2"
            num_match = re.search(r"(?<![.\d])(\d{3,})(?![.\d])", file_path.stem)
            if num_match:
                self.version = num_match.group(1)
                return self.version
            # 4. Fall back to file modification time
            try:
                mtime = file_path.stat().st_mtime
                self.version = date.fromtimestamp(mtime).isoformat()
                return self.version
            except OSError:
                pass

        self.version = date.today().isoformat()
        return self.version

    @abstractmethod
    def parse(self, input_path: Path | str | None) -> MappingSet:
        """Parse the input file(s) and return a MappingSet.

        Args:
            input_path: Path to the input file or directory.

        Returns:
            A MappingSet containing all extracted mappings.
        """

    def _progress(
        self,
        iterable: Iterable[_T],
        desc: str | None = None,
        total: int | None = None,
    ) -> Iterable[_T]:
        """Wrap an iterable with a progress bar if enabled.

        Args:
            iterable: The iterable to wrap.
            desc: Description for the progress bar.
            total: Total number of items (if known).

        Returns:
            The iterable, optionally wrapped in tqdm.
        """
        if self.show_progress:
            return cast(Iterable[_T], tqdm(iterable, desc=desc, total=total))
        return iterable

    def _build_mappings(
        self,
        rows: Iterable[dict[str, Any]],
        fixed_fields: dict[str, Any] | None = None,
        *,
        desc: str = "Building mappings",
        total: int | None = None,
    ) -> list[Mapping]:
        """Build SSSOM Mapping objects from row dicts.

        Automatically injects mapping-level fields from the parser's
        config metadata (e.g. ``confidence``) unless the caller already
        provides them in ``fixed_fields`` or individual row dicts.

        Args:
            rows: Per-row fields as dicts (subject_id, object_id, etc.).
            fixed_fields: Fields shared by all rows (predicate_id, license, etc.).
            desc: Progress bar description.
            total: Total count for the progress bar.

        Returns:
            List of Mapping objects.
        """
        _auto_fields = ("confidence",)
        m_meta = self.get_mapping_metadata()
        auto: dict[str, Any] = {
            k: m_meta[k] for k in _auto_fields if k in m_meta and m_meta[k] is not None
        }

        # Build base: auto fields < fixed_fields (fixed wins over auto)
        base: dict[str, Any] = {**auto, **(fixed_fields or {})}

        if base:
            merged: Iterable[dict[str, Any]] = ({**base, **row} for row in rows)
        else:
            merged = rows
        return [Mapping(**row) for row in self._progress(merged, desc=desc, total=total)]

    def _build_comment(
        self,
        base_comment: str,
        additional: str | None = None,
    ) -> str:
        """Build a comment string with version information.

        Args:
            base_comment: The base comment text.
            additional: Additional text to append.

        Returns:
            The complete comment string.
        """
        parts = [base_comment] if base_comment else []
        if additional:
            parts.append(additional)
        if self.version:
            parts.append(f"Release: {self.version}.")
        return " ".join(parts)

    def _find_merged_column(
        self,
        columns: list[str],
        merged_info_patterns: list[str],
    ) -> str | None:
        """Find the merged info column regardless of naming variant."""
        normalized_patterns = [p.lower() for p in merged_info_patterns]
        for col in columns:
            normalized = self._normalize_column_name(col)
            if normalized in normalized_patterns:
                return col
            # Also check for partial match on key identifying part
            if "merged_into_report" in normalized:
                return col
        return None

    @staticmethod
    def _normalize_column_name(col: str) -> str:
        """Normalize column name for case-insensitive matching."""
        return col.lower().strip()

    @staticmethod
    def _find_column(columns: list[str], name: str) -> str | None:
        """Find column by case-insensitive name."""
        lower_name = name.lower()
        for col in columns:
            if col.lower() == lower_name:
                return col
        return None

    @staticmethod
    def normalize_withdrawn_id(subject_id: str | None) -> str:
        """Normalize a primary ID, converting empty/null to withdrawn.

        Args:
            subject_id: The raw primary identifier from the source file.

        Returns:
            The normalized primary ID, or WITHDRAWN_ENTRY for empty values.
        """
        if not subject_id or subject_id in ("-", ""):
            return WITHDRAWN_ENTRY
        return subject_id

    @staticmethod
    def is_withdrawn(identifier: str) -> bool:
        """Return if is withdrawn."""
        return WITHDRAWN_ENTRY == identifier

    @staticmethod
    def _split_symbols(symbols_str: str) -> list[str]:
        """Split a pipe-separated string of symbols."""
        if not symbols_str:
            return []
        return [s.strip() for s in symbols_str.split("|") if s.strip()]

    @staticmethod
    def is_withdrawn_primary(subject_id: str) -> bool:
        """Check if a primary ID represents a withdrawn/deleted entry.

        Args:
            subject_id: The primary identifier to check.

        Returns:
            True if the primary ID indicates a withdrawn entry.
        """
        return subject_id == WITHDRAWN_ENTRY

    @staticmethod
    def _parse_merged_info(merged_str: str) -> tuple[str, str] | None:
        """Parse merged_into_report to extract hgnc_id and symbol.

        Returns (hgnc_id, symbol) or None if parsing fails.
        """
        if not merged_str or merged_str == "":
            return None
        # Try pipe separator first then slash
        if "|" in merged_str:
            parts = merged_str.split("|")
        else:
            parts = merged_str.split("/")
        if len(parts) >= 2:
            return (parts[0].strip(), parts[1].strip())
        return None

    def create_mapping_set(
        self,
        mappings: list[Mapping],
        mapping_type: str = "id",
    ) -> Sec2PriMappingSet:
        """Create an IdMappingSet or LabelMappingSet with config metadata.

        Common factory method for creating mapping sets with
        all SSSOM metadata populated from the YAML config. It also computes
        cardinalities for mappings.

        Args:
            mappings: List of SSSOM Mapping objects.
            mapping_type: "id" for IdMappingSet (cardinality by ID),
                         "label" for LabelMappingSet (cardinality by label).

        Returns:
            MappingSet with computed cardinalities.
        """
        import curies as _curies

        ms_meta = self.get_mappingset_metadata()
        curie_map = self.get_curie_map()

        # Build a converter from the curie_map
        converter = _curies.Converter.from_prefix_map(curie_map)

        def _compress(val: Any) -> Any:
            """Compress a URI string or list of URI strings to CURIEs."""
            if isinstance(val, str):
                return converter.compress(val) or val
            if isinstance(val, list):
                return [converter.compress(v) or v if isinstance(v, str) else v for v in val]
            return val

        # Choose the appropriate MappingSet class
        if mapping_type == "label":
            mapping_set_class = LabelMappingSet
        else:
            mapping_set_class = IdMappingSet

        # Build description with version if available
        description = ms_meta.get("mapping_set_description", "")
        if self.version and description:
            description = f"{description} Version: {self.version}."

        # Create the mapping set with all SSSOM metadata
        mapping_set = mapping_set_class(
            mappings=mappings,
            curie_map=curie_map,
            mapping_set_id=ms_meta.get("mapping_set_id"),
            mapping_set_version=self.version,
            mapping_set_title=ms_meta.get("mapping_set_title"),
            mapping_set_description=description or None,
            creator_id=_compress(ms_meta.get("creator_id")),
            creator_label=ms_meta.get("creator_label"),
            comment=ms_meta.get("comment"),
            license=_compress(ms_meta.get("license")),
            subject_source=ms_meta.get("subject_source"),
            subject_source_version=self.version,
            object_source=ms_meta.get("object_source"),
            object_source_version=self.version,
            mapping_provider=_compress(ms_meta.get("mapping_provider")),
            mapping_tool=_compress(ms_meta.get("mapping_tool")),
            mapping_tool_version=VERSION,
            mapping_date=date.today().isoformat(),
            see_also=_compress(ms_meta.get("see_also")),
            issue_tracker=_compress(ms_meta.get("issue_tracker")),
            subject_preprocessing=_compress(ms_meta.get("subject_preprocessing")),
            object_preprocessing=_compress(ms_meta.get("object_preprocessing")),
        )

        # Compute cardinalities using Polars vectorised counting
        mapping_set.compute_cardinalities()

        return mapping_set


__all__ = [
    "WITHDRAWN_ENTRY",
    "WITHDRAWN_ENTRY_LABEL",
    "BaseDownloader",
    "BaseParser",
    "DatasourceConfig",
    "IdMappingSet",
    "LabelMappingSet",
    "Sec2PriMappingSet",
    "get_datasource_config",
    "load_config",
]
