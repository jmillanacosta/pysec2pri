"""Base parser class for all datasource parsers."""

from __future__ import annotations

from importlib import resources as _importlib_resources
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

from mapkgsutils.parsers.base import (
    WITHDRAWN_ENTRY,
    WITHDRAWN_ENTRY_LABEL,
    AmbiguousMappingSet,
    BaseMappingSet,
    DatasourceConfig,
    DistributionEra,
    XrefSource,
    _cmp_versions,
    get_datasource_config,
    load_config,
    product_slug_values,
)
from mapkgsutils.parsers.base import BaseDownloader as _MapkgBaseDownloader
from mapkgsutils.parsers.base import BaseParser as _MapkgBaseParser

from pysec2pri.version import VERSION

if TYPE_CHECKING:
    import pandas as pd

# Config directory path

CONFIG_DIR = Path(_importlib_resources.files("pysec2pri.config"))  # type: ignore[arg-type]


#: Species selector accepted by every species-aware method: skip taxon
#: filtering entirely and process every organism in the file together.
ALL_SPECIES = "all"


class BaseDownloader(_MapkgBaseDownloader):
    """Abstract base class for pysec2pri datasource downloaders."""

    config_package = "pysec2pri.config"


class IdMappingSet(BaseMappingSet):
    """Mapping set for ID-based (secondary to primary identifier) mappings."""

    def compute_cardinalities(self) -> None:
        """Compute cardinalities using subject_id and object_id fields."""
        self._compute_cardinalities(on="id")

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
            from pysec2pri.exports import write_pri_ids

            write_pri_ids(self, self._resolve_path(output_path, "_pri_ids.txt"))

        return ids

    def to_sec2pri(self, output_path: Path | str | None = None) -> pd.DataFrame:
        """Return a ``DataFrame`` of secondary to primary ID mappings.

        Columns: ``primary_id`` (object_id), ``secondary_id`` (subject_id),
        ``predicate_id``, ``mapping_cardinality``.

        Args:
            output_path: If given, the DataFrame is also written as a TSV file.

        Returns:
            :class:`pandas.DataFrame` with one row per ID mapping.
        """
        from pysec2pri.exports import _sec2pri_frame, _write_tsv

        df = _sec2pri_frame(self)
        if output_path is not None:
            _write_tsv(df, self._resolve_path(output_path, "_sec2pri.tsv"))
        return df

    def to_secondary(self, output_path: Path | str | None = None) -> list[str]:
        """Return a sorted list of unique secondary IDs, optionally writing to TXT.

        Args:
            output_path: If given, the IDs are also written one-per-line to a
                text file.

        Returns:
            Sorted list of unique secondary ID strings.
        """
        from pysec2pri.exports import _secondary_frame, _write_tsv

        df = _secondary_frame(self)
        if output_path is not None:
            _write_tsv(df, self._resolve_path(output_path, "_secondary_ids.txt"), header=False)
        return [x for x in df["secondary_id"] if isinstance(x, str)]

    def save(
        self,
        fmt: str,
        output_path: Path | str | None = None,
        **kwargs: object,
    ) -> Path:
        """Write to any supported format by name.

        Formats: ``"sssom"``, ``"rdf"``, ``"json"``, ``"owl"``,
        ``"sec2pri"``, ``"pri_ids"``, ``"secondary"``.

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
            from pysec2pri.exports import write_sec2pri

            return write_sec2pri(self, self._resolve_path(output_path, "_sec2pri.tsv"))

        if fmt == "pri_ids":
            from pysec2pri.exports import write_pri_ids

            return write_pri_ids(self, self._resolve_path(output_path, "_pri_ids.txt"))

        if fmt == "secondary":
            from pysec2pri.exports import write_secondary

            write_secondary(self, self._resolve_path(output_path, "_secondary_ids.txt"))
            return self._resolve_path(output_path, "_secondary_ids.txt")

        raise ValueError(
            f"Unknown format {fmt!r}."
            "Choose from: json, owl, pri_ids, rdf, sec2pri, secondary, sssom"
        )


class LabelMappingSet(BaseMappingSet):
    """Mapping set for label-based (previous/alias label to current label) mappings."""

    _ambiguity_mode: ClassVar[str] = "label"

    def compute_cardinalities(self) -> None:
        """Compute cardinalities using subject_label and object_label."""
        self._compute_cardinalities(on="label")

    def to_label_sec2pri(self, output_path: Path | str | None = None) -> pd.DataFrame:
        """Return a ``DataFrame`` of previous/alias label to current label mappings.

        Columns: ``secondary_id`` (subject_id), ``secondary_label``
        (subject_label: alias or previous label), ``primary_id`` (object_id),
        ``primary_label`` (object_label: current approved label),
        ``predicate_id``, ``mapping_cardinality``.

        Args:
            output_path: If given, the DataFrame is also written as a TSV file.

        Returns:
            :class:`pandas.DataFrame` with one row per label mapping.
        """
        from pysec2pri.exports import _label_sec2pri_frame, _write_tsv

        df = _label_sec2pri_frame(self)
        if output_path is not None:
            _write_tsv(df, self._resolve_path(output_path, "_label_sec2pri.tsv"))
        return df

    def to_label2prev(self, output_path: Path | str | None = None) -> pd.DataFrame:
        """Return a ``DataFrame`` of label to previous (deprecated) label mappings.

        Columns: ``primary_id``, ``primary_label``, ``previous_label``,
        ``mapping_cardinality``.

        Only ``IAO:0100001`` (``"term replaced by"``) rows are included;
        synonym rows (``oboInOwl:hasExactSynonym``) belong in
        :meth:`to_name2synonym`, not here.

        Args:
            output_path: If given, the DataFrame is also written as a TSV file.

        Returns:
            :class:`pandas.DataFrame` with deprecation-only label mapping rows.
        """
        from pysec2pri.exports import _label2prev_frame, _write_tsv

        df = _label2prev_frame(self)
        if output_path is not None:
            _write_tsv(df, self._resolve_path(output_path, "_label2prev.tsv"))
        return df

    def to_pri_labels(self, output_path: Path | str | None = None) -> list[tuple[str, str]]:
        r"""Return a sorted list of unique ``(primary_id, label)`` pairs.

        When ``_primary_labels`` is populated (e.g. from the HGNC complete set)
        that dict is used.  Otherwise pairs are derived from the unique
        ``(object_id, object_label)`` values in the mappings.

        Args:
            output_path: If given, the pairs are also written as a two-column
                TSV file (``id\\tlabel``).

        Returns:
            Sorted list of ``(primary_id, label)`` tuples.
        """
        private: dict[str, set[str]] = (
            object.__getattribute__(self, "_primary_labels")
            if hasattr(self, "_primary_labels")
            else {}
        )

        if private:
            # Flatten dict[label, set[id]] : sorted list of (id, label) pairs
            pairs: list[tuple[str, str]] = sorted(
                (pri_id, label) for label, pri_ids in private.items() for pri_id in pri_ids
            )
        else:
            pairs = sorted(
                {
                    (
                        str(getattr(m, "object_id", None) or ""),
                        str(getattr(m, "object_label", None) or ""),
                    )
                    for m in (self.mappings or [])
                }
                - {("", "")}
            )
        if output_path is not None:
            from pysec2pri.exports import write_pri_labels

            write_pri_labels(self, self._resolve_path(output_path, "_pri_labels.txt"))

        return pairs

    def to_name2synonym(self, output_path: Path | str | None = None) -> pd.DataFrame:
        """Return a name to synonym ``DataFrame``, optionally writing to TSV.

        Columns: ``primary_id``, ``name`` (primary / canonical name),
        ``synonym`` (secondary / alternative name).

        Only ``oboInOwl:hasExactSynonym`` rows are included.  Rows with
        ``IAO:0100001`` (``"term replaced by"``) are deprecation mappings and
        belong in the ``label_sec2pri`` output, not here.

        The direction follows the sec:pri structure, where the secondary
        (synonym/alternative) term is the subject and the primary (canonical)
        term is the object.

        Args:
            output_path: If given, the DataFrame is also written as a TSV file.

        Returns:
            :class:`pandas.DataFrame` with synonym-only label mapping rows.
        """
        from pysec2pri.exports import _name2synonym_frame, _write_tsv

        df = _name2synonym_frame(self)
        if output_path is not None:
            _write_tsv(df, self._resolve_path(output_path, "_name2synonym.tsv"))
        return df

    def save(
        self,
        fmt: str,
        output_path: Path | str | None = None,
        **kwargs: object,
    ) -> Path:
        """Write to any supported format by name.

        Formats: ``"sssom"``, ``"rdf"``, ``"json"``, ``"owl"``,
        ``"label_sec2pri"``, ``"label2prev"``, ``"pri_labels"``,
        ``"name2synonym"``.

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

        if fmt == "label_sec2pri":
            from pysec2pri.exports import write_label_sec2pri

            return write_label_sec2pri(self, self._resolve_path(output_path, "_label_sec2pri.tsv"))

        if fmt == "label2prev":
            from pysec2pri.exports import write_label2prev

            return write_label2prev(self, self._resolve_path(output_path, "_label2prev.tsv"))

        if fmt == "pri_labels":
            from pysec2pri.exports import write_pri_labels

            return write_pri_labels(self, self._resolve_path(output_path, "_pri_labels.txt"))

        if fmt == "name2synonym":
            from pysec2pri.exports import write_name2synonym

            return write_name2synonym(self, self._resolve_path(output_path, "_name2synonym.tsv"))

        raise ValueError(
            f"Unknown format {fmt!r}. Choose from: "
            "json, label2prev, name2synonym, owl, pri_labels, rdf, sssom, label_sec2pri"
        )


class BaseParser(_MapkgBaseParser):
    """Abstract base class for all pysec2pri datasource parsers."""

    config_package = "pysec2pri.config"
    mapping_set_classes: ClassVar[dict[str, type[BaseMappingSet]]] = {
        "id": IdMappingSet,
        "label": LabelMappingSet,
    }
    mapping_tool_version = VERSION

    def _source_version(self) -> str | None:
        """Return the source_version."""
        version = super()._source_version()
        if not version or self._config is None:
            return version
        era = self._config.era_for(str(version))
        template = (getattr(era, "version_iri", None) if era else None) or getattr(
            self._config, "version_iri", None
        )
        if not template:
            return version
        return str(template).replace("{version}", str(version))

    def create_mapping_set(self, *args: Any, **kwargs: Any) -> BaseMappingSet:
        """Create a mapping set, tagging it with the run's taxon."""
        species = getattr(self, "species", None)
        if species is not None and str(species) != ALL_SPECIES:
            kwargs.setdefault("extension_metadata", {"taxon": f"NCBITaxon:{species}"})
        return super().create_mapping_set(*args, **kwargs)


__all__ = [
    "ALL_SPECIES",
    "CONFIG_DIR",
    "WITHDRAWN_ENTRY",
    "WITHDRAWN_ENTRY_LABEL",
    "AmbiguousMappingSet",
    "BaseDownloader",
    "BaseMappingSet",
    "BaseParser",
    "DatasourceConfig",
    "DistributionEra",
    "IdMappingSet",
    "LabelMappingSet",
    "XrefSource",
    "_cmp_versions",
    "get_datasource_config",
    "load_config",
    "product_slug_values",
]
