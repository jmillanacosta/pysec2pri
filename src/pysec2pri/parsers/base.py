"""Base parser class for all datasource parsers."""

from __future__ import annotations

from importlib import resources as _importlib_resources
from pathlib import Path
from typing import TYPE_CHECKING, ClassVar

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
)
from mapkgsutils.parsers.base import BaseDownloader as _MapkgBaseDownloader
from mapkgsutils.parsers.base import BaseParser as _MapkgBaseParser

from pysec2pri.version import VERSION

if TYPE_CHECKING:
    import pandas as pd

# Config directory path

CONFIG_DIR = Path(_importlib_resources.files("pysec2pri.config"))  # type: ignore[arg-type]


class BaseDownloader(_MapkgBaseDownloader):
    """Abstract base class for pysec2pri datasource downloaders."""

    config_package = "pysec2pri.config"


class IdMappingSet(BaseMappingSet):
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
        import pandas as pd

        rows = [
            {
                "secondary_id": str(getattr(m, "subject_id", "") or ""),
                "secondary_label": str(getattr(m, "subject_label", "") or ""),
                "primary_id": str(getattr(m, "object_id", "") or ""),
                "primary_label": str(getattr(m, "object_label", "") or ""),
                "predicate_id": str(getattr(m, "predicate_id", "") or ""),
                "mapping_cardinality": str(getattr(m, "mapping_cardinality", "") or ""),
            }
            for m in (self.mappings or [])
        ]
        df = pd.DataFrame(
            rows,
            columns=[
                "secondary_id",
                "secondary_label",
                "primary_id",
                "primary_label",
                "predicate_id",
                "mapping_cardinality",
            ],
        )

        if output_path is not None:
            path = self._resolve_path(output_path, "_label_sec2pri.tsv")
            path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(path, sep="\t", index=False)

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
            path = self._resolve_path(output_path, "_pri_labels.txt")
            path.parent.mkdir(parents=True, exist_ok=True)
            text = "\n".join(f"{pri_id}\t{label}" for pri_id, label in pairs)
            path.write_text("id\tlabel\n" + text + "\n", encoding="utf-8")

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
        import pandas as pd

        rows = [
            {
                "primary_id": str(getattr(m, "object_id", "") or ""),
                "name": str(getattr(m, "object_label", "") or ""),
                "synonym": str(getattr(m, "subject_label", "") or ""),
            }
            for m in (self.mappings or [])
            if getattr(m, "predicate_id", None) == "oboInOwl:hasExactSynonym"
            and (getattr(m, "subject_label", None) or getattr(m, "object_label", None))
        ]
        df = pd.DataFrame(rows, columns=["primary_id", "name", "synonym"])

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
        ``"label_sec2pri"`` (``"label2prev"`` is a deprecated alias),
        ``"pri_labels"``, ``"name2synonym"``.

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

        if fmt in ("label_sec2pri", "label2prev"):
            self.to_label_sec2pri(output_path)
            return self._resolve_path(output_path, "_label_sec2pri.tsv")

        if fmt == "pri_labels":
            self.to_pri_labels(output_path)
            return self._resolve_path(output_path, "_pri_labels.txt")

        if fmt == "name2synonym":
            self.to_name2synonym(output_path)
            return self._resolve_path(output_path, "_name2synonym.tsv")

        raise ValueError(
            f"Unknown format {fmt!r}. Choose from: "
            "json, name2synonym, owl, pri_labels, rdf, sssom, label_sec2pri"
        )


class BaseParser(_MapkgBaseParser):
    """Abstract base class for all pysec2pri datasource parsers."""

    config_package = "pysec2pri.config"
    mapping_set_classes: ClassVar[dict[str, type[BaseMappingSet]]] = {
        "id": IdMappingSet,
        "label": LabelMappingSet,
    }
    mapping_tool_version = VERSION


__all__ = [
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
]
