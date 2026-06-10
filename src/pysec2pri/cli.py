"""Command-line interface for pysec2pri.

Commands are registered automatically from the config files: every
(config_id, kind) pair whose formats list is non-empty gets a subcommand
inside the corresponding datasource group.  Per-datasource extra options
(--subset, --tax-id, --entity-type ...) are declared once in ``_EXTRA_OPTS``
and forwarded as **kwargs to the matching ``generate_*`` function.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any

import click

from pysec2pri.parsers.base import get_datasource_config

if TYPE_CHECKING:
    from pysec2pri.parsers.base import Sec2PriMappingSet

warnings.filterwarnings("ignore", category=FutureWarning, module="sssom")

PathType = click.Path(path_type=Path)  # type: ignore[type-var]
ExistingPathType = click.Path(exists=True, path_type=Path)  # type: ignore[type-var]


# Reusable option decorators


def _opt(*names: str, **kwargs: Any) -> Callable[..., Any]:
    """One-liner Click option factory."""
    return lambda fn: click.option(*names, **kwargs)(fn)


_opt_output = _opt("-o", "--output", type=PathType, help="Output file or directory.")
_opt_version = _opt("--version", "data_version", default=None, help="Datasource release version.")
_opt_no_progress = _opt(
    "--no-progress", is_flag=True, default=False, help="Suppress progress bars."
)

# Datasource-specific options.
# The kwarg name injected by Click MUST match the generate_* function parameter name.
_opt_subset = _opt(
    "--subset",
    default="3star",
    show_default=True,
    type=click.Choice(["3star", "complete"]),
    help="Compound subset.",
)
_opt_tax_id = _opt("--tax-id", default="9606", show_default=True, help="NCBI taxonomy ID.")
_opt_entity_type = _opt(
    "--entity-type",
    default=None,
    type=click.Choice(["metabolites", "chemicals", "genes", "proteins"]),
    help="Wikidata entity type to query. Queries all if omitted.",
)
_opt_test_subset = _opt(
    "--test-subset",
    is_flag=True,
    default=False,
    help="Use test SPARQL queries (LIMIT 10).",
)
_opt_delac_file = _opt(
    "--delac-file",
    type=click.Path(exists=True),
    default=None,
    help="Path to delac_sp.txt (UniProt).",
)


# Shared CLI helpers


def _version_base(ms: object, data_version: str | None, prefix: str) -> tuple[str, str]:
    """Return ``(version, base_name)`` for output filenames."""
    version = getattr(ms, "mapping_set_version", None) or data_version or ""
    base = f"{prefix}{'_' + version if version else ''}"
    return version, base


def _emit(ms: object, fmt: str, output: Path | None, base: str) -> None:
    """Persist *ms* and echo a brief summary."""
    from pysec2pri.api import save

    path = save(ms, fmt, output, base_name=base)  # type: ignore[arg-type]
    if fmt == "all":
        click.echo(f"Wrote all formats -> {path}/")
    elif fmt == "pri_ids":
        n = len(getattr(ms, "_primary_ids", None) or set())
        click.echo(f"Wrote {n} primary IDs -> {path}")
    elif fmt == "pri_labels":
        n = len(getattr(ms, "_primary_labels", None) or set())
        click.echo(f"Wrote {n} primary labels -> {path}")
    else:
        n = len(getattr(ms, "mappings", None) or [])
        click.echo(f"Wrote {n} mappings -> {path}")


def _resolve_and_print(
    input_file: Path,
    ms: Sec2PriMappingSet,
    col_syn_pairs: list[tuple[str, str | None]],
    output_path: Path | None,
    suffix: str,
    sep: str | None,
    label_ms: Sec2PriMappingSet | None = None,
    *,
    mode: str = "ids",
) -> None:
    """Read *input_file* and resolve each (column, synonym_col) pair."""
    import pandas as pd

    from pysec2pri.update_ids import update_ids, update_labels

    inferred_sep = "\t" if input_file.suffix.lower() == ".tsv" else ","
    read_sep = sep if sep is not None else inferred_sep
    df = pd.read_csv(input_file, sep=read_sep, dtype=str)
    assert isinstance(df, pd.DataFrame)
    result: pd.DataFrame = df.copy()
    click.echo(f"Resolving column(s): {', '.join(c for c, _ in col_syn_pairs)}")
    for col, syn_col in col_syn_pairs:
        hint = f" (hints from {syn_col!r})" if syn_col else ""
        click.echo(f"  {col!r}{hint}")
        partial: pd.DataFrame = (
            update_ids(
                result, ms, at=col, suffix=suffix, synonyms=syn_col, label_mapping_set=label_ms
            )
            if mode == "ids"
            else update_labels(result, ms, at=col, suffix=suffix, synonyms=syn_col)
        )
        for new_col in partial.columns:
            if new_col not in result.columns:
                result[new_col] = partial[new_col]
    if output_path is not None:
        out_sep = "\t" if output_path.suffix.lower() == ".tsv" else ","
        result.to_csv(output_path, sep=out_sep, index=False)
        click.echo(f"Wrote {len(result)} rows -> {output_path}")
    else:
        out_sep = "\t" if (sep or "\t") == "\t" else ","
        click.echo(result.to_csv(sep=out_sep, index=False))


# Dynamic generate-command factory


def _make_generate_cmd(
    config_id: str,
    kind: str,
    generate_fn: Callable[..., Any],
    extra_opts: list[Callable[..., Any]],
) -> Callable[..., Any]:
    """Return a decorated (but not yet click.command-wrapped) callable."""
    formats = get_datasource_config(config_id).formats_for(kind)

    def _cmd(
        input_file: Path | None,
        output: Path | None,
        data_version: str | None,
        output_format: str,
        no_progress: bool,
        **extra_kwargs: Any,
    ) -> None:
        ms = generate_fn(
            input_path=input_file,
            version=data_version,
            show_progress=not no_progress,
            **extra_kwargs,
        )
        _, base = _version_base(ms, data_version, f"{config_id}_{kind}")
        _emit(ms, output_format, output, base)

    # Apply Click decorators
    decorators: list[Callable[..., Any]] = [
        click.option(
            "--input-file",
            "input_file",
            type=ExistingPathType,
            default=None,
            show_default=True,
        ),
        _opt_output,
        _opt_version,
        click.option(
            "--format",
            "output_format",
            default="sssom",
            show_default=True,
            type=click.Choice(formats),
            help="Output format.",
        ),
        _opt_no_progress,
        *extra_opts,
    ]
    for dec in reversed(decorators):
        _cmd = dec(_cmd)

    _cmd.__name__ = kind
    return _cmd


# Datasource registry


def _build_registry() -> dict[tuple[str, str], tuple[Callable[..., Any], list[Callable[..., Any]]]]:
    """Return (config_id, kind) -> (generate_fn, [extra_opt_decorators])."""
    from pysec2pri import api

    return {
        ("chebi", "ids"): (api.generate_chebi_ids, [_opt_subset]),
        ("chebi", "labels"): (api.generate_chebi_labels, [_opt_subset]),
        ("hgnc", "ids"): (api.generate_hgnc_ids, []),
        ("hgnc", "labels"): (api.generate_hgnc_labels, []),
        ("ncbi", "ids"): (api.generate_ncbi_ids, [_opt_tax_id]),
        ("ncbi", "labels"): (api.generate_ncbi_labels, [_opt_tax_id]),
        ("hmdb_metabolites", "ids"): (api.generate_hmdb_ids, []),
        ("hmdb_proteins", "ids"): (api.generate_hmdb_proteins_ids, []),
        ("uniprot", "ids"): (api.generate_uniprot_ids, [_opt_delac_file]),
        ("wikidata", "ids"): (api.generate_wikidata_ids, [_opt_entity_type, _opt_test_subset]),
        ("wikidata", "labels"): (
            api.generate_wikidata_labels,
            [_opt_entity_type, _opt_test_subset],
        ),
        # Add new sources here
    }


def _register_datasources(parent: click.Group) -> None:
    """Register one Click group per config_id on *parent*."""
    registry = _build_registry()
    by_config: dict[str, dict[str, tuple[Callable[..., Any], list[Callable[..., Any]]]]] = {}
    for (cfg_id, kind), (fn, opts) in registry.items():
        by_config.setdefault(cfg_id, {})[kind] = (fn, opts)

    for cfg_id, kinds in by_config.items():
        cfg = get_datasource_config(cfg_id)
        grp = click.Group(cfg_id.replace("_", "-"), help=f"{cfg.name} mappings.")
        for kind, (fn, extra_opts) in kinds.items():
            raw_cmd = _make_generate_cmd(cfg_id, kind, fn, extra_opts)
            grp.add_command(click.command(name=kind)(raw_cmd))
        parent.add_command(grp)


# Main entry point


@click.group()
@click.version_option()
def main() -> None:
    """pysec2pri -- secondary-to-primary ID and label mapping."""


_register_datasources(main)


# diff


@main.command()
@click.argument("file1", type=ExistingPathType)
@click.argument("file2", type=ExistingPathType)
@click.option("-o", "--output", type=PathType, help="Output file for diff results (TSV).")
@click.option("--show-all", is_flag=True, default=False, help="Show all differences.")
@click.option("--datasource", default="unknown", help="Datasource name for diff summary.")
def diff(file1: Path, file2: Path, output: Path | None, show_all: bool, datasource: str) -> None:
    """Compare two SSSOM mapping files and show differences."""
    from pysec2pri.api import write_diff_output
    from pysec2pri.diff import diff_sssom_files, summarize_diff

    click.echo(f"Comparing {file1.name} vs {file2.name}...")
    result = diff_sssom_files(file1, file2, datasource=datasource)
    click.echo("")
    click.echo(summarize_diff(result))
    if show_all or result.total_changes <= 50:
        _show_diff_details(result)
    elif result.total_changes > 50:
        click.echo(
            f"\n(Showing summary only -- {result.total_changes} changes. "
            "Use --show-all to see all.)"
        )
    if output:
        if result.total_changes > 0:
            write_diff_output(result, output)
            click.echo(f"\nWrote diff -> {output}")
        else:
            click.echo("\nNo differences to write.")


def _show_diff_details(result: object) -> None:
    for label, cols in [
        ("added", ("object_id", "subject_id")),
        ("removed", ("object_id", "subject_id")),
        ("changed", ("object_id", "old_subject_id", "new_subject_id")),
    ]:
        count = getattr(result, f"{label}_count", 0)
        if count:
            click.echo(f"\n=== {label.capitalize()} mappings ({count}) ===")
            click.echo("\t".join(cols))
            for row in getattr(result, label).iter_rows(named=True):
                click.echo("\t".join(str(row[c]) for c in cols))


# list-versions


@main.command("list-versions")
@click.argument(
    "datasource",
    type=click.Choice(["chebi", "hgnc", "uniprot"], case_sensitive=False),
)
def list_versions_cmd(datasource: str) -> None:
    """List available archive versions for DATASOURCE (chebi, hgnc, uniprot)."""
    from pysec2pri.api import list_versions as _list_versions

    click.echo(f"Fetching available {datasource.upper()} versions...")
    try:
        versions = _list_versions(datasource)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    if not versions:
        click.echo("No versions found.")
        return

    click.echo(f"\n{len(versions)} version(s) for {datasource.upper()}:\n")
    if datasource.lower() == "chebi":  # hacky patch for chebi
        from pysec2pri.parsers.chebi import ChEBIDownloader

        threshold = ChEBIDownloader().new_format_version or 245
        for v in versions:
            click.echo(f"  {v}" + ("  (SDF only)" if int(v) < threshold else ""))
        click.echo(f"\n  Note: versions < {threshold} are SDF only.")
    else:
        for v in versions:
            click.echo(f"  {v}")


# ambiguous


def _ambiguous_choices() -> list[str]:
    return [f"{cfg.replace('_', '-')}-{kind}" for cfg, kind in _build_registry()]


@main.command("ambiguous")
@click.argument("datasource", type=click.Choice(_ambiguous_choices()))
@_opt_output
@_opt_version
@_opt_no_progress
def ambiguous_cmd(
    datasource: str,
    output: Path | None,
    data_version: str | None,
    no_progress: bool,
) -> None:
    """Find ambiguous mappings for DATASOURCE and save as SSSOM.

    DATASOURCE format: ``<config-id>-<kind>``, e.g. ``hgnc-ids``,
    ``chebi-labels``, ``hmdb-metabolites-ids``, ``wikidata-labels``.
    """
    registry = _build_registry()
    match: tuple[str, str] | None = next(
        ((cfg, kind) for cfg, kind in registry if datasource == f"{cfg.replace('_', '-')}-{kind}"),
        None,
    )
    if match is None:
        raise click.ClickException(f"Unknown datasource: {datasource!r}")

    generate_fn, _ = registry[match]
    click.echo(f"Loading {datasource.upper()} mappings...")
    ms = generate_fn(version=data_version, show_progress=not no_progress)

    from pysec2pri.api import find_ambiguous

    click.echo("Detecting ambiguous identifiers...")
    amb = find_ambiguous(ms)
    count = len(amb.mappings or [])
    if count == 0:
        click.echo("No ambiguous identifiers found.")
        return
    click.echo(f"Found {count} ambiguous identifier(s).")
    dest = output or Path(f"{datasource}_ambiguous.sssom.tsv")
    amb.save("sssom", dest)
    click.echo(f"Wrote ambiguous mappings -> {dest}")


# all


@main.command("all")
@click.option("-o", "--output-dir", type=PathType, default=Path("."), help="Output directory.")
@click.option(
    "--datasources",
    default="chebi,hgnc,ncbi,uniprot,wikidata,hmdb_metabolites,hmdb_proteins",
    show_default=True,
    help="Comma-separated config IDs to export.",
)
def export_all(output_dir: Path, datasources: str) -> None:
    """Export all output formats for each listed datasource."""
    from pysec2pri.api import save
    from pysec2pri.download import CloudflareBlockedError

    registry = _build_registry()
    for ds in [d.strip() for d in datasources.split(",")]:
        entry = registry.get((ds, "ids"))
        if entry is None:
            click.echo(f"  Unknown datasource: {ds!r}", err=True)
            continue
        generate_fn, _ = entry
        click.echo(f"\n=== {ds.upper()} ===")
        try:
            ms = generate_fn()
            version = getattr(ms, "mapping_set_version", None)
            base = f"{ds}{'_' + version if version else ''}"
            ds_dir = output_dir / (f"{ds}_{version}" if version else ds)
            path = save(ms, "all", ds_dir, base_name=base)
            click.echo(f"  Wrote all formats -> {path}/")
        except CloudflareBlockedError as e:
            click.echo(f"  WARNING: Skipping {ds} -- {e}", err=True)
        except Exception as e:
            click.echo(f"  WARNING: Skipping {ds} -- {e}", err=True)


# update-ids

_ID_DATASOURCES = sorted(cfg for cfg, kind in _build_registry() if kind == "ids")
_LABEL_GENERATORS_FOR_IDS: dict[str, str] = {
    "chebi": "generate_chebi_labels",
    "hgnc": "generate_hgnc_labels",
    "ncbi": "generate_ncbi_labels",
}


@main.command("update-ids")
@click.argument("input_file", type=ExistingPathType)
@click.argument("datasource", type=click.Choice(_ID_DATASOURCES))
@click.option(
    "--at",
    "columns",
    required=True,
    multiple=True,
    metavar="COLUMN",
    help="Column(s) containing IDs to resolve. Repeat for multiple.",
)
@click.option("-o", "--output", "output_path", type=PathType, help="Output file (TSV or CSV).")
@click.option("--suffix", default="_primary", show_default=True, help="New-column suffix.")
@click.option("--sep", default=None, help="Delimiter (inferred from extension if omitted).")
@click.option(
    "--mapping",
    "mapping_file",
    type=ExistingPathType,
    default=None,
    help="Pre-built sec2pri TSV file (skips download).",
)
@click.option(
    "--synonyms",
    "synonyms_cols",
    default=None,
    multiple=True,
    metavar="COLUMN",
    help="Hint column paired with --at column. Repeat to pair each.",
)
@click.option(
    "--synonyms-mapping",
    "synonyms_mapping_file",
    type=ExistingPathType,
    default=None,
    help="Pre-built label/label mapping file for alias resolution.",
)
@_opt_version
@_opt_no_progress
def update_ids_cmd(
    input_file: Path,
    datasource: str,
    columns: tuple[str, ...],
    output_path: Path | None,
    suffix: str,
    sep: str | None,
    mapping_file: Path | None,
    synonyms_cols: tuple[str, ...],
    synonyms_mapping_file: Path | None,
    data_version: str | None,
    no_progress: bool,
) -> None:
    """Resolve secondary IDs in INPUT_FILE to primary IDs using DATASOURCE.

    Examples::

        pysec2pri update-ids genes.tsv hgnc --at gene_id -o out.tsv
        pysec2pri update-ids genes.tsv hgnc --at gene_id --synonyms label
    """
    from pysec2pri.api import load_label_mapping, load_mapping

    registry = _build_registry()
    if mapping_file is not None:
        click.echo(f"Loading mappings from {mapping_file}...")
        ms = load_mapping(mapping_file)
    else:
        generate_fn, _ = registry[(datasource, "ids")]
        click.echo(f"Loading {datasource.upper()} mappings...")
        ms = generate_fn(version=data_version, show_progress=not no_progress)

    label_ms = None
    if synonyms_cols:
        if synonyms_mapping_file is not None:
            click.echo(f"Loading synonym mappings from {synonyms_mapping_file}...")
            try:
                label_ms = load_label_mapping(synonyms_mapping_file)
            except Exception:
                label_ms = load_mapping(synonyms_mapping_file)
        elif datasource in _LABEL_GENERATORS_FOR_IDS:
            import importlib

            api_mod = importlib.import_module("pysec2pri.api")
            fn_name = _LABEL_GENERATORS_FOR_IDS[datasource]
            fn = getattr(api_mod, fn_name)
            click.echo(f"Loading {datasource.upper()} label mappings for alias resolution...")
            label_ms = fn(version=data_version, show_progress=not no_progress)
        else:
            click.echo(
                f"Warning: --synonyms not supported for {datasource!r} (no label mapping).",
                err=True,
            )

    padded = list(synonyms_cols) + [None] * max(0, len(columns) - len(synonyms_cols))
    _resolve_and_print(
        input_file,
        ms,
        list(zip(columns, padded[: len(columns)], strict=False)),
        output_path,
        suffix,
        sep,
        label_ms,
        mode="ids",
    )


# update-labels

_LABELS_DATASOURCES = sorted(cfg for cfg, kind in _build_registry() if kind == "labels")
_ID_GENERATORS_FOR_LABELS: dict[str, str] = {
    "hgnc": "generate_hgnc_ids",
    "ncbi": "generate_ncbi_ids",
}


@main.command("update-labels")
@click.argument("input_file", type=ExistingPathType)
@click.argument("datasource", type=click.Choice(_LABELS_DATASOURCES))
@click.option(
    "--at",
    "columns",
    required=True,
    multiple=True,
    metavar="COLUMN",
    help="Column(s) containing labels to resolve. Repeat for multiple.",
)
@click.option("-o", "--output", "output_path", type=PathType, help="Output file (TSV or CSV).")
@click.option("--suffix", default="_current", show_default=True, help="New-column suffix.")
@click.option("--sep", default=None, help="Delimiter (inferred from extension if omitted).")
@click.option(
    "--mapping",
    "mapping_file",
    type=ExistingPathType,
    default=None,
    help="Pre-built label2prev TSV file (skips download).",
)
@click.option(
    "--synonyms",
    "synonyms_cols",
    default=None,
    multiple=True,
    metavar="COLUMN",
    help="Hint column paired with --at column. Repeat to pair each.",
)
@click.option(
    "--synonyms-mapping",
    "synonyms_mapping_file",
    type=ExistingPathType,
    default=None,
    help="Pre-built mapping file for alias resolution.",
)
@_opt_tax_id
@_opt_entity_type
@_opt_subset
@_opt_version
@_opt_no_progress
def update_labels_cmd(
    input_file: Path,
    datasource: str,
    columns: tuple[str, ...],
    output_path: Path | None,
    suffix: str,
    sep: str | None,
    mapping_file: Path | None,
    synonyms_cols: tuple[str, ...],
    synonyms_mapping_file: Path | None,
    tax_id: str,
    entity_type: str | None,
    subset: str,
    data_version: str | None,
    no_progress: bool,
) -> None:
    """Resolve previous/alias labels in INPUT_FILE to current labels using DATASOURCE.

    Examples::

        pysec2pri update-labels genes.tsv hgnc --at label -o out.tsv
        pysec2pri update-labels genes.tsv hgnc --at label --mapping labels.tsv
    """
    from pysec2pri.api import (
        generate_chebi_labels,
        generate_hgnc_labels,
        generate_ncbi_labels,
        generate_wikidata_labels,
        load_label_mapping,
        load_mapping,
    )

    if mapping_file is not None:
        click.echo(f"Loading label mappings from {mapping_file}...")
        ms = load_label_mapping(mapping_file)
    else:
        click.echo(f"Loading {datasource.upper()} label mappings...")
        _label_fns: dict[str, Callable[..., Any]] = {
            "chebi": lambda: generate_chebi_labels(
                version=data_version, subset=subset, show_progress=not no_progress
            ),
            "hgnc": lambda: generate_hgnc_labels(
                version=data_version, show_progress=not no_progress
            ),
            "ncbi": lambda: generate_ncbi_labels(
                tax_id=tax_id, version=data_version, show_progress=not no_progress
            ),
            "wikidata": lambda: generate_wikidata_labels(
                entity_type=entity_type, version=data_version, show_progress=not no_progress
            ),
        }
        ms = _label_fns[datasource]()

    label_ms = None
    if synonyms_cols:
        if synonyms_mapping_file is not None:
            click.echo(f"Loading synonym mappings from {synonyms_mapping_file}...")
            try:
                label_ms = load_label_mapping(synonyms_mapping_file)
            except Exception:
                label_ms = load_mapping(synonyms_mapping_file)
        elif datasource in _ID_GENERATORS_FOR_LABELS:
            import importlib

            api_mod = importlib.import_module("pysec2pri.api")
            fn_name = _ID_GENERATORS_FOR_LABELS[datasource]
            fn = getattr(api_mod, fn_name)
            click.echo(f"Loading {datasource.upper()} ID mappings for alias resolution...")
            label_ms = fn(version=data_version, show_progress=not no_progress)
        else:
            click.echo(f"Warning: --synonyms not supported for {datasource!r}.", err=True)

    padded = list(synonyms_cols) + [None] * max(0, len(columns) - len(synonyms_cols))
    _resolve_and_print(
        input_file,
        ms,
        list(zip(columns, padded[: len(columns)], strict=False)),
        output_path,
        suffix,
        sep,
        label_ms,
        mode="labels",
    )


if __name__ == "__main__":
    main()
