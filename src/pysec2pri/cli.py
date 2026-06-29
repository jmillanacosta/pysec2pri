"""Command-line interface for pysec2pri.

Commands are registered automatically from the config files: every
(config_id, kind) pair whose formats list is non-empty gets a subcommand
inside the corresponding datasource group.  Per-datasource extra options
(--subset, --species, --entity-type ...) are declared once in ``_EXTRA_OPTS``
and forwarded as **kwargs to the matching ``generate_*`` function.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import click

from pysec2pri.parsers.base import get_datasource_config

if TYPE_CHECKING:
    from pysec2pri.parsers.base import BaseMappingSet

warnings.filterwarnings("ignore", category=FutureWarning, module="sssom")

_CONFIG_PACKAGE = "pysec2pri.config"

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


_MAX_INLINE_SPECIES = 40


def _species_choices_text(cfg_id: str, limit: int = _MAX_INLINE_SPECIES) -> str:
    """Return a ``"<taxon_id>=<label>, ..."`` list from ``<cfg_id>.yaml``'s ``species.available``.

    These are commonly-used species curated for CLI help/defaults, not a
    hard restriction -- NCBI and VGNC accept any taxon ID they have data
    for, and Ensembl resolves anything beyond this list via its own live
    species index (see ``EnsemblDownloader._resolve_species_token``).

    Truncated to *limit* entries when ``available`` is large (e.g.
    Ensembl's 270+) so ``--help`` output stays readable; the remainder is
    summarized with a pointer to the full config file instead of dumped
    inline.
    """
    cfg = get_datasource_config(cfg_id, config_package=_CONFIG_PACKAGE)
    available = (cfg.species or {}).get("available") or {}
    pairs = sorted(available.items(), key=lambda kv: str(kv[0]))
    shown = pairs[:limit]
    text = ", ".join(f"{tid}={(info or {}).get('label', tid)}" for tid, info in shown)
    remaining = len(pairs) - len(shown)
    if remaining > 0:
        text += f", and {remaining} more. See config/{cfg_id}.yaml for the full list"
    return text


def _opt_species_for(cfg_id: str) -> Callable[..., Any]:
    """Build a ``--species`` option defaulting to ``<cfg_id>.yaml``'s ``species.default``.

    Config-yaml-centric: each species-aware datasource (ensembl, ncbi, vgnc)
    declares its own ``species.default`` in its own config file (see
    ``DatasourceConfig.default_species``), rather than this module
    hardcoding a single literal for every datasource. The help text lists
    every species curated in that datasource's ``species.available``.
    """
    default = str(get_datasource_config(cfg_id, config_package=_CONFIG_PACKAGE).default_species())
    choices = _species_choices_text(cfg_id)
    help_text = "Species as NCBI taxon ID, or 'all' to process every species."
    if choices:
        help_text += f" Known: {choices}."
    if cfg_id == "ensembl":
        help_text += (
            " Other taxon IDs are resolved via Ensembl's own live species list; "
            "'all' downloads and combines all ~276 species (slow, network-heavy)."
        )
    return _opt("--species", default=default, show_default=True, help=help_text)


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

# Disambiguation-context options, shared by update-ids/update-labels.
_opt_xref = _opt(
    "--xref",
    "xref_cols",
    default=None,
    multiple=True,
    metavar="COLUMN",
    help="Column with cross-ref tokens (e.g. Ensembl IDs), paired with --at.",
)
_opt_xref_file = _opt(
    "--xref-file",
    "xref_file",
    type=ExistingPathType,
    default=None,
    help="Crosswalk table (SSSOM or plain TSV) for --xref.",
)
_opt_xref_source = _opt(
    "--xref-source",
    "xref_source",
    default=None,
    help="Use a config-declared suggested crosswalk source (auto-downloaded), e.g. 'hgnc_custom'.",
)
_opt_xref_on = _opt(
    "--xref-on",
    "xref_on",
    default=None,
    help="Which subject column of --xref-source to key on (e.g. ensembl/entrez/refseq/uniprot).",
)
_opt_xref_predicate = _opt(
    "--xref-predicate",
    "xref_predicates",
    multiple=True,
    metavar="PREDICATE",
    help="Accepted equivalence predicate(s) for xref records. Repeat for multiple. "
    "Default accepts any predicate (and unannotated records).",
)
_opt_report = _opt(
    "--report",
    "report_path",
    type=PathType,
    default=None,
    help="Write a per-decision audit log (TSV) for context-based resolutions.",
)


# Shared CLI helpers


def _resolve_xref_mapping(
    xref_file: Path | None,
    xref_source: str | None,
    xref_on: str | None,
    datasource: str,
) -> object | None:
    """Resolve the crosswalk table for --xref-file/--xref-source.

    Returns ``None`` when neither is given (caller decides whether that's an
    error, e.g. when --xref columns were supplied without a table).
    """
    from mapkgsutils.context import load_xref_mapping

    if xref_file is not None:
        return load_xref_mapping(xref_file)

    if xref_source is not None:
        if xref_on is None:
            raise click.ClickException("--xref-on is required together with --xref-source.")
        cfg = get_datasource_config(datasource, config_package=_CONFIG_PACKAGE)
        src = cfg.xref_source(xref_source)
        if src is None:
            known = ", ".join(s.id for s in cfg.xref_sources) or "(none configured)"
            raise click.ClickException(
                f"Unknown --xref-source {xref_source!r} for {datasource!r}. Known: {known}"
            )
        subject_col = src.subject_id_cols.get(xref_on)
        if subject_col is None:
            known_keys = ", ".join(sorted(src.subject_id_cols)) or "(none)"
            raise click.ClickException(
                f"Unknown --xref-on {xref_on!r} for xref-source"
                f" {xref_source!r}. Known: {known_keys}"
            )
        if src.note:
            click.echo(f"Note: {src.note}")
        click.echo(f"Downloading xref source {src.id!r}...")
        from mapkgsutils.context import download_xref_source

        return download_xref_source(src, subject_col, show_progress=False)

    return None


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


def _pad_cols(cols: tuple[str, ...], n: int) -> list[str | None]:
    """Pad/truncate *cols* to length *n*, filling missing entries with ``None``."""
    padded: list[str | None] = list(cols) + [None] * max(0, n - len(cols))
    return padded[:n]


def _report_path_for(base: Path | None, col: str, multi: bool) -> Path | None:
    """Return *base*, or a per-column variant when several --at columns share one --report."""
    if base is None or not multi:
        return base
    return base.with_name(f"{base.stem}_{col}{base.suffix}")


def _resolve_and_print(
    input_file: Path,
    ms: BaseMappingSet,
    col_specs: list[tuple[str, str | None, str | None]],
    output_path: Path | None,
    suffix: str,
    sep: str | None,
    label_ms: BaseMappingSet | None = None,
    *,
    mode: str = "ids",
    xref_mapping: object | None = None,
    xref_predicates: set[str] | None = None,
    report_path: Path | None = None,
) -> None:
    """Read *input_file* and resolve each (column, synonym_col, xref_col) triple."""
    import pandas as pd

    from pysec2pri.update_ids import update_ids, update_labels

    inferred_sep = "\t" if input_file.suffix.lower() == ".tsv" else ","
    read_sep = sep if sep is not None else inferred_sep
    df = pd.read_csv(input_file, sep=read_sep, dtype=str)
    assert isinstance(df, pd.DataFrame)
    result: pd.DataFrame = df.copy()
    click.echo(f"Resolving column(s): {', '.join(c for c, _, _ in col_specs)}")
    multi = len(col_specs) > 1
    for col, syn_col, xref_col in col_specs:
        hints = [f"hints from {syn_col!r}"] if syn_col else []
        if xref_col:
            hints.append(f"xref from {xref_col!r}")
        hint = f" ({', '.join(hints)})" if hints else ""
        click.echo(f"  {col!r}{hint}")
        col_report = _report_path_for(report_path, col, multi)
        kwargs: dict[str, Any] = {
            "at": col,
            "suffix": suffix,
            "synonyms": syn_col,
            "xref": xref_col,
            "xref_mapping": xref_mapping if xref_col else None,
            "xref_predicates": set(xref_predicates) if xref_col and xref_predicates else None,
            "report_path": col_report,
        }
        partial: pd.DataFrame = (
            update_ids(result, ms, label_mapping_set=label_ms, **kwargs)
            if mode == "ids"
            else update_labels(result, ms, **kwargs)
        )
        for new_col in partial.columns:
            if new_col not in result.columns:
                result[new_col] = partial[new_col]
        if col_report is not None:
            click.echo(f"    Wrote decision log -> {col_report}")
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
    formats = get_datasource_config(config_id, config_package=_CONFIG_PACKAGE).formats_for(kind)

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
        prefix = f"{config_id}_{kind}"
        species = extra_kwargs.get("species")
        if species is not None:
            prefix = f"{prefix}_{species}"
        _, base = _version_base(ms, data_version, prefix)
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
        ("ensembl", "ids"): (api.generate_ensembl_ids, [_opt_species_for("ensembl")]),
        ("ensembl", "labels"): (api.generate_ensembl_labels, [_opt_species_for("ensembl")]),
        ("hgnc", "ids"): (api.generate_hgnc_ids, []),
        ("hgnc", "labels"): (api.generate_hgnc_labels, []),
        ("ncbi", "ids"): (api.generate_ncbi_ids, [_opt_species_for("ncbi")]),
        ("ncbi", "labels"): (api.generate_ncbi_labels, [_opt_species_for("ncbi")]),
        ("hmdb_metabolites", "ids"): (api.generate_hmdb_ids, []),
        ("hmdb_proteins", "ids"): (api.generate_hmdb_proteins_ids, []),
        ("uniprot", "ids"): (api.generate_uniprot_ids, [_opt_delac_file]),
        ("vgnc", "ids"): (api.generate_vgnc_ids, [_opt_species_for("vgnc")]),
        ("vgnc", "labels"): (api.generate_vgnc_labels, [_opt_species_for("vgnc")]),
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
        cfg = get_datasource_config(cfg_id, config_package=_CONFIG_PACKAGE)
        group_help = f"{cfg.name} mappings."
        if cfg.species:
            choices = _species_choices_text(cfg_id)
            if choices:
                group_help += f" Species: {choices}."
        grp = click.Group(cfg_id.replace("_", "-"), help=group_help)
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


# ensembl label-history (cross-release; not a single-release generate_* command)


@click.command("label-history")
@_opt_species_for("ensembl")
@click.option("--from-version", default=None, help="Lower bound (inclusive) on the release walk.")
@click.option("--to-version", default=None, help="Upper bound (inclusive) on the release walk.")
@_opt_output
@click.option("--cache-dir", type=PathType, default=None, help="Cache directory.")
@click.option(
    "--force", is_flag=True, default=False, help="Re-walk every release, ignoring resume state."
)
@_opt_no_progress
def ensembl_label_history_cmd(
    species: str,
    from_version: str | None,
    to_version: str | None,
    output: Path | None,
    cache_dir: Path | None,
    force: bool,
    no_progress: bool,
) -> None:
    """Derive cross-release previous-symbol -> current-symbol Ensembl mappings.

    Ensembl's core schema has no previous-gene-symbol table, so this walks
    every historical release (or a bounded --from-version/--to-version
    range) and diffs each release's primary-label snapshot to recover
    true previous -> current symbol transitions. Network-heavy and
    resumable; run on demand, not part of normal mapping generation.
    """
    from pysec2pri.api import generate_ensembl_label_history

    click.echo("Walking Ensembl releases for label history (this may take a while)...")
    ms = generate_ensembl_label_history(
        species=species,
        from_version=from_version,
        to_version=to_version,
        cache_dir=cache_dir,
        show_progress=not no_progress,
        force=force,
    )
    n = len(ms.mappings or [])
    click.echo(f"Found {n} previous -> current label transition(s).")
    if output:
        ms.save("sssom", output)
        click.echo(f"Wrote -> {output}")
    else:
        click.echo("(Also cached under the consolidate cache directory; see --cache-dir.)")


cast(click.Group, main.commands["ensembl"]).add_command(ensembl_label_history_cmd)


# diff


@main.command()
@click.argument("file1", type=ExistingPathType)
@click.argument("file2", type=ExistingPathType)
@click.option("-o", "--output", type=PathType, help="Output file for diff results (TSV).")
@click.option("--show-all", is_flag=True, default=False, help="Show all differences.")
@click.option("--datasource", default="unknown", help="Datasource name for diff summary.")
def diff(file1: Path, file2: Path, output: Path | None, show_all: bool, datasource: str) -> None:
    """Compare two SSSOM mapping files and show differences."""
    from mapkgsutils.diff import diff_sssom_files, summarize_diff

    from pysec2pri.api import write_diff_output

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
    type=click.Choice(["chebi", "ensembl", "hgnc", "uniprot"], case_sensitive=False),
)
def list_versions_cmd(datasource: str) -> None:
    """List available archive versions for DATASOURCE (chebi, ensembl, hgnc, uniprot)."""
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
        from pysec2pri.downloads import ChEBIDownloader

        threshold = ChEBIDownloader().new_format_version or 245
        for v in versions:
            click.echo(f"  {v}" + ("  (SDF only)" if int(v) < threshold else ""))
        click.echo(f"\n  Note: versions < {threshold} are SDF only.")
    else:
        for v in versions:
            click.echo(f"  {v}")


# validate-config


@main.command("validate-config")
@click.argument("datasource", required=False, default=None)
def validate_config_cmd(datasource: str | None) -> None:
    """Validate one (or, if omitted, every) shipped datasource config YAML.

    Examples::

        pysec2pri validate-config
        pysec2pri validate-config hgnc
    """
    from mapkgsutils.config.schema import ConfigValidationError, validate_config_file

    from pysec2pri.parsers.base import CONFIG_DIR

    if datasource:
        paths = [CONFIG_DIR / f"{datasource.lower()}.yaml"]
    else:
        paths = sorted(CONFIG_DIR.glob("*.yaml"))

    failed = False
    for path in paths:
        if not path.exists():
            click.echo(f"  {path.name}: NOT FOUND", err=True)
            failed = True
            continue
        try:
            validate_config_file(path)
        except ConfigValidationError as exc:
            click.echo(f"  {path.name}: INVALID -- {exc}", err=True)
            failed = True
        else:
            click.echo(f"  {path.name}: OK")

    if failed:
        raise click.ClickException("One or more config files failed validation.")


# consolidate (nested per-datasource: only datasources with a versioned
# archive/per-row date get this command, and --subset/--species are only
# attached when the datasource's own config declares that axis)


def _consolidate_extra_opts(cfg_id: str) -> list[Callable[..., Any]]:
    """Return the ``--subset``/``--species`` option(s) *cfg_id*'s config declares.

    Driven by ``DatasourceConfig.subset``/``.species`` (see
    ``config/<cfg_id>.yaml``) rather than a hardcoded per-datasource list, so
    e.g. HGNC never shows an irrelevant ``--subset`` flag.
    """
    cfg = get_datasource_config(cfg_id, config_package=_CONFIG_PACKAGE)
    opts: list[Callable[..., Any]] = []
    if cfg.subset:
        opts.append(_opt_subset)
    if cfg.species:
        opts.append(_opt_species_for(cfg_id))
    return opts


def _make_consolidate_cmd(cfg_id: str, extra_opts: list[Callable[..., Any]]) -> Callable[..., Any]:
    """Return a decorated (but not yet click.command-wrapped) ``consolidate`` callable."""

    def _cmd(
        mode: str,
        mapping_sets: str,
        cache_dir: Path | None,
        force: bool,
        no_progress: bool,
        **extra_kwargs: Any,
    ) -> None:
        from pysec2pri.consolidate import _sssom_output_path, consolidate_mapping_dates

        extras = ", ".join(f"{k}={v}" for k, v in extra_kwargs.items() if v is not None)
        click.echo(
            f"Consolidating {cfg_id.upper()} mapping dates "
            f"(mode={mode}{', ' + extras if extras else ''}, {mapping_sets})..."
        )
        try:
            path, _ = consolidate_mapping_dates(
                cfg_id,
                mode=mode,
                cache_dir=cache_dir,
                mapping_sets=mapping_sets,
                show_progress=not no_progress,
                force=force,
                **extra_kwargs,
            )
        except ValueError as exc:
            raise click.ClickException(str(exc)) from exc
        click.echo(f"Wrote consolidated mapping set -> {_sssom_output_path(path)}")

    decorators: list[Callable[..., Any]] = [
        click.option(
            "--mode",
            default="release",
            show_default=True,
            type=click.Choice(["release", "date"]),
            help=(
                "'release': walk every historical release for the first-seen version/date. "
                "'date': single pass using the source's own per-row date (falls back to "
                "'release' with a warning if the source has none)."
            ),
        ),
        *extra_opts,
        click.option(
            "--mapping-sets",
            default="ids",
            show_default=True,
            type=click.Choice(["ids", "labels"]),
            help="Which mapping-set kind to consolidate dates for.",
        ),
        click.option("--cache-dir", type=PathType, default=None, help="Cache directory."),
        click.option(
            "--force",
            is_flag=True,
            default=False,
            help="Re-scan every release, ignoring resume state.",
        ),
        _opt_no_progress,
    ]
    for dec in reversed(decorators):
        _cmd = dec(_cmd)

    _cmd.__name__ = "consolidate"
    return _cmd


def _register_consolidate_commands(parent: click.Group) -> None:
    """Register a ``consolidate`` subcommand on every supported datasource's group."""
    from pysec2pri.consolidate import SUPPORTED_DATASOURCES

    for cfg_id in SUPPORTED_DATASOURCES:
        group = parent.commands.get(cfg_id)
        if group is None:
            continue
        raw_cmd = _make_consolidate_cmd(cfg_id, _consolidate_extra_opts(cfg_id))
        cast(click.Group, group).add_command(
            click.command(
                name="consolidate",
                help=(
                    "Build/update the per-mapping first-seen-date index, written as a real "
                    "SSSOM mapping set whose mapping_date is each mapping's date of first "
                    "appearance.\n\n"
                    "In 'release' mode (default), walks every historical release once to "
                    "discover the earliest release each mapping appeared in. This is slow "
                    "and network-heavy (~250 releases for ChEBI); meant to be run as a "
                    "manual/one-off operation, not as part of normal mapping generation. "
                    "In 'date' mode, a single current-release pass captures the source's "
                    "own per-row date directly (fast)."
                ),
            )(raw_cmd)
        )


_register_consolidate_commands(main)


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
    "vgnc": "generate_vgnc_labels",
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
@_opt_xref
@_opt_xref_file
@_opt_xref_source
@_opt_xref_on
@_opt_xref_predicate
@_opt_report
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
    xref_cols: tuple[str, ...],
    xref_file: Path | None,
    xref_source: str | None,
    xref_on: str | None,
    xref_predicates: tuple[str, ...],
    report_path: Path | None,
    data_version: str | None,
    no_progress: bool,
) -> None:
    r"""Resolve secondary IDs in INPUT_FILE to primary IDs using DATASOURCE.

    Examples::

        pysec2pri update-ids genes.tsv hgnc --at gene_id -o out.tsv
        pysec2pri update-ids genes.tsv hgnc --at gene_id --synonyms label
        pysec2pri update-ids genes.tsv hgnc --at gene_id --xref ensembl \\
            --xref-source hgnc_custom --xref-on ensembl --report decisions.tsv
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

    xref_mapping = _resolve_xref_mapping(xref_file, xref_source, xref_on, datasource)
    if xref_cols and xref_mapping is None:
        raise click.ClickException("--xref requires --xref-file or --xref-source.")

    syn_padded = _pad_cols(synonyms_cols, len(columns))
    xref_padded = _pad_cols(xref_cols, len(columns))
    _resolve_and_print(
        input_file,
        ms,
        list(zip(columns, syn_padded, xref_padded, strict=True)),
        output_path,
        suffix,
        sep,
        label_ms,
        mode="ids",
        xref_mapping=xref_mapping,
        xref_predicates=set(xref_predicates) if xref_predicates else None,
        report_path=report_path,
    )


# update-labels

_LABELS_DATASOURCES = sorted(cfg for cfg, kind in _build_registry() if kind == "labels")
_ID_GENERATORS_FOR_LABELS: dict[str, str] = {
    "hgnc": "generate_hgnc_ids",
    "ncbi": "generate_ncbi_ids",
    "vgnc": "generate_vgnc_ids",
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
@_opt_xref
@_opt_xref_file
@_opt_xref_source
@_opt_xref_on
@_opt_xref_predicate
@_opt_report
@click.option(
    "--species",
    default=None,
    help=(
        "Species as NCBI taxon ID. Defaults to DATASOURCE's own config default when omitted; "
        "run 'pysec2pri DATASOURCE labels --help' to see its known species."
    ),
)
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
    xref_cols: tuple[str, ...],
    xref_file: Path | None,
    xref_source: str | None,
    xref_on: str | None,
    xref_predicates: tuple[str, ...],
    report_path: Path | None,
    species: str | None,
    entity_type: str | None,
    subset: str,
    data_version: str | None,
    no_progress: bool,
) -> None:
    r"""Resolve previous/alias labels in INPUT_FILE to current labels using DATASOURCE.

    Examples::

        pysec2pri update-labels genes.tsv hgnc --at label -o out.tsv
        pysec2pri update-labels genes.tsv hgnc --at label --mapping labels.tsv
        pysec2pri update-labels genes.tsv hgnc --at label --xref ensembl \\
            --xref-source hgnc_custom --xref-on ensembl --report decisions.tsv
    """
    from pysec2pri.api import (
        generate_chebi_labels,
        generate_hgnc_labels,
        generate_ncbi_labels,
        generate_vgnc_labels,
        generate_wikidata_labels,
        load_label_mapping,
        load_mapping,
    )

    if species is None:
        ds_cfg = get_datasource_config(datasource, config_package=_CONFIG_PACKAGE)
        species = str(ds_cfg.default_species()) if ds_cfg.species else "9606"

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
                species=species, version=data_version, show_progress=not no_progress
            ),
            "vgnc": lambda: generate_vgnc_labels(
                species=species, version=data_version, show_progress=not no_progress
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

    xref_mapping = _resolve_xref_mapping(xref_file, xref_source, xref_on, datasource)
    if xref_cols and xref_mapping is None:
        raise click.ClickException("--xref requires --xref-file or --xref-source.")

    syn_padded = _pad_cols(synonyms_cols, len(columns))
    xref_padded = _pad_cols(xref_cols, len(columns))
    _resolve_and_print(
        input_file,
        ms,
        list(zip(columns, syn_padded, xref_padded, strict=True)),
        output_path,
        suffix,
        sep,
        label_ms,
        mode="labels",
        xref_mapping=xref_mapping,
        xref_predicates=set(xref_predicates) if xref_predicates else None,
        report_path=report_path,
    )


# crosswalk


def _crosswalk_vocab() -> list[str]:
    """``--from`` vocabulary: ``symbol`` plus every key declared in HGNC's xref-source."""
    cfg = get_datasource_config("hgnc", config_package=_CONFIG_PACKAGE)
    src = cfg.xref_source("hgnc_custom")
    keys = sorted(src.subject_id_cols) if src else []
    return ["symbol", *keys]


@main.command("crosswalk")
@click.argument("input_file", type=ExistingPathType, required=False, default=None)
@click.option(
    "--value", default=None, help="A single identifier to look up (instead of INPUT_FILE)."
)
@click.option(
    "--from",
    "frm",
    required=True,
    type=click.Choice(_crosswalk_vocab()),
    help="Source vocabulary.",
)
@click.option(
    "--to",
    default="hgnc_id",
    show_default=True,
    type=click.Choice(["hgnc_id", "symbol"]),
    help="Target vocabulary.",
)
@click.option(
    "--at", "at_col", default=None, metavar="COLUMN", help="Column with --from values (file mode)."
)
@_opt_output
@click.option("--sep", default=None, help="Delimiter (inferred from extension if omitted).")
@_opt_xref_file
@_opt_xref_source
@_opt_report
@_opt_version
@_opt_no_progress
def crosswalk_cmd(
    input_file: Path | None,
    value: str | None,
    frm: str,
    to: str,
    at_col: str | None,
    output: Path | None,
    sep: str | None,
    xref_file: Path | None,
    xref_source: str | None,
    report_path: Path | None,
    data_version: str | None,
    no_progress: bool,
) -> None:
    """Map a gene identifier between vocabularies via HGNC.

    Vocabularies: symbol, hgnc_id, and any cross-reference key HGNC's
    custom-download crosswalk declares (typically ensembl, entrez, refseq,
    uniprot). A ``symbol`` lookup resolves through the temporal-aware,
    ambiguity-safe label resolver; others resolve through HGNC's own
    cross-reference table.

    Examples::

        pysec2pri crosswalk --value TP53 --from symbol --to hgnc_id
        pysec2pri crosswalk genes.tsv --from ensembl --to hgnc_id --at ensembl_id -o out.tsv
    """
    from mapkgsutils.context import load_xref_mapping

    from pysec2pri.api import crosswalk as _crosswalk

    xref_mapping = load_xref_mapping(xref_file) if xref_file is not None else None

    if value is not None:
        looked_up = _crosswalk(
            value,
            frm=frm,
            to=to,
            version=data_version,
            xref_mapping=xref_mapping,
            xref_source=xref_source or "hgnc_custom",
            report_path=report_path,
            show_progress=not no_progress,
        )
        click.echo(looked_up.get(value, "") if isinstance(looked_up, dict) else looked_up)
        return

    if input_file is None:
        raise click.ClickException("Provide INPUT_FILE or --value.")
    if at_col is None:
        raise click.ClickException("--at is required with INPUT_FILE.")

    import pandas as pd

    read_sep = sep if sep is not None else ("\t" if input_file.suffix.lower() == ".tsv" else ",")
    df = pd.read_csv(input_file, sep=read_sep, dtype=str)
    result_df = _crosswalk(
        df,
        frm=frm,
        to=to,
        at=at_col,
        version=data_version,
        xref_mapping=xref_mapping,
        xref_source=xref_source or "hgnc_custom",
        report_path=report_path,
        show_progress=not no_progress,
    )
    assert isinstance(result_df, pd.DataFrame)  # guaranteed: input_data was a DataFrame
    if output is not None:
        out_sep = "\t" if output.suffix.lower() == ".tsv" else ","
        result_df.to_csv(output, sep=out_sep, index=False)
        click.echo(f"Wrote {len(result_df)} rows -> {output}")
    else:
        out_sep = "\t" if (sep or "\t") == "\t" else ","
        click.echo(result_df.to_csv(sep=out_sep, index=False))


if __name__ == "__main__":
    main()
