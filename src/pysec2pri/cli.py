"""Command-line interface for pysec2pri.

Nothing here names a datasource. Every command is registered from the config
files: each source's ``mapping_sets`` gives its subcommands, and each command's
options come from the same config -- one ``--<input>`` per declared input, one
``--<product>`` per name in ``products``, ``--format`` from that kind's
``formats``, and ``--consolidate`` where the kind declares it. Adding a source
is a config plus parser change; the CLI needs no edit.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable
from importlib import resources as _importlib_resources
from pathlib import Path
from typing import TYPE_CHECKING, Any

import click
from mapkgsutils.parsers.config import get_datasource_config, product_dimensions

if TYPE_CHECKING:
    from pysec2pri.parsers.base import BaseMappingSet

warnings.filterwarnings("ignore", category=FutureWarning, module="sssom")

_CONFIG_PACKAGE = "pysec2pri.config"

#: Config directory, resolved without importing the parser stack: command
#: registration runs at import and must stay free of sssom/rdflib/polars.
CONFIG_DIR = Path(_importlib_resources.files(_CONFIG_PACKAGE))  # type: ignore[arg-type]

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


def _opt_subset_for(cfg_id: str) -> Callable[..., Any]:
    """Build a ``--subset`` option from ``<cfg_id>.yaml``'s ``subset`` block.

    Choices come from ``subset.available`` and the default from
    ``subset.default``, so a datasource's entry subsets are declared once.
    """
    cfg = get_datasource_config(cfg_id, config_package=_CONFIG_PACKAGE)
    available = cfg.subset.get("available") or {}
    help_text = "Entry subset."
    labels = [
        f"{key} ({entry['label']})"
        for key, entry in available.items()
        if isinstance(entry, dict) and entry.get("label")
    ]
    if labels:
        help_text += f" {'; '.join(labels)}."
    return _opt(
        "--subset",
        default=cfg.default_subset(),
        show_default=True,
        type=click.Choice(sorted(available)),
        help=help_text,
    )


_MAX_INLINE_SPECIES = 40


def _species_choices_text(cfg_id: str, limit: int = _MAX_INLINE_SPECIES) -> str:
    """Return a ``"<taxon_id>=<label>, ..."`` list from the config's ``species.available``."""
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
    """Build a ``--species`` option defaulting to ``<cfg_id>.yaml``'s ``species.default``."""
    cfg = get_datasource_config(cfg_id, config_package=_CONFIG_PACKAGE)
    help_text = "Species as NCBI taxon ID, or 'all' to process every species."
    choices = _species_choices_text(cfg_id)
    if choices:
        help_text += f" Known: {choices}."
    note = (cfg.species or {}).get("note")
    if note:
        help_text += f" {' '.join(str(note).split())}"
    return _opt("--species", default=str(cfg.default_species()), show_default=True, help=help_text)


#: ``--subset``
_opt_subset_any = _opt(
    "--subset",
    default=None,
    help="Entry subset. Defaults to DATASOURCE's configured subset when omitted.",
)


def _opt_entity_type_for(cfg_id: str) -> Callable[..., Any]:
    """Build an ``--entity-type`` option from ``<cfg_id>.yaml``'s ``queries`` keys."""
    cfg = get_datasource_config(cfg_id, config_package=_CONFIG_PACKAGE)
    return _opt(
        "--entity-type",
        default=None,
        type=click.Choice(list(cfg.queries)),
        help=f"{cfg.name} entity type to query. Queries all if omitted.",
    )


_opt_test_subset = _opt(
    "--test-subset",
    is_flag=True,
    default=False,
    help="Use test SPARQL queries (LIMIT 10).",
)
# Hints for update-ids/update-labels. A value is ambiguous when it is both
# retired and current; these say which one a row means. All optional: without
# them an ambiguous row gets an empty cell.
_opt_xref = _opt(
    "--xref",
    "xref_cols",
    default=None,
    multiple=True,
    metavar="COLUMN",
    help="Column holding an identifier from another vocabulary for the same row, "
    "as a hint. Pairs with --at, and needs --xref-file or --xref-source.",
)
_opt_xref_file = _opt(
    "--xref-file",
    "xref_file",
    type=ExistingPathType,
    default=None,
    help="Table saying which of this source's entries each --xref value belongs to "
    "(SSSOM or plain TSV).",
)
_opt_xref_source = _opt(
    "--xref-source",
    "xref_source",
    default=None,
    help="Same as --xref-file, but downloaded: name one of the crosswalks DATASOURCE's "
    "config lists, e.g. 'hgnc_custom'.",
)
_opt_xref_on = _opt(
    "--xref-on",
    "xref_on",
    default=None,
    help="Which vocabulary the --xref values are from, e.g. ensembl/entrez/refseq/uniprot. "
    "Required with --xref-source.",
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

    from pysec2pri.update import update_ids, update_labels

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


def _opt_inputs_for(config_id: str, kind: str) -> list[Callable[..., Any]]:
    """Build one ``--<input>`` option per input the kind declares.

    Each of a mapping set's inputs is a real file, so each gets its own named
    option (``--withdrawn``, ``--complete`` ...) instead of one overloaded
    ``--input-file``. Anything not passed is downloaded.
    """
    cfg = get_datasource_config(config_id, config_package=_CONFIG_PACKAGE)
    spec = cfg.mapping_sets.get(kind) or {}
    return [
        _opt(
            f"--{key.replace('_', '-')}",
            key,
            type=ExistingPathType,
            default=None,
            help=f"Local {key} file. Downloaded when omitted.",
        )
        for key in (spec.get("inputs") or {})
    ]


_opt_consolidate = _opt(
    "--consolidate",
    is_flag=True,
    default=False,
    help=(
        "Also recover mappings the current release no longer states, by walking the "
        "source's past releases, and stamp each mapping with the release it first "
        "appeared in. Slow and network-heavy; resumable."
    ),
)
_opt_cache_dir = _opt(
    "--cache-dir",
    type=PathType,
    default=None,
    help="Where to keep the resumable --consolidate index. Defaults to this OS's "
    "per-user cache directory, or $PYSEC2PRI_CACHE_DIR when set.",
)
_opt_force = _opt(
    "--force",
    is_flag=True,
    default=False,
    help="With --consolidate, re-walk every release, ignoring resume state.",
)
_opt_from_version = _opt(
    "--from-version",
    default=None,
    help="With --consolidate, lower bound (inclusive) on the release walk.",
)
_opt_to_version = _opt(
    "--to-version",
    default=None,
    help="With --consolidate, upper bound (inclusive) on the release walk.",
)


def _make_generate_cmd(
    config_id: str,
    kind: str,
    extra_opts: list[Callable[..., Any]],
) -> Callable[..., Any]:
    """Return a decorated (but not yet click.command-wrapped) callable.

    Everything the command needs comes from *config_id*'s config: the inputs it
    accepts, the formats it can emit, and whether it supports --consolidate.
    """
    cfg = get_datasource_config(config_id, config_package=_CONFIG_PACKAGE)
    spec = cfg.mapping_sets.get(kind) or {}
    formats = cfg.formats_for(kind)
    input_keys = list(spec.get("inputs") or {})
    can_consolidate = bool(spec.get("consolidate"))

    def _cmd(
        output: Path | None,
        data_version: str | None,
        output_format: str,
        no_progress: bool,
        **extra_kwargs: Any,
    ) -> None:
        from pysec2pri import api

        inputs = {key: extra_kwargs.pop(key) for key in input_keys if extra_kwargs.get(key)}
        for key in input_keys:
            extra_kwargs.pop(key, None)
        ms = api._generate(
            config_id,
            kind,
            version=data_version,
            show_progress=not no_progress,
            inputs=inputs or None,
            **extra_kwargs,
        )
        prefix = f"{config_id}_{kind}"
        species = extra_kwargs.get("species")
        if species is not None:
            prefix = f"{prefix}_{species}"
        _, base = _version_base(ms, data_version, prefix)
        if extra_kwargs.get("consolidate"):
            # A consolidated walk is a distinct data product at the same
            # release, so it must not overwrite the plain run's output.
            base = f"{base}_consolidate"
        _emit(ms, output_format, output, base)

    # Apply Click decorators
    decorators: list[Callable[..., Any]] = [
        *_opt_inputs_for(config_id, kind),
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
    if can_consolidate:
        decorators += [
            _opt_consolidate,
            _opt_cache_dir,
            _opt_force,
            _opt_from_version,
            _opt_to_version,
        ]
    for dec in reversed(decorators):
        _cmd = dec(_cmd)

    _cmd.__name__ = kind
    return _cmd


# Datasource registry


#: Mapping-set kinds the CLI exposes. Primary ID/label listings are an API
#: concern, outside this package's command-line focus.
_CLI_KINDS = ("ids", "labels")


def _archive_sources() -> list[str]:
    """Return the sources whose config declares an archive of past releases."""
    return [
        path.stem
        for path in sorted(CONFIG_DIR.glob("*.yaml"))
        if get_datasource_config(path.stem, config_package=_CONFIG_PACKAGE).archive_url
    ]


_ARCHIVE_SOURCES = _archive_sources()


def _opt_product_generic(name: str) -> Callable[..., Any]:
    """Build a plain ``--<name>`` option for a dimension with no richer factory."""
    return _opt(
        f"--{name.replace('_', '-')}",
        default=None,
        help=f"Which {name.replace('_', ' ')} to build the mapping set for.",
    )


#: Dimension name -> the option factory that reads its choices and default out
#: of the config block naming them. A dimension with no entry here still gets
#: an option, just without choices; see :func:`_opt_product_generic`.
_PRODUCT_OPTS: dict[str, Callable[[str], Callable[..., Any]]] = {
    "subset": _opt_subset_for,
    "species": _opt_species_for,
    "entity_type": _opt_entity_type_for,
}


def _extra_opts_for(config_id: str) -> list[Callable[..., Any]]:
    """Return the option decorators *config_id*'s config calls for.

    A source gets one option per name in its ``products``, plus
    ``--test-subset`` when it queries a SPARQL endpoint. So a new source needs
    no CLI change: declare the dimension, get the option.
    """
    cfg = get_datasource_config(config_id, config_package=_CONFIG_PACKAGE)
    opts: list[Callable[..., Any]] = []
    for name in product_dimensions(cfg):
        factory = _PRODUCT_OPTS.get(name)
        opts.append(factory(config_id) if factory else _opt_product_generic(name))
    if cfg.queries:
        opts.append(_opt_test_subset)
    return opts


def _build_registry() -> dict[tuple[str, str], list[Callable[..., Any]]]:
    """Return (config_id, kind) -> [extra option decorators], read from the configs.

    Both the datasources and the mapping-set kinds they publish come from the
    config files, so adding a source is a config + parser change only.
    """
    registry: dict[tuple[str, str], list[Callable[..., Any]]] = {}
    for path in sorted(CONFIG_DIR.glob("*.yaml")):
        cfg_id = path.stem
        cfg = get_datasource_config(cfg_id, config_package=_CONFIG_PACKAGE)
        for kind in cfg.mapping_sets:
            if kind in _CLI_KINDS:
                registry[(cfg_id, kind)] = _extra_opts_for(cfg_id)
    return registry


def _has_kind(config_id: str, kind: str) -> bool:
    """Whether *config_id* publishes a *kind* mapping set."""
    cfg = get_datasource_config(config_id, config_package=_CONFIG_PACKAGE)
    return kind in cfg.mapping_sets


def _load_mapping_set(config_id: str, kind: str, **options: Any) -> Any:
    """Build *config_id*'s *kind* mapping set, passing only options it accepts."""
    from pysec2pri import api

    return api._generate(config_id, kind, **{k: v for k, v in options.items() if v is not None})


def _register_datasources(parent: click.Group) -> None:
    """Register one Click group per config_id on *parent*."""
    by_config: dict[str, dict[str, list[Callable[..., Any]]]] = {}
    for (cfg_id, kind), opts in _build_registry().items():
        by_config.setdefault(cfg_id, {})[kind] = opts

    for cfg_id, kinds in by_config.items():
        cfg = get_datasource_config(cfg_id, config_package=_CONFIG_PACKAGE)
        group_help = f"{cfg.name} mappings."
        if cfg.species:
            choices = _species_choices_text(cfg_id)
            if choices:
                group_help += f" Species: {choices}."
        grp = click.Group(cfg_id.replace("_", "-"), help=group_help)
        for kind, extra_opts in kinds.items():
            raw_cmd = _make_generate_cmd(cfg_id, kind, extra_opts)
            grp.add_command(click.command(name=kind)(raw_cmd))
        parent.add_command(grp)


# Main entry point


@click.group()
@click.version_option()
def main() -> None:
    """pysec2pri -- secondary-to-primary ID and label mapping.

    Each source command returns `ids` and, where available, `labels`: the retired
    identifiers and superseded or alias labels that its current release states,
    as an SSSOM mapping set.

    Sources vary in how much of their own history they keep, and some drop
    retired entries entirely. `--consolidate` goes over all available releases
    to recover the mappings the current release no longer states, and stamp every
    mapping with the release (date) it first appeared in.
    """


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
@click.argument("datasource", type=click.Choice(_ARCHIVE_SOURCES, case_sensitive=False))
def list_versions_cmd(datasource: str) -> None:
    """List the releases DATASOURCE keeps in its archive.

    Only sources whose config declares an ``archive_url`` have one; the rest
    publish their current release only.
    """
    from pysec2pri.api import list_versions as _list_versions

    click.echo(f"Fetching available {datasource.upper()} versions...")
    try:
        versions = _list_versions(datasource)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    if not versions:
        click.echo("No versions found.")
        return

    # A source may have distributed different formats over its lifetime; its
    # config's distribution_eras say which applies to each release.
    cfg = get_datasource_config(datasource.lower(), config_package=_CONFIG_PACKAGE)
    click.echo(f"\n{len(versions)} version(s) for {datasource.upper()}:\n")
    for v in versions:
        era = cfg.era_for(v)
        fmt = era.format if era is not None else None
        click.echo(f"  {v}" + (f"  ({fmt})" if fmt else ""))


# validate-config


@main.command("validate-config")
@click.argument("datasource", required=False, default=None)
def validate_config_cmd(datasource: str | None) -> None:
    """Validate one (or, if omitted, every) datasource config YAML.

    Examples::

        pysec2pri validate-config
        pysec2pri validate-config hgnc
    """
    from mapkgsutils.config.schema import ConfigValidationError, validate_config_file

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

    click.echo(f"Loading {datasource.upper()} mappings...")
    from pysec2pri import api

    cfg_id, kind = match
    ms = api._generate(cfg_id, kind, version=data_version, show_progress=not no_progress)

    from pysec2pri.api import DEFAULT_OUTPUT_DIR, find_ambiguous

    click.echo("Detecting ambiguous identifiers...")
    amb = find_ambiguous(ms)
    count = len(amb.mappings or [])
    if count == 0:
        click.echo("No ambiguous identifiers found.")
        return
    click.echo(f"Found {count} ambiguous identifier(s).")
    dest = output or (DEFAULT_OUTPUT_DIR / f"{datasource}_ambiguous.sssom.tsv")
    amb.save("sssom", dest)
    click.echo(f"Wrote ambiguous mappings -> {dest}")


# update-ids

_ID_DATASOURCES = sorted(cfg for cfg, kind in _build_registry() if kind == "ids")


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

    Reads the --at column and writes a new column with each ID's current one.

    An ID that is both retired and still current is ambiguous: there is no
    safe answer, so the new column is left empty for that row. To resolve
    those, give a hint: another column of the same row saying which entry it
    means. All hints are optional.

    --synonyms names a column of names for the row. --xref names a column of
    identifiers from another vocabulary, and needs a table to read them
    against: --xref-file for your own, or --xref-source to download one
    DATASOURCE's config lists. --report writes down every decision.

    Examples::

        pysec2pri update-ids genes.tsv hgnc --at gene_id -o out.tsv
        pysec2pri update-ids genes.tsv hgnc --at gene_id --synonyms label
        pysec2pri update-ids genes.tsv hgnc --at gene_id --xref ensembl \\
            --xref-source hgnc_custom --xref-on ensembl --report decisions.tsv
    """
    from pysec2pri.api import load_label_mapping, load_mapping

    if mapping_file is not None:
        click.echo(f"Loading mappings from {mapping_file}...")
        ms = load_mapping(mapping_file)
    else:
        click.echo(f"Loading {datasource.upper()} mappings...")
        ms = _load_mapping_set(
            datasource, "ids", version=data_version, show_progress=not no_progress
        )

    label_ms = None
    if synonyms_cols:
        if synonyms_mapping_file is not None:
            click.echo(f"Loading synonym mappings from {synonyms_mapping_file}...")
            try:
                label_ms = load_label_mapping(synonyms_mapping_file)
            except Exception:
                label_ms = load_mapping(synonyms_mapping_file)
        elif _has_kind(datasource, "labels"):
            click.echo(f"Loading {datasource.upper()} label mappings for alias resolution...")
            label_ms = _load_mapping_set(
                datasource, "labels", version=data_version, show_progress=not no_progress
            )
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
        "Species as NCBI taxon ID. Defaults to DATASOURCE config default when omitted; "
        "run 'pysec2pri DATASOURCE labels --help' to see its known species."
    ),
)
@_opt_entity_type_for("wikidata")
@_opt_subset_any
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

    Reads the --at column and writes a new column with each label's current one.

    A label that is both old and still in use for something else is
    ambiguous: there is no safe answer, so the new column is left empty for
    that row. To resolve those, give a hint: another column of the same row
    saying which entry it means. All hints are optional.

    --synonyms names a column of names for the row. --xref names a column of
    identifiers from another vocabulary, and needs a table to read them
    against: --xref-file for your own, or --xref-source to download one
    DATASOURCE's config lists. --report writes down every decision.

    Examples::

        pysec2pri update-labels genes.tsv hgnc --at label -o out.tsv
        pysec2pri update-labels genes.tsv hgnc --at label --mapping labels.tsv
        pysec2pri update-labels genes.tsv hgnc --at label --xref ensembl \\
            --xref-source hgnc_custom --xref-on ensembl --report decisions.tsv
    """
    from pysec2pri.api import load_label_mapping, load_mapping

    if species is None:
        ds_cfg = get_datasource_config(datasource, config_package=_CONFIG_PACKAGE)
        species = str(ds_cfg.default_species()) if ds_cfg.species else "9606"

    if mapping_file is not None:
        click.echo(f"Loading label mappings from {mapping_file}...")
        ms = load_label_mapping(mapping_file)
    else:
        click.echo(f"Loading {datasource.upper()} label mappings...")
        if not _has_kind(datasource, "labels"):
            raise click.ClickException(f"{datasource!r} has no label mapping set.")
        ms = _load_mapping_set(
            datasource,
            "labels",
            version=data_version,
            show_progress=not no_progress,
            species=species,
            subset=subset,
            entity_type=entity_type,
        )

    label_ms = None
    if synonyms_cols:
        if synonyms_mapping_file is not None:
            click.echo(f"Loading synonym mappings from {synonyms_mapping_file}...")
            try:
                label_ms = load_label_mapping(synonyms_mapping_file)
            except Exception:
                label_ms = load_mapping(synonyms_mapping_file)
        elif _has_kind(datasource, "ids"):
            click.echo(f"Loading {datasource.upper()} ID mappings for alias resolution...")
            label_ms = _load_mapping_set(
                datasource, "ids", version=data_version, show_progress=not no_progress
            )
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
