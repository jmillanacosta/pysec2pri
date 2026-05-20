"""Command-line interface for pysec2pri."""

from __future__ import annotations

import warnings
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, TypeVar

import click

if TYPE_CHECKING:
    import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning, module="sssom")

PathType = click.Path(path_type=Path)  # type: ignore[type-var]
ExistingPathType = click.Path(exists=True, path_type=Path)  # type: ignore[type-var]

_F = TypeVar("_F", bound=Callable[..., object])


# Reusable options

_ID_FORMATS = ["sssom", "sec2pri", "pri_ids", "rdf", "json", "owl", "all"]
_LABEL_FORMATS = ["sssom", "symbol_sec2pri", "pri_symbols", "rdf", "json", "owl", "all"]
_SYNONYM_FORMATS = [
    "sssom",
    "symbol_sec2pri",
    "name2synonym",
    "pri_symbols",
    "rdf",
    "json",
    "owl",
    "all",
]


def _opt_output(fn: _F) -> _F:
    return click.option("-o", "--output", type=PathType, help="Output file or directory")(fn)


def _opt_version(fn: _F) -> _F:
    return click.option(
        "--version", "data_version", default=None, help="Datasource release version."
    )(fn)


def _opt_id_format(fn: _F) -> _F:
    return click.option(
        "--format",
        "output_format",
        default="sssom",
        show_default=True,
        type=click.Choice(_ID_FORMATS),
    )(fn)


def _opt_label_format(fn: _F) -> _F:
    return click.option(
        "--format",
        "output_format",
        default="sssom",
        show_default=True,
        type=click.Choice(_LABEL_FORMATS),
    )(fn)


def _opt_synonym_format(fn: _F) -> _F:
    return click.option(
        "--format",
        "output_format",
        default="sssom",
        show_default=True,
        type=click.Choice(_SYNONYM_FORMATS),
    )(fn)


def _opt_subset(fn: _F) -> _F:
    return click.option(
        "--subset",
        default="3star",
        show_default=True,
        type=click.Choice(["3star", "complete"]),
        help="Compound subset.",
    )(fn)


def _opt_tax_id(fn: _F) -> _F:
    return click.option("--tax-id", default="9606", show_default=True, help="Taxonomy ID.")(fn)


def _opt_entity_type(fn: _F) -> _F:
    return click.option(
        "--entity-type",
        default=None,
        type=click.Choice(["metabolites", "chemicals", "genes", "proteins"]),
        help="Entity type to query. Queries all if omitted.",
    )(fn)


def _opt_no_progress(fn: _F) -> _F:
    return click.option(
        "--no-progress", is_flag=True, default=False, help="Suppress progress bars."
    )(fn)


# Shared helpers


def _version_base(ms: object, data_version: str | None, prefix: str) -> tuple[str, str]:
    """Return ``(version, base_name)`` for output filenames."""
    version = getattr(ms, "mapping_set_version", None) or data_version or ""
    base = f"{prefix}{'_' + version if version else ''}"
    return version, base


def _emit(ms: object, fmt: str, output: Path | None, base: str) -> None:
    """Save *ms* and echo a count summary."""
    from pysec2pri.api import save

    path = save(ms, fmt, output, base_name=base)  # type: ignore[arg-type]
    if fmt == "pri_ids":
        count = len(getattr(ms, "_primary_ids", set()) or set())
        click.echo(f"Wrote {count} primary IDs to {path}")
    elif fmt == "pri_symbols":
        count = len(getattr(ms, "_primary_symbols", set()) or set())
        click.echo(f"Wrote {count} primary symbols to {path}")
    else:
        click.echo(f"Wrote {len(getattr(ms, 'mappings', None) or [])} mappings to {path}")


def _resolve_and_print(
    input_file: Path,
    ms: object,
    resolve_fn: Callable[..., pd.DataFrame | str | list[str]],
    columns: tuple[str, ...],
    output_path: Path | None,
    suffix: str,
    sep: str | None,
) -> None:
    """Compute ``at``, call *resolve_fn*, then print or write results."""
    at: str | list[str] = list(columns) if len(columns) > 1 else columns[0]
    click.echo(f"Resolving in column(s): {', '.join(columns)}")
    df = resolve_fn(input_file, ms, at=at, output_path=output_path, suffix=suffix, sep=sep)
    assert not isinstance(df, (str, list)), "Expected DataFrame in CLI mode"
    if output_path:
        click.echo(f"Wrote {len(df)} rows to {output_path}")
    else:
        click.echo(df.to_csv(sep="\t" if (sep or "\t") == "\t" else ",", index=False))


@click.group()
@click.version_option()
def main() -> None:
    """pysec2pri: Secondary to Primary ID mapping tool."""


@main.command()
@click.argument("file1", type=ExistingPathType)
@click.argument("file2", type=ExistingPathType)
@click.option("-o", "--output", type=PathType, help="Output file for diff results (TSV)")
@click.option("--show-all", is_flag=True, default=False, help="Show all differences")
@click.option("--datasource", default="unknown", help="Datasource name for diff summary")
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
            f"\n(Showing summary only - {result.total_changes} changes. Use --show-all to see all.)"
        )

    if output:
        if result.total_changes > 0:
            write_diff_output(result, output)
            click.echo(f"\nWrote diff to {output}")
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


@main.group()
def chebi() -> None:
    """Parse ChEBI data and generate mappings."""


@chebi.command("ids")
@click.argument("input_file", type=ExistingPathType, required=False)
@_opt_output
@_opt_version
@_opt_subset
@_opt_id_format
def chebi_ids(
    input_file: Path | None,
    output: Path | None,
    data_version: str | None,
    subset: str,
    output_format: str,
) -> None:
    """Generate ChEBI ID mappings (secondary  to primary ChEBI IDs).

    Formats: sec2pri (dict), pri_ids (set of current IDs), sssom/rdf/json/owl/all.
    """
    from pysec2pri.api import generate_chebi

    ms = generate_chebi(
        input_path=input_file, version=data_version, subset=subset, mapping_sets="ids"
    )
    subset_sfx = "_3star" if subset == "3star" else ""
    _, base = _version_base(ms, data_version, f"chebi{subset_sfx}")
    _emit(ms, output_format, output, base)


@chebi.command("synonyms")
@click.argument("input_file", type=ExistingPathType, required=False)
@_opt_output
@_opt_version
@_opt_subset
@_opt_synonym_format
def chebi_synonyms(
    input_file: Path | None,
    output: Path | None,
    data_version: str | None,
    subset: str,
    output_format: str,
) -> None:
    """Generate ChEBI synonym mappings (previous/alias name  to current name).

    Formats: symbol_sec2pri (dict), name2synonym, pri_symbols (set of current names),
    sssom/rdf/json/owl/all.
    """
    from pysec2pri.api import generate_chebi_synonyms

    ms = generate_chebi_synonyms(input_path=input_file, version=data_version, subset=subset)
    subset_sfx = "_3star" if subset == "3star" else ""
    _, base = _version_base(ms, data_version, f"chebi{subset_sfx}_synonyms")
    _emit(ms, output_format, output, base)


@main.command()
@click.option("--metabolites-file", type=ExistingPathType, default=None)
@click.option("--proteins-file", type=ExistingPathType, default=None)
@click.option("--metabolites-only", is_flag=True, default=False)
@click.option("--proteins-only", is_flag=True, default=False)
@_opt_output
@_opt_version
@_opt_id_format
def hmdb(
    metabolites_file: Path | None,
    proteins_file: Path | None,
    metabolites_only: bool,
    proteins_only: bool,
    output: Path | None,
    data_version: str | None,
    output_format: str,
) -> None:
    """Parse HMDB XML files and generate secondary-to-primary mappings."""
    from pysec2pri.api import generate_hmdb, generate_hmdb_proteins

    if metabolites_only and proteins_only:
        raise click.UsageError("--metabolites-only and --proteins-only are mutually exclusive.")

    want_metabolites = not proteins_only
    want_proteins = not metabolites_only
    if proteins_file is not None and metabolites_file is None and not metabolites_only:
        want_metabolites = False
    if metabolites_file is not None and proteins_file is None and not proteins_only:
        want_proteins = False

    if want_metabolites:
        click.echo("Parsing HMDB metabolites...")
        ms = generate_hmdb(metabolites_file, version=data_version)
        _, base = _version_base(ms, data_version, "hmdb_metabolites")
        _emit(ms, output_format, output, base)

    if want_proteins:
        click.echo("Parsing HMDB proteins...")
        ms = generate_hmdb_proteins(proteins_file, version=data_version)
        _, base = _version_base(ms, data_version, "hmdb_proteins")
        _emit(ms, output_format, output, base)


@main.group()
def hgnc() -> None:
    """Parse HGNC files and generate mappings."""


@hgnc.command("ids")
@click.argument("input_file", type=ExistingPathType, required=False)
@_opt_output
@_opt_version
@_opt_id_format
def hgnc_ids(
    input_file: Path | None,
    output: Path | None,
    data_version: str | None,
    output_format: str,
) -> None:
    """Generate HGNC ID mappings (secondary to primary HGNC IDs).

    Formats: sec2pri (dict), pri_ids (set of all current IDs), sssom/rdf/json/owl/all.
    """
    from pysec2pri.api import generate_hgnc, generate_hgnc_primary_ids

    if output_format == "pri_ids":
        ms = generate_hgnc_primary_ids(input_file, version=data_version)
    else:
        ms = generate_hgnc(input_file, version=data_version)

    _, base = _version_base(ms, data_version, "hgnc")
    _emit(ms, output_format, output, base)


@hgnc.command("symbols")
@click.argument("complete_set_file", type=ExistingPathType, required=False)
@_opt_output
@_opt_version
@_opt_label_format
def hgnc_symbols(
    complete_set_file: Path | None,
    output: Path | None,
    data_version: str | None,
    output_format: str,
) -> None:
    """Generate HGNC symbol mappings (previous/alias to current symbol).

    Formats: symbol_sec2pri (dict), pri_symbols (set of current symbols),
    sssom/rdf/json/owl/all.
    """
    from pysec2pri.api import generate_hgnc_symbols

    ms = generate_hgnc_symbols(complete_set_file, version=data_version)
    _, base = _version_base(ms, data_version, "hgnc_symbols")
    _emit(ms, output_format, output, base)


@main.group()
def ncbi() -> None:
    """Parse NCBI Gene files and generate mappings."""


@ncbi.command("ids")
@click.argument("input_file", type=ExistingPathType, required=False)
@_opt_output
@_opt_tax_id
@_opt_version
@_opt_id_format
def ncbi_ids(
    input_file: Path | None,
    output: Path | None,
    tax_id: str,
    data_version: str | None,
    output_format: str,
) -> None:
    """Generate NCBI Gene ID mappings (discontinued to current Gene IDs).

    Formats: sec2pri (dict), pri_ids (set of all current IDs), sssom/rdf/json/owl/all.
    """
    from pysec2pri.api import generate_ncbi

    ms = generate_ncbi(input_file, tax_id=tax_id, version=data_version)
    _, base = _version_base(ms, data_version, "ncbi")
    _emit(ms, output_format, output, base)


@ncbi.command("symbols")
@click.argument("input_file", type=ExistingPathType, required=False)
@_opt_output
@_opt_tax_id
@_opt_version
@_opt_label_format
def ncbi_symbols(
    input_file: Path | None,
    output: Path | None,
    tax_id: str,
    data_version: str | None,
    output_format: str,
) -> None:
    """Generate NCBI Gene symbol mappings (previous to current gene symbols).

    Formats: symbol_sec2pri (dict), pri_symbols (set of current symbols),
    sssom/rdf/json/owl/all.
    """
    from pysec2pri.api import generate_ncbi_symbols

    ms = generate_ncbi_symbols(input_file, tax_id=tax_id, version=data_version)
    _, base = _version_base(ms, data_version, "ncbi_symbols")
    _emit(ms, output_format, output, base)


@main.command()
@click.argument("input_file", type=ExistingPathType, required=False)
@_opt_output
@_opt_version
@_opt_id_format
@click.option("--delac-file", type=ExistingPathType)
def uniprot(
    input_file: Path | None,
    output: Path | None,
    data_version: str | None,
    output_format: str,
    delac_file: Path | None,
) -> None:
    """Parse UniProt secondary accessions and generate mappings."""
    from pysec2pri.api import generate_uniprot

    ms = generate_uniprot(input_path=input_file, delac_file=delac_file, version=data_version)
    _, base = _version_base(ms, data_version, "uniprot")
    _emit(ms, output_format, output, base)


@main.group()
def wikidata() -> None:
    """Query Wikidata SPARQL for redirect mappings."""


@wikidata.command("ids")
@click.argument("input_file", type=ExistingPathType, required=False)
@_opt_output
@_opt_id_format
@_opt_entity_type
@click.option("--test-subset", is_flag=True, help="Use test queries (LIMIT 10)")
def wikidata_ids(
    input_file: Path | None,
    output: Path | None,
    output_format: str,
    entity_type: str | None,
    test_subset: bool,
) -> None:
    """Generate Wikidata ID mappings (redirected  to current Wikidata QIDs).

    Formats: sec2pri (dict), pri_ids (set of current QIDs), sssom/rdf/json/owl/all.
    """
    from pysec2pri.api import generate_wikidata

    for etype in [entity_type] if entity_type else ["metabolites", "genes", "proteins"]:
        action = f"Parsing {input_file}" if input_file else "Querying Wikidata"
        click.echo(f"{action} for {etype}{'  (test subset)' if test_subset else ''}...")
        ms = generate_wikidata(input_path=input_file, entity_type=etype, test_subset=test_subset)
        _, base = _version_base(ms, None, f"wikidata_{etype}")
        _emit(ms, output_format, output, base)


@wikidata.command("symbols")
@click.argument("input_file", type=ExistingPathType, required=False)
@_opt_output
@_opt_label_format
@_opt_entity_type
@click.option("--test-subset", is_flag=True, help="Use test queries (LIMIT 10)")
def wikidata_symbols(
    input_file: Path | None,
    output: Path | None,
    output_format: str,
    entity_type: str | None,
    test_subset: bool,
) -> None:
    """Generate Wikidata label mappings (previous label  to current label).

    Formats: symbol_sec2pri (dict), pri_symbols (set of current labels),
    sssom/rdf/json/owl/all.
    """
    from pysec2pri.api import generate_wikidata_symbols

    for etype in [entity_type] if entity_type else ["metabolites", "genes", "proteins"]:
        action = f"Parsing {input_file}" if input_file else "Querying Wikidata"
        click.echo(f"{action} for {etype}{'  (test subset)' if test_subset else ''}...")
        ms = generate_wikidata_symbols(
            input_path=input_file, entity_type=etype, test_subset=test_subset
        )
        _, base = _version_base(ms, None, f"wikidata_{etype}_symbols")
        _emit(ms, output_format, output, base)


@main.command(name="all")
@click.option("-o", "--output-dir", type=PathType, default=Path("."), help="Output directory")
@click.option(
    "--datasources",
    default="chebi,hgnc,hmdb,ncbi,uniprot,wikidata",
    help="Comma-separated list of datasources",
)
def export_all(output_dir: Path, datasources: str) -> None:
    """Export all formats for specified datasources."""
    from pysec2pri.api import (
        combine_mapping_sets,
        generate_chebi,
        generate_hgnc,
        generate_hmdb,
        generate_hmdb_proteins,
        generate_ncbi,
        generate_uniprot,
        generate_wikidata,
        save,
    )
    from pysec2pri.download import CloudflareBlockedError
    from pysec2pri.parsers.base import Sec2PriMappingSet

    _generators: dict[str, Callable[[], Sec2PriMappingSet]] = {
        "chebi": generate_chebi,
        "hgnc": generate_hgnc,
        "ncbi": generate_ncbi,
        "uniprot": generate_uniprot,
        "wikidata": generate_wikidata,
    }

    for ds in [d.strip().lower() for d in datasources.split(",")]:
        click.echo(f"\n=== Processing {ds.upper()} ===")
        try:
            if ds == "hmdb":
                ms = combine_mapping_sets(generate_hmdb(), generate_hmdb_proteins())
            elif ds in _generators:
                ms = _generators[ds]()
            else:
                click.echo(f"  Unknown datasource: {ds}", err=True)
                continue
        except CloudflareBlockedError as e:
            click.echo(f"  WARNING: Skipping {ds.upper()} - {e}", err=True)
            continue
        except Exception as e:
            click.echo(f"  WARNING: Skipping {ds.upper()} - {e}", err=True)
            continue

        version = getattr(ms, "mapping_set_version", None)
        base = f"{ds}{'_' + version if version else ''}"
        ds_dir = output_dir / (f"{ds}_{version}" if version else ds)
        path = save(ms, "all", ds_dir, base_name=base)
        click.echo(f"  Wrote all formats to {path}/")


@main.command("update-ids")
@click.argument("input_file", type=ExistingPathType)
@click.argument(
    "datasource", type=click.Choice(["chebi", "hgnc", "hmdb", "ncbi", "uniprot", "wikidata"])
)
@click.option(
    "--at",
    "columns",
    required=True,
    multiple=True,
    metavar="COLUMN",
    help="Column name(s) containing identifiers to resolve. Repeat for multiple columns.",
)
@click.option("-o", "--output", "output_path", type=PathType, help="Output file path (TSV or CSV).")
@click.option("--suffix", default="_primary", show_default=True, help="Suffix for new columns.")
@click.option("--sep", default=None, help="Delimiter (inferred from extension if omitted).")
@click.option(
    "--mapping",
    "mapping_file",
    type=ExistingPathType,
    default=None,
    help="Pre-built sec2pri TSV mapping file to use instead of regenerating.",
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
    data_version: str | None,
    no_progress: bool,
) -> None:
    """Resolve secondary IDs in INPUT_FILE to primary IDs using DATASOURCE mappings.

    For each column specified with --at, a new column ``<col><suffix>`` is
    added to the output containing the resolved primary identifiers.
    Identifiers not found in the mapping are kept unchanged.

    Pass --mapping to skip downloading/regenerating the mapping set and use
    an existing sec2pri TSV file instead.

    Example::

        pysec2pri update-ids my_genes.tsv hgnc --at gene_id -o my_genes_updated.tsv
        pysec2pri update-ids my_genes.tsv hgnc --at gene_id --mapping hgnc_sec2pri.tsv
    """
    from pysec2pri.api import (
        generate_chebi,
        generate_hgnc,
        generate_hmdb,
        generate_ncbi,
        generate_uniprot,
        generate_wikidata,
        load_mapping,
        resolve_ids,
    )

    if mapping_file is not None:
        click.echo(f"Loading mappings from {mapping_file}...")
        ms = load_mapping(mapping_file)
    else:
        _generators = {
            "chebi": lambda: generate_chebi(version=data_version, show_progress=not no_progress),
            "hgnc": lambda: generate_hgnc(version=data_version, show_progress=not no_progress),
            "hmdb": lambda: generate_hmdb(version=data_version, show_progress=not no_progress),
            "ncbi": lambda: generate_ncbi(version=data_version, show_progress=not no_progress),
            "uniprot": lambda: generate_uniprot(
                version=data_version, show_progress=not no_progress
            ),
            "wikidata": lambda: generate_wikidata(
                version=data_version, show_progress=not no_progress
            ),
        }
        click.echo(f"Loading {datasource.upper()} mappings...")
        ms = _generators[datasource]()
    _resolve_and_print(input_file, ms, resolve_ids, columns, output_path, suffix, sep)


@main.command("update-symbols")
@click.argument("input_file", type=ExistingPathType)
@click.argument("datasource", type=click.Choice(["chebi", "hgnc", "ncbi", "wikidata"]))
@click.option(
    "--at",
    "columns",
    required=True,
    multiple=True,
    metavar="COLUMN",
    help="Column name(s) containing symbols to resolve. Repeat for multiple columns.",
)
@click.option("-o", "--output", "output_path", type=PathType, help="Output file path (TSV or CSV).")
@click.option("--suffix", default="_current", show_default=True, help="Suffix for new columns.")
@click.option("--sep", default=None, help="Delimiter (inferred from extension if omitted).")
@click.option(
    "--mapping",
    "mapping_file",
    type=ExistingPathType,
    default=None,
    help="Pre-built symbol2prev TSV mapping file to use instead of regenerating.",
)
@_opt_tax_id
@_opt_entity_type
@_opt_subset
@_opt_version
@_opt_no_progress
def update_symbols_cmd(
    input_file: Path,
    datasource: str,
    columns: tuple[str, ...],
    output_path: Path | None,
    suffix: str,
    sep: str | None,
    mapping_file: Path | None,
    tax_id: str,
    entity_type: str | None,
    subset: str,
    data_version: str | None,
    no_progress: bool,
) -> None:
    """Resolve previous/alias labels in INPUT_FILE to current labels using DATASOURCE.

    For each column specified with --at, a new column ``<col><suffix>`` is
    added containing the resolved current labels.  Labels not found in the
    mapping are kept unchanged.

    Pass --mapping to skip downloading/regenerating the mapping set and use
    an existing symbol2prev TSV file instead.

    Example::

        pysec2pri update-symbols my_genes.tsv hgnc --at symbol -o my_genes_updated.tsv
        pysec2pri update-symbols my_genes.tsv hgnc --at symbol --mapping hgnc_symbol2prev.tsv
    """
    from pysec2pri.api import (
        generate_chebi_synonyms,
        generate_hgnc_symbols,
        generate_ncbi_symbols,
        generate_wikidata_symbols,
        load_label_mapping,
        resolve_symbols,
    )

    if mapping_file is not None:
        click.echo(f"Loading symbol mappings from {mapping_file}...")
        ms = load_label_mapping(mapping_file)
    else:
        click.echo(f"Loading {datasource.upper()} symbol mappings...")
        if datasource == "chebi":
            ms = generate_chebi_synonyms(
                version=data_version, subset=subset, show_progress=not no_progress
            )
        elif datasource == "hgnc":
            ms = generate_hgnc_symbols(version=data_version, show_progress=not no_progress)
        elif datasource == "ncbi":
            ms = generate_ncbi_symbols(
                tax_id=tax_id, version=data_version, show_progress=not no_progress
            )
        else:
            ms = generate_wikidata_symbols(
                entity_type=entity_type, version=data_version, show_progress=not no_progress
            )
    _resolve_and_print(input_file, ms, resolve_symbols, columns, output_path, suffix, sep)


if __name__ == "__main__":
    main()
