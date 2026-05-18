"""Command-line interface for pysec2pri."""

from __future__ import annotations

import warnings
from pathlib import Path

import click

warnings.filterwarnings("ignore", category=FutureWarning, module="sssom")

PathType = click.Path(path_type=Path)  # type: ignore[type-var]
ExistingPathType = click.Path(exists=True, path_type=Path)  # type: ignore[type-var]


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


@main.command()
@click.argument("input_file", type=ExistingPathType, required=False)
@click.option("-o", "--output", type=PathType, help="Output file or directory")
@click.option("--version", "data_version", help="Release number (e.g. 245)")
@click.option(
    "--subset",
    default="3star",
    type=click.Choice(["3star", "complete"]),
    help="Compound subset",
)
@click.option(
    "--format",
    "output_format",
    default="sssom",
    type=click.Choice(["sssom", "sec2pri", "pri_ids", "name2synonym", "rdf", "json", "owl", "all"]),
)
@click.option(
    "--mapping-sets",
    default="ids",
    type=click.Choice(["ids", "synonyms", "all"]),
    help="Which mappings to include",
)
def chebi(
    input_file: Path | None,
    output: Path | None,
    data_version: str | None,
    subset: str,
    output_format: str,
    mapping_sets: str,
) -> None:
    """Parse ChEBI data and generate mappings."""
    from pysec2pri.api import generate_chebi, save

    ms = generate_chebi(
        input_path=input_file,
        version=data_version,
        subset=subset,
        mapping_sets=mapping_sets,
    )
    version = getattr(ms, "mapping_set_version", None) or data_version
    subset_sfx = "_3star" if subset == "3star" else ""
    base = f"chebi{subset_sfx}{'_' + version if version else ''}"
    path = save(ms, output_format, output, base_name=base)
    click.echo(f"Wrote {len(ms.mappings or [])} mappings to {path}")


@main.command()
@click.option("--metabolites-file", type=ExistingPathType, default=None)
@click.option("--proteins-file", type=ExistingPathType, default=None)
@click.option("--metabolites-only", is_flag=True, default=False)
@click.option("--proteins-only", is_flag=True, default=False)
@click.option("-o", "--output", type=PathType, help="Output file or directory")
@click.option("--version", "data_version")
@click.option(
    "--format",
    "output_format",
    default="sssom",
    type=click.Choice(["sssom", "sec2pri", "pri_ids", "rdf", "json", "owl", "all"]),
)
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
    from pysec2pri.api import generate_hmdb, generate_hmdb_proteins, save

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
        version = getattr(ms, "mapping_set_version", None) or data_version
        base = f"hmdb_metabolites{'_' + version if version else ''}"
        path = save(ms, output_format, output, base_name=base)
        click.echo(f"  Wrote {len(ms.mappings or [])} metabolite mappings to {path}")

    if want_proteins:
        click.echo("Parsing HMDB proteins...")
        ms = generate_hmdb_proteins(proteins_file, version=data_version)
        version = getattr(ms, "mapping_set_version", None) or data_version
        base = f"hmdb_proteins{'_' + version if version else ''}"
        path = save(ms, output_format, output, base_name=base)
        click.echo(f"  Wrote {len(ms.mappings or [])} protein mappings to {path}")


@main.command()
@click.argument("input_file", type=ExistingPathType, required=False)
@click.option("-o", "--output", type=PathType)
@click.option("--version", "data_version")
@click.option(
    "--format",
    "output_format",
    default="sssom",
    type=click.Choice(["sssom", "sec2pri", "symbol2prev", "pri_ids", "rdf", "json", "owl", "all"]),
)
@click.option("--symbols-file", type=ExistingPathType)
@click.option(
    "--status",
    "statuses",
    multiple=True,
    default=["Approved"],
    show_default=True,
    help="Entry status to include. Repeat for multiple. Use --status='' for all.",
)
def hgnc(
    input_file: Path | None,
    output: Path | None,
    data_version: str | None,
    output_format: str,
    symbols_file: Path | None,
    statuses: tuple[str, ...],
) -> None:
    """Parse HGNC files and generate mappings."""
    from pysec2pri.api import generate_hgnc, generate_hgnc_symbols, save

    status_filter: list[str] | None = None if not statuses or statuses == ("",) else list(statuses)

    if output_format == "symbol2prev":
        ms = generate_hgnc_symbols(symbols_file, version=data_version, statuses=status_filter)
        if status_filter is None:
            status_sfx, note = "_all_statuses", "Includes all entry statuses."
        elif status_filter != ["Approved"]:
            joined = "-".join(s.lower().replace(" ", "_") for s in status_filter)
            status_sfx, note = f"_{joined}", f"Filtered to statuses: {', '.join(status_filter)}."
        else:
            status_sfx, note = "", "Filtered to Approved entries only."
        existing = getattr(ms, "comment", None) or ""
        ms.comment = f"{existing} {note}".strip() if existing else note
    else:
        ms = generate_hgnc(input_file, version=data_version)
        status_sfx = ""

    version = getattr(ms, "mapping_set_version", None) or data_version
    base = f"hgnc{status_sfx}{'_' + version if version else ''}"
    path = save(ms, output_format, output, base_name=base)
    click.echo(f"Wrote {len(ms.mappings or [])} mappings to {path}")


@main.command()
@click.argument("input_file", type=ExistingPathType, required=False)
@click.option("-o", "--output", type=PathType)
@click.option("--tax-id", default="9606", help="Taxonomy ID (default: 9606)")
@click.option("--version", "data_version")
@click.option(
    "--format",
    "output_format",
    default="sssom",
    type=click.Choice(["sssom", "sec2pri", "symbol2prev", "pri_ids", "rdf", "json", "owl", "all"]),
)
@click.option("--symbols-file", type=ExistingPathType)
def ncbi(
    input_file: Path | None,
    output: Path | None,
    tax_id: str,
    data_version: str | None,
    output_format: str,
    symbols_file: Path | None,
) -> None:
    """Parse NCBI Gene files and generate mappings."""
    from pysec2pri.api import generate_ncbi, generate_ncbi_symbols, save

    if output_format == "symbol2prev":
        ms = generate_ncbi_symbols(symbols_file, tax_id=tax_id, version=data_version)
    else:
        ms = generate_ncbi(input_file, tax_id=tax_id, version=data_version)

    version = getattr(ms, "mapping_set_version", None) or data_version
    base = f"ncbi{'_' + version if version else ''}"
    path = save(ms, output_format, output, base_name=base)
    click.echo(f"Wrote {len(ms.mappings or [])} mappings to {path}")


@main.command()
@click.argument("input_file", type=ExistingPathType, required=False)
@click.option("-o", "--output", type=PathType)
@click.option("--version", "data_version")
@click.option(
    "--format",
    "output_format",
    default="sssom",
    type=click.Choice(["sssom", "sec2pri", "pri_ids", "rdf", "json", "owl", "all"]),
)
@click.option("--delac-file", type=ExistingPathType)
def uniprot(
    input_file: Path | None,
    output: Path | None,
    data_version: str | None,
    output_format: str,
    delac_file: Path | None,
) -> None:
    """Parse UniProt secondary accessions and generate mappings."""
    from pysec2pri.api import generate_uniprot, save

    ms = generate_uniprot(input_path=input_file, delac_file=delac_file, version=data_version)
    version = getattr(ms, "mapping_set_version", None) or data_version
    base = f"uniprot{'_' + version if version else ''}"
    path = save(ms, output_format, output, base_name=base)
    click.echo(f"Wrote {len(ms.mappings or [])} mappings to {path}")


@main.command()
@click.argument("input_file", type=ExistingPathType, required=False)
@click.option("-o", "--output", type=PathType)
@click.option(
    "--format",
    "output_format",
    default="sssom",
    type=click.Choice(["sssom", "sec2pri", "pri_ids", "rdf", "json", "owl", "all"]),
)
@click.option(
    "--entity-type",
    default=None,
    type=click.Choice(["metabolites", "chemicals", "genes", "proteins"]),
    help="Entity type to query. Queries all if omitted.",
)
@click.option("--test-subset", is_flag=True, help="Use test queries (LIMIT 10)")
def wikidata(
    input_file: Path | None,
    output: Path | None,
    output_format: str,
    entity_type: str | None,
    test_subset: bool,
) -> None:
    """Query Wikidata SPARQL for redirect mappings."""
    from pysec2pri.api import generate_wikidata, save

    entity_types = [entity_type] if entity_type else ["metabolites", "genes", "proteins"]

    for etype in entity_types:
        action = f"Parsing {input_file}" if input_file else "Querying Wikidata"
        click.echo(f"{action} for {etype}{'  (test subset)' if test_subset else ''}...")
        ms = generate_wikidata(input_path=input_file, entity_type=etype, test_subset=test_subset)
        version = getattr(ms, "mapping_set_version", None)
        base = f"wikidata_{etype}{'_' + version if version else ''}"
        path = save(ms, output_format, output, base_name=base)
        click.echo(f"  Wrote {len(ms.mappings or [])} {etype} mappings to {path}")


@main.command(name="all")
@click.option("-o", "--output-dir", type=PathType, default=Path("."), help="Output directory")
@click.option(
    "--datasources",
    default="chebi,hgnc,hmdb,ncbi,uniprot,wikidata",
    help="Comma-separated list of datasources",
)
def export_all(output_dir: Path, datasources: str) -> None:
    """Export all formats for specified datasources."""
    from collections.abc import Callable

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
@click.option("--version", "data_version", default=None, help="Datasource release version.")
@click.option("--no-progress", is_flag=True, default=False, help="Suppress progress bars.")
def update_ids_cmd(
    input_file: Path,
    datasource: str,
    columns: tuple[str, ...],
    output_path: Path | None,
    suffix: str,
    sep: str | None,
    data_version: str | None,
    no_progress: bool,
) -> None:
    """Resolve secondary IDs in INPUT_FILE to primary IDs using DATASOURCE mappings.

    For each column specified with --at, a new column ``<col><suffix>`` is
    added to the output containing the resolved primary identifiers.
    Identifiers not found in the mapping are kept unchanged.

    Example::

        pysec2pri update-ids my_genes.tsv hgnc --at gene_id -o my_genes_updated.tsv
    """
    from pysec2pri.api import (
        generate_chebi,
        generate_hgnc,
        generate_hmdb,
        generate_ncbi,
        generate_uniprot,
        generate_wikidata,
        resolve_ids,
    )

    _generators = {
        "chebi": lambda: generate_chebi(version=data_version, show_progress=not no_progress),
        "hgnc": lambda: generate_hgnc(version=data_version, show_progress=not no_progress),
        "hmdb": lambda: generate_hmdb(version=data_version, show_progress=not no_progress),
        "ncbi": lambda: generate_ncbi(version=data_version, show_progress=not no_progress),
        "uniprot": lambda: generate_uniprot(version=data_version, show_progress=not no_progress),
        "wikidata": lambda: generate_wikidata(version=data_version, show_progress=not no_progress),
    }

    click.echo(f"Loading {datasource.upper()} mappings...")
    ms = _generators[datasource]()

    at: str | list[str] = list(columns) if len(columns) > 1 else columns[0]
    click.echo(f"Resolving IDs in column(s): {', '.join(columns)}")

    df = resolve_ids(input_file, ms, at=at, output_path=output_path, suffix=suffix, sep=sep)

    if output_path:
        click.echo(f"Wrote {len(df)} rows to {output_path}")
    else:
        click.echo(df.to_csv(sep="\t" if (sep or "\t") == "\t" else ",", index=False))


if __name__ == "__main__":
    main()
