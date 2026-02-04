"""Command-line interface for pysec2pri."""

from __future__ import annotations

import tempfile
import warnings
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any

import click

if TYPE_CHECKING:
    from sssom_schema import MappingSet

warnings.filterwarnings("ignore", category=FutureWarning, module="sssom")

PathType: Any = click.Path(path_type=Path)  # type: ignore[type-var]
ExistingPathType: Any = click.Path(
    exists=True,
    path_type=Path,  # type: ignore[type-var]
)


def _download_if_needed(
    datasource: str,
    input_file: Path | None,
    file_key: str = "main",
    version: str | None = None,
) -> Path:
    """Download datasource file if not provided."""
    if input_file is not None:
        return input_file

    from pysec2pri.download import download_datasource, get_download_urls

    version_str = f" version {version}" if version else ""
    click.echo(f"No input file provided. Downloading {datasource}{version_str}...")

    # Show URLs that will be downloaded
    urls = get_download_urls(datasource, version=version)
    for key, url in urls.items():
        click.echo(f"  {key}: {url}")

    tmpdir = Path(tempfile.mkdtemp(prefix=f"pysec2pri_{datasource}_"))
    downloaded = download_datasource(datasource, tmpdir, version=version)

    if file_key in downloaded:
        return downloaded[file_key]
    # Return first available file
    return next(iter(downloaded.values()))


def _get_output_filename(
    datasource: str,
    output_format: str,
    version: str | None = None,
) -> Path:
    """Generate output filename with optional version."""
    if version:
        return Path(f"{datasource}_{version}_{output_format}.tsv")
    return Path(f"{datasource}_{output_format}.tsv")


def _get_output_dir(
    datasource: str,
    version: str | None = None,
) -> Path:
    """Generate output directory name for 'all' format."""
    if version:
        return Path(f"{datasource}_{version}")
    return Path(datasource)


def _resolve_all_format_dir(
    output: Path | None,
    datasource: str,
    version: str | None,
) -> Path:
    """Resolve output directory for 'all' format.

    If output is provided and looks like a file (has .tsv extension),
    use its parent directory plus versioned folder.
    If output is a directory, use it directly.
    If output is None, generate default directory name.
    """
    if output is not None:
        # If output looks like a file path (has extension), use parent dir
        if output.suffix:
            return output.parent / _get_output_dir(datasource, version)
        # Output is a directory - use it directly without adding subfolder
        return output
    # No output specified - create versioned folder in current directory
    return _get_output_dir(datasource, version)


@click.group()
@click.version_option()
def main() -> None:
    """pysec2pri: Secondary to Primary ID mapping tool."""


@main.command()
@click.argument("input_file", type=ExistingPathType, required=False)
@click.option("-o", "--output", type=PathType, help="Output file path")
@click.option(
    "--version",
    "data_version",
    help="Data version string (e.g. 245)",
)
@click.option(
    "--format",
    "output_format",
    default="sssom",
    type=click.Choice(["sssom", "sec2pri", "pri_ids", "name2synonym", "all"]),
    help="Output format",
)
def chebi(
    input_file: Path | None,
    output: Path | None,
    data_version: str | None,
    output_format: str,
) -> None:
    """Parse ChEBI SDF file and generate mappings."""
    from pysec2pri.api import parse_chebi, parse_chebi_synonyms
    from pysec2pri.exports import (
        write_name2synonym,
        write_pri_ids,
        write_sec2pri,
        write_sssom,
    )

    input_file = _download_if_needed("chebi", input_file, "sdf", version=data_version)

    if output_format == "name2synonym":
        result = parse_chebi_synonyms(input_file, version=data_version)
    else:
        result = parse_chebi(input_file, version=data_version)

    # Get version from result for output naming
    version = getattr(result, "mapping_set_version", None) or data_version
    v_suffix = f"_{version}" if version else ""

    if output_format == "all":
        out_dir = _resolve_all_format_dir(output, "chebi", version)
        out_dir.mkdir(parents=True, exist_ok=True)
        write_sssom(result, out_dir / f"chebi{v_suffix}_sssom.tsv")
        write_sec2pri(result, out_dir / f"chebi{v_suffix}_sec2pri.tsv")
        write_pri_ids(result, out_dir / f"chebi{v_suffix}_pri_ids.tsv")
        name2syn = out_dir / f"chebi{v_suffix}_name2synonym.tsv"
        write_name2synonym(result, name2syn)
        click.echo(f"Wrote all formats to {out_dir}/")
    else:
        out_path = output or _get_output_filename("chebi", output_format, version)
        if output_format == "sssom":
            write_sssom(result, out_path)
        elif output_format == "sec2pri":
            write_sec2pri(result, out_path)
        elif output_format == "pri_ids":
            write_pri_ids(result, out_path)
        elif output_format == "name2synonym":
            write_name2synonym(result, out_path)
        click.echo(f"Wrote {len(result.mappings or [])} mappings to {out_path}")


@main.command()
@click.argument("input_file", type=ExistingPathType, required=False)
@click.option("-o", "--output", type=PathType, help="Output file path")
@click.option(
    "--format",
    "output_format",
    default="sssom",
    type=click.Choice(["sssom", "sec2pri", "pri_ids", "all"]),
    help="Output format",
)
def hmdb(
    input_file: Path | None,
    output: Path | None,
    output_format: str,
) -> None:
    """Parse HMDB XML file and generate mappings."""
    from pysec2pri.api import parse_hmdb
    from pysec2pri.exports import write_pri_ids, write_sec2pri, write_sssom

    input_file = _download_if_needed("hmdb", input_file, "metabolites")

    result = parse_hmdb(input_file)
    version = getattr(result, "mapping_set_version", None)
    v_suffix = f"_{version}" if version else ""

    if output_format == "all":
        out_dir = _resolve_all_format_dir(output, "hmdb", version)
        out_dir.mkdir(parents=True, exist_ok=True)
        write_sssom(result, out_dir / f"hmdb{v_suffix}_sssom.tsv")
        write_sec2pri(result, out_dir / f"hmdb{v_suffix}_sec2pri.tsv")
        write_pri_ids(result, out_dir / f"hmdb{v_suffix}_pri_ids.tsv")
        click.echo(f"Wrote all formats to {out_dir}/")
    else:
        out_path = output or _get_output_filename("hmdb", output_format, version)
        if output_format == "sssom":
            write_sssom(result, out_path)
        elif output_format == "sec2pri":
            write_sec2pri(result, out_path)
        elif output_format == "pri_ids":
            write_pri_ids(result, out_path)
        click.echo(f"Wrote {len(result.mappings or [])} mappings to {out_path}")


@main.command()
@click.argument("input_file", type=ExistingPathType, required=False)
@click.option("-o", "--output", type=PathType, help="Output file path")
@click.option(
    "--version",
    "data_version",
    help="Data version string (eg 2023-07-01)",
)
@click.option(
    "--format",
    "output_format",
    default="sssom",
    type=click.Choice(["sssom", "sec2pri", "symbol2prev", "pri_ids", "all"]),
    help="Output format",
)
@click.option(
    "--symbols-file",
    type=ExistingPathType,
    help="HGNC complete set file for symbol mappings",
)
def hgnc(
    input_file: Path | None,
    output: Path | None,
    data_version: str | None,
    output_format: str,
    symbols_file: Path | None,
) -> None:
    """Parse HGNC withdrawn file and generate mappings."""
    from pysec2pri.api import parse_hgnc, parse_hgnc_symbols
    from pysec2pri.exports import (
        write_pri_ids,
        write_sec2pri,
        write_sssom,
        write_symbol2prev,
    )

    if output_format == "symbol2prev":
        if symbols_file is None:
            symbols_file = _download_if_needed("hgnc", None, "complete", version=data_version)
        result = parse_hgnc_symbols(symbols_file, version=data_version)
    else:
        input_file = _download_if_needed("hgnc", input_file, "withdrawn", version=data_version)
        result = parse_hgnc(input_file, version=data_version)

    version = getattr(result, "mapping_set_version", None) or data_version
    v_suffix = f"_{version}" if version else ""

    if output_format == "all":
        out_dir = _resolve_all_format_dir(output, "hgnc", version)
        out_dir.mkdir(parents=True, exist_ok=True)
        write_sssom(result, out_dir / f"hgnc{v_suffix}_sssom.tsv")
        write_sec2pri(result, out_dir / f"hgnc{v_suffix}_sec2pri.tsv")
        write_pri_ids(result, out_dir / f"hgnc{v_suffix}_pri_ids.tsv")
        write_symbol2prev(result, out_dir / f"hgnc{v_suffix}_symbol2prev.tsv")
        click.echo(f"Wrote all formats to {out_dir}/")
    else:
        out_path = output or _get_output_filename("hgnc", output_format, version)
        if output_format == "sssom":
            write_sssom(result, out_path)
        elif output_format == "sec2pri":
            write_sec2pri(result, out_path)
        elif output_format == "symbol2prev":
            write_symbol2prev(result, out_path)
        elif output_format == "pri_ids":
            write_pri_ids(result, out_path)
        click.echo(f"Wrote {len(result.mappings or [])} mappings to {out_path}")


@main.command()
@click.argument("input_file", type=ExistingPathType, required=False)
@click.option("-o", "--output", type=PathType, help="Output file path")
@click.option("--tax-id", default="9606", help="Taxonomy ID (default: 9606)")
@click.option(
    "--format",
    "output_format",
    default="sssom",
    type=click.Choice(["sssom", "sec2pri", "symbol2prev", "pri_ids", "all"]),
    help="Output format",
)
@click.option(
    "--symbols-file",
    type=ExistingPathType,
    help="Gene info file for symbol mappings",
)
def ncbi(
    input_file: Path | None,
    output: Path | None,
    tax_id: str,
    output_format: str,
    symbols_file: Path | None,
) -> None:
    """Parse NCBI Gene history file and generate mappings."""
    from pysec2pri.api import parse_ncbi, parse_ncbi_symbols
    from pysec2pri.exports import (
        write_pri_ids,
        write_sec2pri,
        write_sssom,
        write_symbol2prev,
    )

    if output_format == "symbol2prev":
        if symbols_file is None:
            symbols_file = _download_if_needed("ncbi", None, "gene_info")
        result = parse_ncbi_symbols(symbols_file, tax_id=tax_id)
    else:
        input_file = _download_if_needed("ncbi", input_file, "gene_history")
        result = parse_ncbi(input_file, tax_id=tax_id)

    version = getattr(result, "mapping_set_version", None)
    v_suffix = f"_{version}" if version else ""

    if output_format == "all":
        out_dir = _resolve_all_format_dir(output, "ncbi", version)
        out_dir.mkdir(parents=True, exist_ok=True)
        write_sssom(result, out_dir / f"ncbi{v_suffix}_sssom.tsv")
        write_sec2pri(result, out_dir / f"ncbi{v_suffix}_sec2pri.tsv")
        write_pri_ids(result, out_dir / f"ncbi{v_suffix}_pri_ids.tsv")
        write_symbol2prev(result, out_dir / f"ncbi{v_suffix}_symbol2prev.tsv")
        click.echo(f"Wrote all formats to {out_dir}/")
    else:
        out_path = output or _get_output_filename("ncbi", output_format, version)
        if output_format == "sssom":
            write_sssom(result, out_path)
        elif output_format == "sec2pri":
            write_sec2pri(result, out_path)
        elif output_format == "symbol2prev":
            write_symbol2prev(result, out_path)
        elif output_format == "pri_ids":
            write_pri_ids(result, out_path)
        click.echo(f"Wrote {len(result.mappings or [])} mappings to {out_path}")


@main.command()
@click.argument("input_file", type=ExistingPathType, required=False)
@click.option("-o", "--output", type=PathType, help="Output file path")
@click.option(
    "--version",
    "data_version",
    help="Release version (e.g., 2024_01)",
)
@click.option(
    "--format",
    "output_format",
    default="sssom",
    type=click.Choice(["sssom", "sec2pri", "pri_ids", "all"]),
    help="Output format",
)
@click.option(
    "--delac-file",
    type=ExistingPathType,
    help="Deleted accessions file (delac_sp.txt)",
)
def uniprot(
    input_file: Path | None,
    output: Path | None,
    data_version: str | None,
    output_format: str,
    delac_file: Path | None,
) -> None:
    """Parse UniProt secondary accessions file and generate mappings."""
    from pysec2pri.api import parse_uniprot
    from pysec2pri.exports import write_pri_ids, write_sec2pri, write_sssom

    input_file = _download_if_needed("uniprot", input_file, "sec_ac", version=data_version)

    result = parse_uniprot(input_file, delac_file=delac_file)
    version = getattr(result, "mapping_set_version", None) or data_version
    v_suffix = f"_{version}" if version else ""

    if output_format == "all":
        out_dir = _resolve_all_format_dir(output, "uniprot", version)
        out_dir.mkdir(parents=True, exist_ok=True)
        write_sssom(result, out_dir / f"uniprot{v_suffix}_sssom.tsv")
        write_sec2pri(result, out_dir / f"uniprot{v_suffix}_sec2pri.tsv")
        write_pri_ids(result, out_dir / f"uniprot{v_suffix}_pri_ids.tsv")
        click.echo(f"Wrote all formats to {out_dir}/")
    else:
        out_path = output or _get_output_filename("uniprot", output_format, version)
        if output_format == "sssom":
            write_sssom(result, out_path)
        elif output_format == "sec2pri":
            write_sec2pri(result, out_path)
        elif output_format == "pri_ids":
            write_pri_ids(result, out_path)
        click.echo(f"Wrote {len(result.mappings or [])} mappings to {out_path}")


@main.command()
@click.argument("input_file", type=ExistingPathType, required=False)
@click.option("-o", "--output", type=PathType, help="Output file path")
@click.option(
    "--format",
    "output_format",
    default="sssom",
    type=click.Choice(["sssom", "sec2pri", "pri_ids", "all"]),
    help="Output format",
)
@click.option(
    "--entity-type",
    default=None,
    type=click.Choice(["metabolites", "chemicals", "genes", "proteins"]),
    help="Entity type to query. If not specified, queries all types.",
)
@click.option(
    "--test-subset",
    is_flag=True,
    help="Use test queries with LIMIT 10 (for testing)",
)
def wikidata(
    input_file: Path | None,
    output: Path | None,
    output_format: str,
    entity_type: str | None,
    test_subset: bool,
) -> None:
    """Query Wikidata SPARQL endpoint for redirect mappings.

    By default, queries ALL entity types (metabolites, genes, proteins)
    defined in the config. Use --entity-type to query a specific type.

    Each entity type is saved to a separate file with appropriate naming.

    If INPUT_FILE is provided, parses pre-downloaded TSV results instead
    of querying the SPARQL endpoint.

    Uses the QLever endpoint for better performance.
    """
    from pysec2pri.api import parse_wikidata
    from pysec2pri.exports import write_pri_ids, write_sec2pri, write_sssom

    # Determine which entity types to process
    if entity_type:
        entity_types = [entity_type]
    else:
        entity_types = ["metabolites", "genes", "proteins"]

    subset_msg = " (test subset)" if test_subset else ""

    for etype in entity_types:
        if input_file is None:
            click.echo(f"Querying Wikidata for {etype}{subset_msg}...")
        else:
            click.echo(f"Parsing {input_file} for {etype}...")

        result = parse_wikidata(
            input_file,
            entity_type=etype,
            test_subset=test_subset,
        )

        version = getattr(result, "mapping_set_version", None)
        v_suffix = f"_{version}" if version else ""
        out_prefix = f"wikidata_{etype}"

        if output_format == "all":
            out_dir = _resolve_all_format_dir(output, "wikidata", version)
            out_dir.mkdir(parents=True, exist_ok=True)
            sssom_f = out_dir / f"{out_prefix}{v_suffix}_sssom.tsv"
            sec2pri_f = out_dir / f"{out_prefix}{v_suffix}_sec2pri.tsv"
            pri_ids_f = out_dir / f"{out_prefix}{v_suffix}_pri_ids.tsv"
            write_sssom(result, sssom_f)
            write_sec2pri(result, sec2pri_f)
            write_pri_ids(result, pri_ids_f)
            n = len(result.mappings or [])
            click.echo(f"  Wrote {n} {etype} mappings to {out_dir}/")
        else:
            # Determine output path
            if output is None:
                out_path = _get_output_filename(out_prefix, output_format, version)
            elif output.suffix == ".tsv":
                # User provided a specific file - use parent as dir
                out_dir = output.parent
                out_dir.mkdir(parents=True, exist_ok=True)
                fname = f"{out_prefix}{v_suffix}_{output_format}.tsv"
                out_path = out_dir / fname
            else:
                # User provided a directory
                out_dir = output
                out_dir.mkdir(parents=True, exist_ok=True)
                fname = f"{out_prefix}{v_suffix}_{output_format}.tsv"
                out_path = out_dir / fname

            if output_format == "sssom":
                write_sssom(result, out_path)
            elif output_format == "sec2pri":
                write_sec2pri(result, out_path)
            elif output_format == "pri_ids":
                write_pri_ids(result, out_path)
            n = len(result.mappings or [])
            click.echo(f"  Wrote {n} {etype} mappings to {out_path}")


def _parse_datasource(ds: str) -> MappingSet | None:
    """Download and parse a datasource, returning the MappingSet or None."""
    from pysec2pri.api import (
        parse_chebi,
        parse_hgnc,
        parse_hmdb,
        parse_ncbi,
        parse_uniprot,
        parse_wikidata,
    )
    from pysec2pri.constants import ALL_DATASOURCES

    config = ALL_DATASOURCES.get(ds.lower())
    if config is None:
        return None

    # Wikidata uses SPARQL queries, no file download
    if ds.lower() == "wikidata":
        return parse_wikidata()

    parsers: dict[str, Callable[..., MappingSet]] = {
        "chebi": parse_chebi,
        "hgnc": parse_hgnc,
        "hmdb": parse_hmdb,
        "ncbi": parse_ncbi,
        "uniprot": parse_uniprot,
    }

    parser = parsers.get(ds.lower())
    if parser is None:
        return None

    file_key = config.primary_file_key or next(iter(config.download_urls), "")
    input_file = _download_if_needed(ds, None, file_key)
    return parser(input_file)


def _get_formats_for_datasource(
    ds: str,
) -> list[tuple[str, Callable[[MappingSet, Path], None]]]:
    """Get the list of export formats applicable to a datasource."""
    from pysec2pri.exports import (
        write_name2synonym,
        write_pri_ids,
        write_sec2pri,
        write_sssom,
        write_symbol2prev,
    )

    base_formats: list[tuple[str, Any]] = [
        ("sssom", write_sssom),
        ("sec2pri", write_sec2pri),
        ("pri_ids", write_pri_ids),
    ]

    if ds == "chebi":
        base_formats.append(("name2synonym", write_name2synonym))
    elif ds in ("hgnc", "ncbi"):
        base_formats.append(("symbol2prev", write_symbol2prev))

    return base_formats


@main.command(name="all")
@click.option(
    "-o",
    "--output-dir",
    type=PathType,
    default=Path("."),
    help="Output directory",
)
@click.option(
    "--datasources",
    default="chebi,hgnc,hmdb,ncbi,uniprot,wikidata",
    help="Comma-separated list of datasources",
)
def export_all(
    output_dir: Path,
    datasources: str,
) -> None:
    """Export all formats for specified datasources."""
    ds_list = [d.strip().lower() for d in datasources.split(",")]

    for ds in ds_list:
        click.echo(f"\n=== Processing {ds.upper()} ===")

        result = _parse_datasource(ds)
        if result is None:
            click.echo(f"Unknown datasource: {ds}")
            continue

        # Get version from result if available
        version = getattr(result, "mapping_set_version", None)

        # Create output folder: datasource_version or just datasource
        if version:
            ds_dir = output_dir / f"{ds}_{version}"
        else:
            ds_dir = output_dir / ds
        ds_dir.mkdir(parents=True, exist_ok=True)

        # Export all formats
        for fmt, writer in _get_formats_for_datasource(ds):
            out_file = ds_dir / _get_output_filename(ds, fmt, version)
            writer(result, out_file)
            click.echo(f"  Wrote {fmt}: {out_file}")


if __name__ == "__main__":
    main()
