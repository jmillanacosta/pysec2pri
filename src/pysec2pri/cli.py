"""Command line interface for :mod:`pysec2pri`.

This CLI provides commands for parsing various biological database formats
and generating SSSOM (Simple Standard for Sharing Ontology Mappings) files
as well as legacy output formats.

By default, commands auto-download source files from upstream databases.
Use input options to specify local files instead.

Output files are placed in a date-stamped folder: {database}{ddmmyy}/
"""

from pathlib import Path

import click

from pysec2pri.constants import OutputType, get_output_directory_name
from pysec2pri.logging import set_log_level

__all__ = [
    "main",
]


def get_available_outputs_for_db(db_name: str) -> list[str]:
    """Get available output types for a database."""
    from pysec2pri.constants import get_datasource_config

    config = get_datasource_config(db_name)
    if config:
        return [ot.value for ot in config.available_outputs]
    return [OutputType.SSSOM.value]


def validate_outputs(ctx, param, value: tuple[str, ...]) -> list[str] | None:
    """Validate and convert output type values."""
    if not value:
        return None
    valid_types = {ot.value for ot in OutputType}
    for v in value:
        if v not in valid_types:
            raise click.BadParameter(
                f"Invalid output type '{v}'. "
                f"Valid types: {', '.join(valid_types)}"
            )
    return list(value)


# Common options as decorators
def common_options(f):
    """Add common options to database commands."""
    f = click.option(
        "--output-dir",
        type=click.Path(path_type=Path),
        default=None,
        help="Base output directory. Default: ./{database}{ddmmyy}/",
    )(f)
    f = click.option(
        "--outputs",
        multiple=True,
        callback=validate_outputs,
        help=(
            "Output types to generate. Can be specified multiple times. "
            "Options: sssom, priIDs, secID2priID, name2synonym, symbol2prev. "
            "Default: all available for the database."
        ),
    )(f)
    f = click.option(
        "-v",
        "--version",
        "source_version",
        help="Version/release of the source database.",
    )(f)
    f = click.option(
        "-d",
        "--date",
        "mapping_date",
        help="Mapping date (YYYY-MM-DD format).",
    )(f)
    f = click.option(
        "--keep-download",
        is_flag=True,
        help="Keep downloaded files after processing.",
    )(f)
    return f


@click.group()
@click.version_option()
@click.option(
    "--log-level",
    type=click.Choice(
        ["debug", "info", "warning", "error", "critical"],
        case_sensitive=False,
    ),
    default="critical",
    help="Set logging verbosity (default: critical).",
)
def main(log_level: str) -> None:
    """pysec2pri: Secondary to Primary Identifier Mapping.

    Convert biological database secondary identifier mappings to SSSOM files
    and legacy output formats.

    By default, commands auto-download source files from upstream databases.
    Output is written to {database}{ddmmyy}/ folder.
    """
    set_log_level(log_level)


# =============================================================================
# ChEBI command
# =============================================================================


@main.command()
@click.option(
    "-i",
    "--input",
    "input_file",
    type=click.Path(exists=True, path_type=Path),
    help="Local ChEBI SDF file. If not provided, downloads from upstream.",
)
@common_options
def chebi(
    input_file: Path | None,
    output_dir: Path | None,
    outputs: list[str] | None,
    source_version: str | None,
    mapping_date: str | None,
    keep_download: bool,
) -> None:
    """Parse ChEBI SDF file and generate mapping files.

    By default, downloads the latest ChEBI SDF from upstream.
    Use --input to specify a local file instead.

    Available outputs: sssom, priIDs, secID2priID, name2synonym

    \b
    Examples:
        # Auto-download and process (all outputs)
        pysec2pri chebi

        # Only generate SSSOM and priIDs
        pysec2pri chebi --outputs sssom --outputs priIDs

        # Use local file with custom output directory
        pysec2pri chebi --input ChEBI_complete_3star.sdf --output-dir ./results
    """
    from pysec2pri.api import parse_chebi, write_outputs
    from pysec2pri.download import download_datasource

    # Determine output directory
    if output_dir is None:
        output_dir = Path(".") / get_output_directory_name("chebi")

    downloaded_files = None

    # Download if no input file provided
    if input_file is None:
        click.echo("Downloading ChEBI SDF file...")
        downloaded_files = download_datasource("chebi", output_dir)
        input_file = downloaded_files.get("sdf")
        if input_file is None:
            raise click.ClickException("Failed to download ChEBI SDF file")
        click.echo(f"Downloaded: {input_file}")

    click.echo(f"Parsing {input_file}...")
    mapping_set = parse_chebi(input_file, version=source_version)

    # Write outputs
    written = write_outputs(
        mapping_set,
        output_dir,
        output_types=outputs,
        datasource_name="chebi",
        mapping_date=mapping_date,
    )

    click.echo(f"Generated {len(written)} output file(s) in {output_dir}:")
    for output_type, path in written.items():
        click.echo(f"  {output_type}: {path.name}")
    click.echo(f"Total mappings: {len(mapping_set)}")

    # Clean up downloaded files if not keeping
    if not keep_download and downloaded_files:
        for path in downloaded_files.values():
            if path.exists():
                path.unlink()
                click.echo(f"Cleaned up: {path}")


# =============================================================================
# HMDB command
# =============================================================================


@main.command()
@click.option(
    "-i",
    "--input",
    "input_file",
    type=click.Path(exists=True, path_type=Path),
    help="Local HMDB ZIP file. If not provided, downloads from upstream.",
)
@common_options
def hmdb(
    input_file: Path | None,
    output_dir: Path | None,
    outputs: list[str] | None,
    source_version: str | None,
    mapping_date: str | None,
    keep_download: bool,
) -> None:
    """Parse HMDB XML files and generate mapping files.

    By default, downloads the latest HMDB metabolites ZIP from upstream.
    Use --input to specify a local file instead.

    Available outputs: sssom, priIDs, secID2priID, name2synonym

    \b
    Examples:
        # Auto-download and process (all outputs)
        pysec2pri hmdb

        # Only generate SSSOM
        pysec2pri hmdb --outputs sssom

        # Use local file
        pysec2pri hmdb --input hmdb_metabolites.zip
    """
    from pysec2pri.api import parse_hmdb, write_outputs
    from pysec2pri.download import download_datasource

    # Determine output directory
    if output_dir is None:
        output_dir = Path(".") / get_output_directory_name("hmdb")

    downloaded_files = None

    # Download if no input file provided
    if input_file is None:
        click.echo("Downloading HMDB metabolites file...")
        downloaded_files = download_datasource("hmdb", output_dir)
        input_file = downloaded_files.get("metabolites")
        if input_file is None:
            raise click.ClickException("Failed to download HMDB file")
        click.echo(f"Downloaded: {input_file}")

    click.echo(f"Parsing {input_file}...")
    mapping_set = parse_hmdb(input_file, version=source_version)

    # Write outputs
    written = write_outputs(
        mapping_set,
        output_dir,
        output_types=outputs,
        datasource_name="hmdb",
        mapping_date=mapping_date,
    )

    click.echo(f"Generated {len(written)} output file(s) in {output_dir}:")
    for output_type, path in written.items():
        click.echo(f"  {output_type}: {path.name}")
    click.echo(f"Total mappings: {len(mapping_set)}")

    # Clean up downloaded files if not keeping
    if not keep_download and downloaded_files:
        for path in downloaded_files.values():
            if path.exists():
                path.unlink()
                click.echo(f"Cleaned up: {path}")


# =============================================================================
# HGNC command
# =============================================================================


@main.command()
@click.option(
    "-w",
    "--withdrawn",
    "withdrawn_file",
    type=click.Path(exists=True, path_type=Path),
    help="Local HGNC withdrawn file. Downloads from upstream if not set.",
)
@click.option(
    "-c",
    "--complete-set",
    type=click.Path(exists=True, path_type=Path),
    help="Local HGNC complete set file. Downloaded if --withdrawn is not set.",
)
@click.option(
    "--include-unmapped",
    is_flag=True,
    help="Include genes with no alias/previous symbols (primary ID only).",
    default=True,
)
@common_options
def hgnc(
    withdrawn_file: Path | None,
    complete_set: Path | None,
    include_unmapped: bool,
    output_dir: Path | None,
    outputs: list[str] | None,
    source_version: str | None,
    mapping_date: str | None,
    keep_download: bool,
) -> None:
    """Parse HGNC files and generate mapping files.

    By default, downloads the latest HGNC files from upstream.
    Use --withdrawn and --complete-set to specify local files.

    Available outputs: sssom, priIDs, secID2priID, symbol2prev

    \b
    Examples:
        # Auto-download and process (all outputs)
        pysec2pri hgnc

        # Only generate SSSOM
        pysec2pri hgnc --outputs sssom

        # Use local files
        pysec2pri hgnc --withdrawn withdrawn.txt --complete-set hgnc_complete.txt
    """
    from pysec2pri.api import parse_hgnc, write_outputs
    from pysec2pri.download import download_datasource

    # Determine output directory
    if output_dir is None:
        output_dir = Path(".") / get_output_directory_name("hgnc")

    downloaded_files = None

    # Download if no input files provided
    if withdrawn_file is None:
        click.echo("Downloading HGNC files...")
        downloaded_files = download_datasource("hgnc", output_dir)
        withdrawn_file = downloaded_files.get("withdrawn")
        complete_set = downloaded_files.get("complete")
        if withdrawn_file is None:
            raise click.ClickException("Failed to download HGNC files")
        click.echo(f"Downloaded withdrawn: {withdrawn_file}")
        if complete_set:
            click.echo(f"Downloaded complete set: {complete_set}")

    click.echo(f"Parsing {withdrawn_file}...")
    mapping_set = parse_hgnc(
        withdrawn_file,
        complete_set_file=complete_set,
        version=source_version,
        include_unmapped_genes=include_unmapped,
    )

    # Write outputs
    written = write_outputs(
        mapping_set,
        output_dir,
        output_types=outputs,
        datasource_name="hgnc",
        mapping_date=mapping_date,
    )

    click.echo(f"Generated {len(written)} output file(s) in {output_dir}:")
    for output_type, path in written.items():
        click.echo(f"  {output_type}: {path.name}")
    click.echo(f"Total mappings: {len(mapping_set)}")

    # Clean up downloaded files if not keeping
    if not keep_download and downloaded_files:
        for path in downloaded_files.values():
            if path.exists():
                path.unlink()
                click.echo(f"Cleaned up: {path}")


# =============================================================================
# NCBI command
# =============================================================================


@main.command()
@click.option(
    "-H",
    "--history",
    "history_file",
    type=click.Path(exists=True, path_type=Path),
    help="Local gene_history file. If not provided, downloads from upstream.",
)
@click.option(
    "-I",
    "--gene-info",
    type=click.Path(exists=True, path_type=Path),
    help="Local gene_info file. Downloaded if --history is not set.",
)
@click.option(
    "-t",
    "--tax-id",
    default="9606",
    help="Taxonomy ID to filter (default: 9606 for human).",
)
@common_options
def ncbi(
    history_file: Path | None,
    gene_info: Path | None,
    tax_id: str,
    output_dir: Path | None,
    outputs: list[str] | None,
    source_version: str | None,
    mapping_date: str | None,
    keep_download: bool,
) -> None:
    """Parse NCBI Gene files and generate mapping files.

    By default, downloads the latest NCBI Gene files from upstream.
    Use --history and --gene-info to specify local files.

    Available outputs: sssom, priIDs, secID2priID, symbol2prev

    \b
    Examples:
        # Auto-download and process (human genes, all outputs)
        pysec2pri ncbi

        # Different organism
        pysec2pri ncbi --tax-id 10090

        # Use local files
        pysec2pri ncbi --history gene_history.gz --gene-info gene_info.gz
    """
    from pysec2pri.api import parse_ncbi, write_outputs
    from pysec2pri.download import download_datasource

    # Determine output directory
    if output_dir is None:
        output_dir = Path(".") / get_output_directory_name("ncbi")

    downloaded_files = None

    # Download if no input files provided
    if history_file is None:
        click.echo("Downloading NCBI Gene files...")
        downloaded_files = download_datasource("ncbi", output_dir)
        history_file = downloaded_files.get("gene_history")
        gene_info = downloaded_files.get("gene_info")
        if history_file is None:
            raise click.ClickException("Failed to download NCBI files")
        click.echo(f"Downloaded history: {history_file}")
        if gene_info:
            click.echo(f"Downloaded gene_info: {gene_info}")

    click.echo(f"Parsing {history_file} (tax_id={tax_id})...")
    mapping_set = parse_ncbi(
        history_file,
        gene_info_file=gene_info,
        tax_id=tax_id,
        version=source_version,
    )

    # Write outputs
    written = write_outputs(
        mapping_set,
        output_dir,
        output_types=outputs,
        datasource_name="ncbi",
        mapping_date=mapping_date,
    )

    click.echo(f"Generated {len(written)} output file(s) in {output_dir}:")
    for output_type, path in written.items():
        click.echo(f"  {output_type}: {path.name}")
    click.echo(f"Total mappings: {len(mapping_set)}")

    # Clean up downloaded files if not keeping
    if not keep_download and downloaded_files:
        for path in downloaded_files.values():
            if path.exists():
                path.unlink()
                click.echo(f"Cleaned up: {path}")


# =============================================================================
# UniProt command
# =============================================================================


@main.command()
@click.option(
    "-s",
    "--sec-ac",
    "sec_ac_file",
    type=click.Path(exists=True, path_type=Path),
    help="Local sec_ac.txt file. If not provided, downloads from upstream.",
)
@click.option(
    "--delac",
    type=click.Path(exists=True, path_type=Path),
    help="Local delac_sp.txt file. Downloaded if --sec-ac is not set.",
)
@common_options
def uniprot(
    sec_ac_file: Path | None,
    delac: Path | None,
    output_dir: Path | None,
    outputs: list[str] | None,
    source_version: str | None,
    mapping_date: str | None,
    keep_download: bool,
) -> None:
    """Parse UniProt files and generate mapping files.

    By default, downloads the latest UniProt files from upstream.
    Use --sec-ac and --delac to specify local files.

    Available outputs: sssom, priIDs, secID2priID

    \b
    Examples:
        # Auto-download and process (all outputs)
        pysec2pri uniprot

        # Only generate SSSOM
        pysec2pri uniprot --outputs sssom

        # Use local files
        pysec2pri uniprot --sec-ac sec_ac.txt --delac delac_sp.txt
    """
    from pysec2pri.api import parse_uniprot, write_outputs
    from pysec2pri.download import download_datasource

    # Determine output directory
    if output_dir is None:
        output_dir = Path(".") / get_output_directory_name("uniprot")

    downloaded_files = None

    # Download if no input files provided
    if sec_ac_file is None:
        click.echo("Downloading UniProt files...")
        downloaded_files = download_datasource("uniprot", output_dir)
        sec_ac_file = downloaded_files.get("secondary")
        delac = downloaded_files.get("deleted")
        if sec_ac_file is None:
            raise click.ClickException("Failed to download UniProt files")
        click.echo(f"Downloaded sec_ac: {sec_ac_file}")
        if delac:
            click.echo(f"Downloaded delac: {delac}")

    click.echo(f"Parsing {sec_ac_file}...")
    mapping_set = parse_uniprot(
        sec_ac_file,
        delac_file=delac,
        version=source_version,
    )

    # Write outputs
    written = write_outputs(
        mapping_set,
        output_dir,
        output_types=outputs,
        datasource_name="uniprot",
        mapping_date=mapping_date,
    )

    click.echo(f"Generated {len(written)} output file(s) in {output_dir}:")
    for output_type, path in written.items():
        click.echo(f"  {output_type}: {path.name}")
    click.echo(f"Total mappings: {len(mapping_set)}")

    # Clean up downloaded files if not keeping
    if not keep_download and downloaded_files:
        for path in downloaded_files.values():
            if path.exists():
                path.unlink()
                click.echo(f"Cleaned up: {path}")


# =============================================================================
# Wikidata command
# =============================================================================


@main.command()
@click.option(
    "-t",
    "--type",
    "entity_type",
    type=click.Choice(["metabolites", "genes", "proteins"]),
    default="metabolites",
    help="Type of Wikidata entities to query.",
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Output directory. Default: ./wikidata{ddmmyy}/",
)
@click.option(
    "--outputs",
    multiple=True,
    callback=validate_outputs,
    help=(
        "Output types to generate. Can be specified multiple times. "
        "Options: sssom, secID2priID. Default: all available."
    ),
)
@click.option(
    "-v",
    "--version",
    "source_version",
    help="Version/date string for the mappings.",
)
@click.option(
    "-d",
    "--date",
    "mapping_date",
    help="Mapping date (YYYY-MM-DD format).",
)
def wikidata(
    entity_type: str,
    output_dir: Path | None,
    outputs: list[str] | None,
    source_version: str | None,
    mapping_date: str | None,
) -> None:
    """Query Wikidata for redirect mappings.

    Queries the QLever Wikidata endpoint for owl:sameAs redirect mappings.

    Available outputs: sssom, secID2priID

    \b
    Examples:
        # Query metabolite redirects
        pysec2pri wikidata --type metabolites

        # Query gene redirects
        pysec2pri wikidata --type genes
    """
    from pysec2pri.api import write_outputs
    from pysec2pri.parsers.wikidata import parse_wikidata

    # Determine output directory
    if output_dir is None:
        output_dir = Path(".") / get_output_directory_name("wikidata")

    click.echo(f"Querying Wikidata for {entity_type} redirects...")
    mapping_set = parse_wikidata(
        entity_type=entity_type,
        version=source_version,
    )

    # Write outputs (None means all available for the datasource)
    written = write_outputs(
        mapping_set,
        output_dir,
        output_types=outputs,
        datasource_name="wikidata",
        mapping_date=mapping_date,
    )

    click.echo(f"Generated {len(written)} output file(s) in {output_dir}:")
    for output_type, path in written.items():
        click.echo(f"  {output_type}: {path.name}")
    click.echo(f"Total mappings: {len(mapping_set)}")


# =============================================================================
# Download command
# =============================================================================


@main.command()
@click.argument(
    "datasource",
    type=click.Choice(["chebi", "hmdb", "hgnc", "ncbi", "uniprot", "all"]),
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(path_type=Path),
    default=Path("."),
    help="Directory to save downloaded files.",
)
@click.option(
    "--no-decompress",
    is_flag=True,
    help="Don't decompress .gz files.",
)
def download(
    datasource: str,
    output_dir: Path,
    no_decompress: bool,
) -> None:
    """Download source files for a datasource (without processing).

    Downloads the required input files from the upstream database.

    \b
    Examples:
        # Download ChEBI files
        pysec2pri download chebi -o ./data

        # Download all datasources
        pysec2pri download all -o ./data
    """
    from pysec2pri.download import download_datasource

    datasources = (
        ["chebi", "hmdb", "hgnc", "ncbi", "uniprot"]
        if datasource == "all"
        else [datasource]
    )

    for ds in datasources:
        click.echo(f"Downloading {ds} files to {output_dir}...")
        try:
            files = download_datasource(
                ds,
                output_dir,
                decompress=not no_decompress,
            )
            click.echo("Downloaded files:")
            for key, path in files.items():
                click.echo(f"  {key}: {path}")
        except Exception as e:
            click.echo(f"Error downloading {ds}: {e}", err=True)


# =============================================================================
# Check release command
# =============================================================================


@main.command("check-release")
@click.argument(
    "datasource",
    type=click.Choice(["chebi", "hmdb", "hgnc", "ncbi", "uniprot", "all"]),
)
@click.option(
    "--current-version",
    help="Current version to compare against.",
)
def check_release_cmd(
    datasource: str,
    current_version: str | None,
) -> None:
    """Check if a new release is available for a datasource.

    \b
    Examples:
        # Check ChEBI for updates
        pysec2pri check-release chebi

        # Check all datasources
        pysec2pri check-release all
    """
    from pysec2pri.download import get_latest_release_info

    datasources = (
        ["chebi", "hmdb", "hgnc", "ncbi", "uniprot"]
        if datasource == "all"
        else [datasource]
    )

    for ds in datasources:
        click.echo(f"\n=== {ds.upper()} ===")
        try:
            info = get_latest_release_info(ds)
            click.echo(f"Latest version: {info.version or 'unknown'}")
            if info.release_date:
                click.echo(f"Release date: {info.release_date.strftime('%Y-%m-%d')}")

            if current_version:
                if info.version and info.version != current_version:
                    click.echo(
                        f"New release available! ({current_version} -> {info.version})"
                    )
                else:
                    click.echo("No new release available.")
        except Exception as e:
            click.echo(f"Error checking {ds}: {e}", err=True)


# =============================================================================
# Diff command
# =============================================================================


@main.command()
@click.argument("old_file", type=click.Path(exists=True, path_type=Path))
@click.argument("new_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-d",
    "--datasource",
    default="unknown",
    help="Name of the datasource.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    help="Output file for diff report.",
)
def diff(
    old_file: Path,
    new_file: Path,
    datasource: str,
    output: Path | None,
) -> None:
    """Compare two SSSOM files and show differences.

    OLD_FILE: Path to the older SSSOM file.
    NEW_FILE: Path to the newer SSSOM file.

    \b
    Examples:
        pysec2pri diff old_mappings.sssom.tsv new_mappings.sssom.tsv
    """
    from pysec2pri.diff import diff_sssom_files, summarize_diff

    result = diff_sssom_files(old_file, new_file, datasource=datasource)
    summary = summarize_diff(result)

    click.echo(summary)

    if output:
        output.write_text(summary)
        click.echo(f"\nDiff report saved to: {output}")


if __name__ == "__main__":
    main()
