"""Command-line interface for pysec2pri."""

from __future__ import annotations

import warnings
from collections.abc import Callable
from pathlib import Path
from typing import Any

import click

from pysec2pri.constants import OutputType, get_output_directory_name
from pysec2pri.logging import set_log_level
from pysec2pri.models import MappingSet

# Suppress FutureWarning from sssom library (pandas downcasting deprecation)
warnings.filterwarnings("ignore", category=FutureWarning, module="sssom")

# Type alias for click.Path with path_type - workaround for mypy type-var issue
# See: https://github.com/pallets/click/issues/2558
PathType: Any = click.Path(path_type=Path)  # type: ignore[type-var]
ExistingPathType: Any = click.Path(  # type: ignore[type-var]
    exists=True, path_type=Path
)


# === Helpers for CLI commands ===
def _subset_gzip_file(path: Path, output_dir: Path, max_lines: int) -> Path:
    """Create a subset of a gzip file with limited lines."""
    import gzip
    from pathlib import Path

    subset_path = Path(output_dir) / path.name.replace(".gz", "_testsubset.gz")
    with (
        gzip.open(path, "rt", encoding="utf-8") as fin,
        gzip.open(subset_path, "wt", encoding="utf-8") as fout,
    ):
        for i, line in enumerate(fin):
            if i >= max_lines:
                break
            fout.write(line)
    return subset_path


def _subset_zip_file(path: Path, output_dir: Path, max_lines: int) -> Path:
    """Create a subset of a zip file with limited lines per entry."""
    import zipfile
    from pathlib import Path

    subset_path = Path(output_dir) / path.name.replace(".zip", "_testsubset.zip")
    with (
        zipfile.ZipFile(path, "r") as zin,
        zipfile.ZipFile(subset_path, "w") as zout,
    ):
        for name in zin.namelist():
            with zin.open(name) as fin:
                data = b"".join(line for i, line in enumerate(fin) if i < max_lines)
                zout.writestr(name, data)
    return subset_path


def _subset_text_file(path: Path, output_dir: Path, max_lines: int) -> Path:
    """Create a subset of a text file with limited lines."""
    from pathlib import Path

    subset_path = Path(output_dir) / path.name.replace(path.suffix, f"_testsubset{path.suffix}")
    with (
        path.open("r", encoding="utf-8") as fin,
        subset_path.open("w", encoding="utf-8") as fout,
    ):
        for i, line in enumerate(fin):
            if i >= max_lines:
                break
            fout.write(line)
    return subset_path


def _download_and_subset(
    datasource: str,
    output_dir: Path,
    test_subset: bool = False,
    subset_lines: int = 1000,
    subset_gz_lines: int = 1000,
    subset_zip_lines: int = 1000,
    file_keys: list[str] | None = None,
) -> dict[str, Path]:
    """Download files for a datasource and optionally subset them."""
    from pysec2pri.download import download_datasource

    if file_keys is None:
        file_keys = []
    downloaded_files = download_datasource(datasource, output_dir)

    if not (test_subset and downloaded_files):
        return downloaded_files

    subsetted = {}
    for key, path in downloaded_files.items():
        if file_keys and key not in file_keys:
            continue

        if path.suffix == ".gz":
            subsetted[key] = _subset_gzip_file(path, output_dir, subset_gz_lines)
        elif path.suffix == ".zip":
            subsetted[key] = _subset_zip_file(path, output_dir, subset_zip_lines)
        else:
            subsetted[key] = _subset_text_file(path, output_dir, subset_lines)

    downloaded_files.update(subsetted)
    return downloaded_files


def _write_outputs_and_report(
    mapping_set: MappingSet,
    output_dir: Path,
    outputs: list[str] | None,
    datasource_name: str,
    mapping_date: str | None,
    click_echo: Callable[..., Any] = print,
) -> dict[str, Path]:
    from pysec2pri.api import write_outputs

    written = write_outputs(
        mapping_set,
        output_dir,
        output_types=outputs,
        datasource_name=datasource_name,
        mapping_date=mapping_date,
    )
    click_echo(f"Generated {len(written)} output file(s) in {output_dir}:")
    for output_type, path in written.items():
        click_echo(f"  {output_type}: {path.name}")
    click_echo(f"Total mappings: {len(mapping_set)}")
    return written


def _cleanup_downloaded(
    downloaded_files: dict[str, Path] | None,
    keep_download: bool,
    click_echo: Callable[..., Any] = print,
) -> None:
    if not keep_download and downloaded_files:
        for path in downloaded_files.values():
            if path.exists():
                path.unlink()
                click_echo(f"Cleaned up: {path}")


"""Command line interface for :mod:`pysec2pri`.

This CLI provides commands for parsing various biological database formats
and generating SSSOM (Simple Standard for Sharing Ontology Mappings) files
as well as legacy output formats.

By default, commands auto-download source files from upstream databases.
Use input options to specify local files instead.

Output files are placed in a date-stamped folder: {database}{ddmmyy}/
"""


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


def validate_outputs(
    ctx: click.Context, param: click.Parameter, value: tuple[str, ...]
) -> list[str] | None:
    """Validate and convert output type values."""
    if not value:
        return None
    valid_types = {ot.value for ot in OutputType}
    for v in value:
        if v not in valid_types:
            raise click.BadParameter(
                f"Invalid output type '{v}'. Valid types: {', '.join(valid_types)}"
            )
    return list(value)


# Common options as decorators
def common_options(f: Callable[..., Any]) -> Callable[..., Any]:
    """Add common options to database commands."""
    f = click.option(
        "--output-dir",
        type=PathType,
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
    f = click.option(
        "--test-subset",
        is_flag=True,
        help="Process only a small test subset (first lines/records).",
        default=False,
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
    type=ExistingPathType,
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
    test_subset: bool,
) -> None:
    r"""Parse ChEBI SDF file and generate mapping files.

    By default, downloads the latest ChEBI SDF from upstream.
    Use --input to specify a local file instead.

    Available outputs: sssom, priIDs, secID2priID, name2synonym


    Examples:
        # Auto-download and process (all outputs)
        pysec2pri chebi

        # Only generate SSSOM and priIDs
        pysec2pri chebi --outputs sssom --outputs priIDs

        # Use local file with custom output directory
        pysec2pri chebi --input ChEBI_complete_3star.sdf --output-dir ./results
    """
    from pysec2pri.api import parse_chebi

    # Determine output directory
    if output_dir is None:
        output_dir = Path(".") / get_output_directory_name("chebi")
    downloaded_files = None
    if input_file is None:
        click.echo("Downloading ChEBI SDF file...")
        downloaded_files = _download_and_subset(
            "chebi", output_dir, test_subset, subset_lines=50000, file_keys=["sdf"]
        )
        input_file = downloaded_files.get("sdf")
        if input_file is None:
            raise click.ClickException("Failed to download ChEBI SDF file")
        click.echo(f"Downloaded: {input_file}")
    click.echo(f"Parsing {input_file}...")
    mapping_set = parse_chebi(input_file, version=source_version)
    _write_outputs_and_report(
        mapping_set, output_dir, outputs, "chebi", mapping_date, click_echo=click.echo
    )
    _cleanup_downloaded(downloaded_files, keep_download, click_echo=click.echo)


# =============================================================================
# HMDB command
# =============================================================================


@main.command()
@click.option(
    "-i",
    "--input",
    "input_file",
    type=ExistingPathType,
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
    test_subset: bool,
) -> None:
    r"""Parse HMDB XML files and generate mapping files.

    By default, downloads the latest HMDB metabolites ZIP from upstream.
    Use --input to specify a local file instead.

    Available outputs: sssom, priIDs, secID2priID, name2synonym


    Examples:
        # Auto-download and process (all outputs)
        pysec2pri hmdb

        # Only generate SSSOM
        pysec2pri hmdb --outputs sssom

        # Use local file
        pysec2pri hmdb --input hmdb_metabolites.zip
    """
    from pysec2pri.api import parse_hmdb

    if output_dir is None:
        output_dir = Path(".") / get_output_directory_name("hmdb")
    downloaded_files = None
    if input_file is None:
        click.echo("Downloading HMDB metabolites file...")
        downloaded_files = _download_and_subset(
            "hmdb", output_dir, test_subset, file_keys=["metabolites"]
        )
        input_file = downloaded_files.get("metabolites")
        if input_file is None:
            raise click.ClickException("Failed to download HMDB file")
        click.echo(f"Downloaded: {input_file}")
    click.echo(f"Parsing {input_file}...")
    mapping_set = parse_hmdb(input_file, version=source_version)
    _write_outputs_and_report(
        mapping_set, output_dir, outputs, "hmdb", mapping_date, click_echo=click.echo
    )
    _cleanup_downloaded(downloaded_files, keep_download, click_echo=click.echo)


# =============================================================================
# HGNC command
# =============================================================================


@main.command()
@click.option(
    "-w",
    "--withdrawn",
    "withdrawn_file",
    type=ExistingPathType,
    help="Local HGNC withdrawn file. Downloads from upstream if not set.",
)
@click.option(
    "-c",
    "--complete-set",
    type=ExistingPathType,
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
    test_subset: bool,
) -> None:
    r"""Parse HGNC files and generate mapping files.

    By default, downloads the latest HGNC files from upstream.
    Use --withdrawn and --complete-set to specify local files.

    Available outputs: sssom, priIDs, secID2priID, symbol2prev


    Examples:
        # Auto-download and process (all outputs)
        pysec2pri hgnc

        # Only generate SSSOM
        pysec2pri hgnc --outputs sssom

        # Use local files
        pysec2pri hgnc --withdrawn withdrawn.txt --complete-set hgnc_complete.txt
    """
    from pysec2pri.api import parse_hgnc

    if output_dir is None:
        output_dir = Path(".") / get_output_directory_name("hgnc")
    downloaded_files = None
    if withdrawn_file is None:
        click.echo("Downloading HGNC files...")
        downloaded_files = _download_and_subset(
            "hgnc", output_dir, test_subset, file_keys=["withdrawn", "complete"]
        )
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
    _write_outputs_and_report(
        mapping_set, output_dir, outputs, "hgnc", mapping_date, click_echo=click.echo
    )
    _cleanup_downloaded(downloaded_files, keep_download, click_echo=click.echo)


# =============================================================================
# NCBI command
# =============================================================================


@main.command()
@click.option(
    "-H",
    "--history",
    "history_file",
    type=ExistingPathType,
    help="Local gene_history file. If not provided, downloads from upstream.",
)
@click.option(
    "-I",
    "--gene-info",
    type=ExistingPathType,
    help="Local gene_info file. Downloaded if --history is not set.",
)
@click.option(
    "-t",
    "--tax-id",
    default="9606",
    help="Taxonomy ID to filter (default: 9606 for human, others unsupported).",
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
    test_subset: bool,
) -> None:
    r"""Parse NCBI Gene files and generate mapping files.

    By default, downloads the latest NCBI Gene files from upstream.
    Use --history and --gene-info to specify local files.

    Available outputs: sssom, priIDs, secID2priID, symbol2prev


    Examples:
        # Auto-download and process (human genes, all outputs)
        pysec2pri ncbi

        # Different organism
        pysec2pri ncbi --tax-id 10090

        # Use local files
        pysec2pri ncbi --history gene_history.gz --gene-info gene_info.gz
    """
    from pysec2pri.api import parse_ncbi

    if output_dir is None:
        output_dir = Path(".") / get_output_directory_name("ncbi")
    downloaded_files = None
    if history_file is None:
        click.echo("Downloading NCBI Gene files...")
        downloaded_files = _download_and_subset(
            "ncbi", output_dir, test_subset, file_keys=["gene_history", "gene_info"]
        )
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
    _write_outputs_and_report(
        mapping_set, output_dir, outputs, "ncbi", mapping_date, click_echo=click.echo
    )
    _cleanup_downloaded(downloaded_files, keep_download, click_echo=click.echo)


# =============================================================================
# UniProt command
# =============================================================================


@main.command()
@click.option(
    "-s",
    "--sec-ac",
    "sec_ac_file",
    type=ExistingPathType,
    help="Local sec_ac.txt file. If not provided, downloads from upstream.",
)
@click.option(
    "--delac",
    type=ExistingPathType,
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
    test_subset: bool,
) -> None:
    r"""Parse UniProt files and generate mapping files.

    By default, downloads the latest UniProt files from upstream.
    Use --sec-ac and --delac to specify local files.

    Available outputs: sssom, priIDs, secID2priID


    Examples:
        # Auto-download and process (all outputs)
        pysec2pri uniprot

        # Only generate SSSOM
        pysec2pri uniprot --outputs sssom

        # Use local files
        pysec2pri uniprot --sec-ac sec_ac.txt --delac delac_sp.txt
    """
    from pysec2pri.api import parse_uniprot

    if output_dir is None:
        output_dir = Path(".") / get_output_directory_name("uniprot")
    downloaded_files = None
    if sec_ac_file is None:
        click.echo("Downloading UniProt files...")
        downloaded_files = _download_and_subset(
            "uniprot", output_dir, test_subset, file_keys=["secondary", "deleted"]
        )
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
    _write_outputs_and_report(
        mapping_set, output_dir, outputs, "uniprot", mapping_date, click_echo=click.echo
    )
    _cleanup_downloaded(downloaded_files, keep_download, click_echo=click.echo)


# =============================================================================
# Wikidata command
# =============================================================================


@main.command()
@click.option(
    "-t",
    "--type",
    "entity_type",
    type=click.Choice(["all", "metabolites", "genes", "proteins"]),
    default="all",
    help="Type of Wikidata entities to query. Use 'all' to query all types.",
)
@click.option(
    "--output-dir",
    type=PathType,
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
@click.option(
    "--test-subset",
    is_flag=True,
    help="Process only a small test subset (uses _test.rq query).",
    default=False,
)
def wikidata(
    entity_type: str,
    output_dir: Path | None,
    outputs: list[str] | None,
    source_version: str | None,
    mapping_date: str | None,
    test_subset: bool,
) -> None:
    r"""Query Wikidata for redirect mappings.

    Queries the QLever Wikidata endpoint for owl:sameAs redirect mappings.

    Available outputs: sssom, secID2priID


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

    if entity_type == "all":
        types_to_query = ["metabolites", "genes", "proteins"]
    else:
        types_to_query = [entity_type]
    all_written = {}
    all_mappings = []
    click.echo(f"Will query Wikidata for {', '.join(types_to_query)}")
    for etype in types_to_query:
        click.echo(f"Querying Wikidata for {etype}...{' (test subset)' if test_subset else ''}")
        mapping_set = parse_wikidata(
            entity_type=etype,
            version=source_version,
            test_subset=test_subset,
        )
        written = write_outputs(
            mapping_set,
            f"{output_dir}_{etype}",
            output_types=outputs,
            datasource_name="wikidata",
            mapping_date=mapping_date,
        )
        all_written.update({f"{etype}_{k}": v for k, v in written.items()})
        all_mappings.append(mapping_set)
    click.echo(f"Generated {len(all_written)} output file(s) in {output_dir}:")
    for output_type, path in all_written.items():
        click.echo(f"  {output_type}: {path.name}")
    total = sum(len(ms) for ms in all_mappings)
    click.echo(f"Total mappings: {total}")


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
    type=PathType,
    default=".",
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
    r"""Download source files for a datasource (no processing).

    Downloads the required input files from the upstream database.


    Examples:
        # Download ChEBI files
        pysec2pri download chebi -o ./data

        # Download all datasources
        pysec2pri download all -o ./data
    """
    from pysec2pri.download import download_datasource

    datasources = (
        ["chebi", "hmdb", "hgnc", "ncbi", "uniprot"] if datasource == "all" else [datasource]
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
    r"""Check if a new release is available for a datasource.

    Examples:
        # Check ChEBI for updates
        pysec2pri check-release chebi

        # Check all datasources
        pysec2pri check-release all
    """
    from pysec2pri.download import get_latest_release_info

    datasources = (
        ["chebi", "hmdb", "hgnc", "ncbi", "uniprot"] if datasource == "all" else [datasource]
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
                    click.echo(f"New release available! ({current_version} -> {info.version})")
                else:
                    click.echo("No new release available.")
        except Exception as e:
            click.echo(f"Error checking {ds}: {e}", err=True)


# =============================================================================
# Diff command
# =============================================================================


@main.command()
@click.argument("old_file", type=ExistingPathType)
@click.argument("new_file", type=ExistingPathType)
@click.option(
    "-d",
    "--datasource",
    default="unknown",
    help="Name of the datasource.",
)
@click.option(
    "-o",
    "--output",
    type=PathType,
    help="Output file for diff report.",
)
def diff(
    old_file: Path,
    new_file: Path,
    datasource: str,
    output: Path | None,
) -> None:
    r"""Compare two SSSOM files and show differences.

    OLD_FILE: Path to the older SSSOM file.
    NEW_FILE: Path to the newer SSSOM file.


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
