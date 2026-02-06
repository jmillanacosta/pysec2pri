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

    from pysec2pri.diff import MappingDiff
    from pysec2pri.parsers.base import Sec2PriMappingSet

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


def _download_chebi_files(
    input_file: Path | None,
    version: str | None,
    subset: str,
    use_sdf: bool = False,
) -> dict[str, Any]:
    """Download ChEBI files based on version.

    For version >= 245 (and use_sdf=False): downloads TSV flat files
    For version < 245 (or use_sdf=True): downloads SDF file

    Args:
        input_file: User-provided input file (optional).
        version: Version/release number (optional).
        subset: "3star" or "complete".
        use_sdf: Force SDF format even for releases >= 245.

    Returns:
        Dictionary with keys:
            - "format": "tsv" or "sdf"
            - For TSV: "secondary_ids", "names", "compounds" paths
            - For SDF: "sdf" path
    """
    from pysec2pri.download import check_chebi_release
    from pysec2pri.parsers.base import ChEBIDownloader

    # If user provided a file, assume it's SDF (legacy behavior)
    if input_file is not None:
        return {"format": "sdf", "sdf": input_file}

    # Determine version to use
    if version is None:
        release_info = check_chebi_release()
        version = release_info.version or "245"
        click.echo(f"Latest ChEBI release: {version}")

    # Use the ChEBIDownloader to handle format detection and download
    downloader = ChEBIDownloader(
        version=version,
        subset=subset,
        use_sdf=use_sdf,
    )

    # Get format and URLs
    data_format = downloader.get_format(version)
    urls = downloader.get_download_urls(version)

    click.echo(f"Downloading ChEBI {version} ({data_format.upper()} format)...")
    for key, url in urls.items():
        click.echo(f"  {key}: {url}")

    tmpdir = Path(tempfile.mkdtemp(prefix=f"pysec2pri_chebi_{version}_"))
    downloaded = downloader.download(tmpdir, version=version)

    if data_format == "tsv":
        return {
            "format": "tsv",
            "secondary_ids": downloaded.get("secondary_ids"),
            "names": downloaded.get("names"),
            "compounds": downloaded.get("compounds"),
        }
    else:
        return {"format": "sdf", "sdf": downloaded.get("sdf")}


@click.group()
@click.version_option()
def main() -> None:
    """pysec2pri: Secondary to Primary ID mapping tool."""


@main.command()
@click.argument("file1", type=ExistingPathType)
@click.argument("file2", type=ExistingPathType)
@click.option(
    "-o",
    "--output",
    type=PathType,
    help="Output file for diff results (TSV format)",
)
@click.option(
    "--show-all",
    is_flag=True,
    default=False,
    help="Show all differences (not just summary)",
)
@click.option(
    "--datasource",
    default="unknown",
    help="Datasource name for the diff summary",
)
def diff(
    file1: Path,
    file2: Path,
    output: Path | None,
    show_all: bool,
    datasource: str,
) -> None:
    """Compare two SSSOM mapping files and show differences.

    FILE1 is the first (old/baseline) SSSOM file.
    FILE2 is the second (new/comparison) SSSOM file.

    Shows added, removed, and changed mappings between the two files.
    """
    from pysec2pri.api import write_diff_output
    from pysec2pri.diff import diff_sssom_files, summarize_diff

    click.echo(f"Comparing {file1.name} vs {file2.name}...")

    result = diff_sssom_files(file1, file2, datasource=datasource)

    # Print summary
    click.echo("")
    click.echo(summarize_diff(result))

    # Show detailed differences if requested or if counts are manageable
    if show_all or result.total_changes <= 50:
        _show_diff_details(result)
    elif result.total_changes > 50:
        click.echo("")
        click.echo(
            f"(Showing summary only - {result.total_changes} changes. "
            "Use --show-all to see all differences)"
        )

    # Write to output file if specified
    if output:
        if result.total_changes > 0:
            write_diff_output(result, output)
            click.echo(f"\nWrote diff to {output}")
        else:
            click.echo("\nNo differences to write.")


def _show_diff_details(result: MappingDiff) -> None:
    """Show detailed diff output for added/removed/changed mappings."""
    if result.added_count > 0:
        click.echo("")
        click.echo(f"=== Added mappings ({result.added_count}) ===")
        click.echo("object_id\tsubject_id")
        for row in result.added.iter_rows(named=True):
            click.echo(f"{row['object_id']}\t{row['subject_id']}")

    if result.removed_count > 0:
        click.echo("")
        click.echo(f"=== Removed mappings ({result.removed_count}) ===")
        click.echo("object_id\tsubject_id")
        for row in result.removed.iter_rows(named=True):
            click.echo(f"{row['object_id']}\t{row['subject_id']}")

    if result.changed_count > 0:
        click.echo("")
        click.echo(f"=== Changed mappings ({result.changed_count}) ===")
        click.echo("object_id\told_subject_id\tnew_subject_id")
        for row in result.changed.iter_rows(named=True):
            click.echo(f"{row['object_id']}\t{row['old_subject_id']}\t{row['new_subject_id']}")


@main.command()
@click.argument("input_file", type=ExistingPathType, required=False)
@click.option("-o", "--output", type=PathType, help="Output file path")
@click.option(
    "--version",
    "data_version",
    help="Data version string (e.g. 245)",
)
@click.option(
    "--subset",
    default="3star",
    type=click.Choice(["3star", "complete"]),
    help="Compound subset: 3star (default) or complete",
)
@click.option(
    "--use-sdf",
    is_flag=True,
    default=False,
    help="Force SDF format even for releases >= 245 (default: auto-detect)",
)
@click.option(
    "--format",
    "output_format",
    default="sssom",
    type=click.Choice(["sssom", "sec2pri", "pri_ids", "name2synonym", "all"]),
    help="Output format",
)
@click.option(
    "--mapping-sets",
    default="all",
    type=click.Choice(["ids", "synonyms", "all"]),
    help="Which mapping sets to include: ids (sec2pri), synonyms, or all",
)
def chebi(
    input_file: Path | None,
    output: Path | None,
    data_version: str | None,
    subset: str,
    use_sdf: bool,
    output_format: str,
    mapping_sets: str,
) -> None:
    """Parse ChEBI data files and generate mappings.

    Automatically detects format based on version:
    - Releases >= 245: Uses TSV flat files (secondary_ids.tsv, names.tsv)
    - Releases < 245: Uses SDF files (chebi_3_stars.sdf or ChEBI_complete.sdf)

    Use --use-sdf to force SDF format even for releases >= 245.

    Use --mapping-sets to control which mappings are included:
    - ids: Only secondary ID → primary ID mappings
    - synonyms: Only name → synonym mappings
    - all: Both mapping sets combined (default)
    """
    # Download files and determine format
    files = _download_chebi_files(input_file, data_version, subset, use_sdf)

    # Parse requested mapping sets using API functions
    id_mappings, synonym_mappings = _parse_chebi_mapping_sets(
        files, data_version, subset, mapping_sets
    )

    # Combine mapping sets if needed
    result = _combine_chebi_results(id_mappings, synonym_mappings, mapping_sets)

    # Get version for output naming
    version = getattr(result, "mapping_set_version", None) or data_version
    base_name = _get_chebi_base_name(subset, version)

    # Write output in requested format
    _write_chebi_output(result, output, output_format, base_name, version)


def _parse_chebi_mapping_sets(
    files: dict[str, Any],
    version: str | None,
    subset: str,
    mapping_sets: str,
) -> tuple[Sec2PriMappingSet | None, Sec2PriMappingSet | None]:
    """Parse ChEBI mapping sets based on format and mapping_sets option."""
    from pysec2pri.api import parse_chebi, parse_chebi_synonyms

    id_mappings = None
    synonym_mappings = None

    if files["format"] == "tsv":
        # New TSV format (>= 245)
        if mapping_sets in ("ids", "all"):
            id_mappings = parse_chebi(
                version=version,
                subset=subset,
                secondary_ids_path=files["secondary_ids"],
                compounds_path=files.get("compounds"),
            )
        if mapping_sets in ("synonyms", "all"):
            synonym_mappings = parse_chebi_synonyms(
                version=version,
                subset=subset,
                names_path=files["names"],
                compounds_path=files.get("compounds"),
            )
    else:
        # Legacy SDF format (< 245)
        if mapping_sets in ("ids", "all"):
            id_mappings = parse_chebi(files["sdf"], version=version, subset=subset)
        if mapping_sets in ("synonyms", "all"):
            synonym_mappings = parse_chebi_synonyms(files["sdf"], version=version, subset=subset)

    return id_mappings, synonym_mappings


def _combine_chebi_results(
    id_mappings: Sec2PriMappingSet | None,
    synonym_mappings: Sec2PriMappingSet | None,
    mapping_sets: str,
) -> Sec2PriMappingSet:
    """Combine mapping sets and report counts."""
    from pysec2pri.api import combine_mapping_sets

    if mapping_sets == "all" and id_mappings and synonym_mappings:
        id_count = len(id_mappings.mappings or [])
        syn_count = len(synonym_mappings.mappings or [])
        result = combine_mapping_sets(id_mappings, synonym_mappings)
        total = len(result.mappings or [])
        click.echo(
            f"Combined {id_count} ID mappings + {syn_count} synonym mappings = {total} total"
        )
        return result
    elif id_mappings:
        return id_mappings
    elif synonym_mappings:
        return synonym_mappings
    else:
        raise click.ClickException("No mapping sets generated")


def _get_chebi_base_name(subset: str, version: str | None) -> str:
    """Generate base name for ChEBI output files."""
    subset_sfx = "_3star" if subset == "3star" else ""
    v_suffix = f"_{version}" if version else ""
    return f"chebi{subset_sfx}{v_suffix}"


def _write_chebi_output(
    result: Sec2PriMappingSet,
    output: Path | None,
    output_format: str,
    base_name: str,
    version: str | None,
) -> None:
    """Write ChEBI output in requested format."""
    from pysec2pri.api import write_all_formats
    from pysec2pri.exports import (
        write_name2synonym,
        write_pri_ids,
        write_sec2pri,
        write_sssom,
    )

    if output_format == "all":
        datasource = base_name.split("_")[0]
        out_dir = _resolve_all_format_dir(output, datasource, version)
        write_all_formats(result, out_dir, base_name)
        click.echo(f"Wrote all formats to {out_dir}/")
    else:
        out_path = output or Path(f"{base_name}_{output_format}.tsv")
        if output_format == "sssom":
            write_sssom(result, out_path)
        elif output_format == "sec2pri":
            write_sec2pri(result, out_path)
        elif output_format == "pri_ids":
            write_pri_ids(result, out_path)
        elif output_format == "name2synonym":
            write_name2synonym(result, out_path)
        n = len(result.mappings or [])
        click.echo(f"Wrote {n} mappings to {out_path}")


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
