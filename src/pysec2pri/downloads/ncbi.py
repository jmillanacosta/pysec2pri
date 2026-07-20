"""Downloading and release detection for NCBI Gene."""

from __future__ import annotations

from typing import TYPE_CHECKING

from mapkgsutils.download import ReleaseInfo, get_file_last_modified

from pysec2pri.constants import ALL_DATASOURCES

if TYPE_CHECKING:
    from pathlib import Path

__all__ = ["check_ncbi_release", "discover_ncbi_species"]


def check_ncbi_release() -> ReleaseInfo:
    """Check for the latest NCBI Gene release.

    Returns:
        ReleaseInfo with the latest NCBI release details.
    """
    history_url = ALL_DATASOURCES["ncbi"].download_urls["gene_history"]
    last_modified = get_file_last_modified(history_url)
    version = last_modified.strftime("%Y-%m-%d") if last_modified else None

    return ReleaseInfo(
        datasource="ncbi",
        version=version,
        release_date=last_modified,
        is_new=True,
        files=dict(ALL_DATASOURCES["ncbi"].download_urls),
    )


def discover_ncbi_species(gene_info_path: Path | str) -> list[str]:
    """Return every distinct NCBI taxon ID present in a downloaded gene_info file.

    Args:
        gene_info_path: Path to a downloaded ``gene_info.gz`` (or
            decompressed ``gene_info``) file.

    Returns:
        Sorted (numerically) list of taxon ID strings found in ``#tax_id``.
    """
    import polars as pl

    df = (
        pl.scan_csv(
            gene_info_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=["-"],
        )
        .select(pl.col("#tax_id").cast(pl.Utf8).unique())
        .collect()
    )
    return sorted(df["#tax_id"].to_list(), key=int)
