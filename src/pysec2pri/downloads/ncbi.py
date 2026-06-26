"""Downloading and release detection for NCBI Gene."""

from __future__ import annotations

from mapkgsutils.download import ReleaseInfo, get_file_last_modified

from pysec2pri.constants import ALL_DATASOURCES

__all__ = ["check_ncbi_release"]


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
