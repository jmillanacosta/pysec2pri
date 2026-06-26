"""Download and release detection for biological database sources.

The datasource-agnostic dispatch logic lives in :mod:`mapkgsutils.download`;
this module just binds it to pysec2pri's own datasource registries (see
:mod:`pysec2pri.downloads`).
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

from mapkgsutils.download import (
    CloudflareBlockedError,
    ReleaseInfo,
    download_file,
    get_file_last_modified,
)
from mapkgsutils.download import check_release as _check_release
from mapkgsutils.download import download_datasource as _download_datasource
from mapkgsutils.download import (
    download_datasource_with_release as _download_datasource_with_release,
)
from mapkgsutils.download import get_download_urls as _get_download_urls
from mapkgsutils.download import get_latest_release_info as _get_latest_release_info
from mapkgsutils.download import list_versions as _list_versions
from mapkgsutils.download import resolve_release_date as _resolve_release_date

from pysec2pri.constants import ALL_DATASOURCES
from pysec2pri.downloads import CHECK_RELEASE, DOWNLOADERS, TAR_EXTRACTORS, URLS_AND_DATE
from pysec2pri.downloads.chebi import check_chebi_release
from pysec2pri.downloads.ensembl import check_ensembl_release
from pysec2pri.downloads.hgnc import check_hgnc_release
from pysec2pri.downloads.hmdb import check_hmdb_release
from pysec2pri.downloads.ncbi import check_ncbi_release
from pysec2pri.downloads.uniprot import check_uniprot_release

__all__ = [
    "CloudflareBlockedError",
    "ReleaseInfo",
    "check_chebi_release",
    "check_ensembl_release",
    "check_hgnc_release",
    "check_hmdb_release",
    "check_ncbi_release",
    "check_release",
    "check_uniprot_release",
    "download_datasource",
    "download_datasource_with_release",
    "download_file",
    "get_download_urls",
    "get_file_last_modified",
    "get_latest_release_info",
    "list_versions",
    "resolve_release_date",
]


def get_latest_release_info(datasource: str) -> ReleaseInfo:
    """Get release information for a datasource.

    Args:
        datasource: Name of the datasource (chebi, hmdb, hgnc, ncbi, uniprot).

    Returns:
        ReleaseInfo with the latest release details.

    Raises:
        ValueError: If the datasource is not supported.
    """
    return _get_latest_release_info(datasource, checkers=CHECK_RELEASE)


def check_release(
    datasource: str,
    current_version: str | None = None,
    current_date: datetime | None = None,
) -> ReleaseInfo:
    """Check if a new release is available for a datasource.

    Args:
        datasource: Name of the datasource.
        current_version: Current version string to compare against.
        current_date: Current release date to compare against.

    Returns:
        ReleaseInfo with is_new indicating if update is available.
    """
    return _check_release(datasource, current_version, current_date, checkers=CHECK_RELEASE)


def list_versions(datasource: str) -> Any:
    """List all available archive versions for a datasource.

    Delegates to the datasource's downloader class ``list_versions()``
    method, which contains all source-specific retrieval logic.

    For datasources that publish versioned archives (ChEBI, HGNC, UniProt),
    returns all available version strings sorted in ascending order.

    NCBI and HMDB do not maintain versioned archives; calling this function
    for those datasources raises :class:`ValueError`.

    Args:
        datasource: Datasource name (``"chebi"``, ``"hgnc"``, or
            ``"uniprot"``).

    Returns:
        Sorted list of version strings. Format depends on the datasource:

        - **chebi**: integer release numbers, e.g. ``["200", "201", ..., "245"]``
        - **hgnc**: ISO dates, e.g. ``["2023-01-01", ..., "2026-04-07"]``
        - **uniprot**: release identifiers, e.g. ``["2024_01", "2024_02", ...]``

    Raises:
        ValueError: If the datasource is unknown or has no versioned archive.
    """
    return _list_versions(datasource, known_datasources=ALL_DATASOURCES, downloaders=DOWNLOADERS)


def get_download_urls(
    datasource: str,
    version: str | None = None,
    **kwargs: Any,
) -> dict[str, str]:
    """Get download URLs for a datasource.

    Args:
        datasource: Name of the datasource.
        version: Specific version to get URLs for.
        **kwargs: Datasource-specific knobs (e.g. ``subset`` for ChEBI,
            ``species`` for Ensembl/NCBI) -- which ones apply depends on the
            datasource's own config (see ``DatasourceConfig.subset``/
            ``.species``); irrelevant kwargs are simply ignored.

    Returns:
        Dictionary mapping file keys to URLs.
    """
    return _get_download_urls(
        datasource, version, all_datasources=ALL_DATASOURCES, urls_and_date=URLS_AND_DATE, **kwargs
    )


def resolve_release_date(
    datasource: str,
    version: str | None = None,
    **kwargs: Any,
) -> datetime | None:
    """Resolve the upstream release date for a datasource/version.

    This is the date used for the SSSOM ``mapping_date`` of generated mapping
    sets. It does not download the data files (it may issue a lightweight
    ``HEAD`` request to read a ``Last-Modified`` header). Prefer
    :func:`download_datasource_with_release` when you are downloading anyway,
    to avoid an extra round-trip.

    Args:
        datasource: Name of the datasource.
        version: Specific version, when applicable.
        **kwargs: Datasource-specific knobs (``subset`` for ChEBI, ``species``
            for Ensembl/NCBI); see :func:`get_download_urls`.

    Returns:
        The release date, or None when it cannot be determined.
    """
    return _resolve_release_date(
        datasource, version, all_datasources=ALL_DATASOURCES, urls_and_date=URLS_AND_DATE, **kwargs
    )


def download_datasource(
    datasource: str,
    output_dir: Path,
    decompress: bool = True,
    version: str | None = None,
    keys: list[str] | None = None,
    **kwargs: Any,
) -> dict[str, Path]:
    """Download all files for a datasource.

    For datasources with dynamic URLs (like HGNC quarterly archive),
    this function first checks for the latest release and uses those URLs.
    If a version is specified, it downloads that specific version.

    Args:
        datasource: Name of the datasource.
        output_dir: Directory to save files.
        decompress: Whether to decompress .gz files.
        version: Specific version to download. Format depends on datasource.
        keys: Optional list of file-key names to download. When given, only
            URLs whose key is in this list are fetched. Defaults to all keys.
        **kwargs: Datasource-specific knobs (``subset`` for ChEBI, ``species``
            for Ensembl); see :func:`get_download_urls`.

    Returns:
        Dictionary mapping file keys to downloaded paths.
    """
    return _download_datasource(
        datasource,
        output_dir,
        all_datasources=ALL_DATASOURCES,
        urls_and_date=URLS_AND_DATE,
        decompress=decompress,
        version=version,
        keys=keys,
        tar_extractors=TAR_EXTRACTORS,
        **kwargs,
    )


def download_datasource_with_release(
    datasource: str,
    output_dir: Path,
    decompress: bool = True,
    version: str | None = None,
    keys: list[str] | None = None,
    **kwargs: Any,
) -> tuple[dict[str, Path], datetime | None]:
    """Download all files for a datasource and report its release date.

    Same as :func:`download_datasource`, but also returns the resolved source
    release date (used for the SSSOM ``mapping_date``).

    Args:
        datasource: Name of the datasource.
        output_dir: Directory to save files.
        decompress: Whether to decompress .gz files.
        version: Specific version to download. Format depends on datasource.
        keys: Optional list of file-key names to download.
        **kwargs: Datasource-specific knobs (``subset`` for ChEBI, ``species``
            for Ensembl); see :func:`get_download_urls`.

    Returns:
        Tuple of (file-key -> downloaded path mapping, release date or None).
    """
    return _download_datasource_with_release(
        datasource,
        output_dir,
        all_datasources=ALL_DATASOURCES,
        urls_and_date=URLS_AND_DATE,
        decompress=decompress,
        version=version,
        keys=keys,
        tar_extractors=TAR_EXTRACTORS,
        **kwargs,
    )
