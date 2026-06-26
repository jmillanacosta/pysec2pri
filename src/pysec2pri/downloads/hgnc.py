"""Downloading and release detection for HGNC."""

from __future__ import annotations

import re
from datetime import datetime
from typing import TYPE_CHECKING, Any

import httpx
from mapkgsutils.download import ReleaseInfo, get_file_last_modified

from pysec2pri.logging import logger
from pysec2pri.parsers.base import BaseDownloader

if TYPE_CHECKING:
    from pathlib import Path

    from pysec2pri.parsers.base import DatasourceConfig

__all__ = ["HGNCDownloader", "check_hgnc_release", "urls_and_date"]

_GCS_API_URL = (
    "https://storage.googleapis.com/storage/v1/b/public-download-files/o"
    "?prefix=hgnc/archive/archive/quarterly/tsv/"
)
_ARCHIVE_BASE_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv"
)


class HGNCDownloader(BaseDownloader):
    """Downloader for HGNC data files from the quarterly archive."""

    datasource_name = "hgnc"

    def get_download_urls(
        self,
        version: str | None = None,
        **kwargs: object,
    ) -> dict[str, str]:
        """Get HGNC download URLs for *version* (``YYYY-MM-DD``), or latest."""
        if version:
            return _get_hgnc_urls_for_version(version)
        return check_hgnc_release().files

    def download(
        self,
        output_dir: Path,
        version: str | None = None,
        decompress: bool = True,
        **kwargs: object,
    ) -> dict[str, Path]:
        """Download HGNC files into *output_dir*."""
        urls = self.get_download_urls(version)
        return self._download_urls(urls, output_dir, decompress)

    def list_versions(self) -> list[str]:
        """List all available HGNC quarterly archive versions.

        Queries the Google Cloud Storage API for all complete-set files and
        returns their dates in ascending order.

        Returns:
            Sorted list of version strings in ``YYYY-MM-DD`` format.
        """
        with httpx.Client(follow_redirects=True, timeout=30.0) as client:
            response = client.get(_GCS_API_URL)
            response.raise_for_status()
            data = response.json()

        items = data.get("items", [])
        versions: list[str] = []
        for item in items:
            name = item.get("name", "")
            match = re.search(r"hgnc_complete_set_(\d{4}-\d{2}-\d{2})\.txt$", name)
            if match:
                if f"withdrawn_{match.group(1)}" not in response.text:
                    continue
                versions.append(match.group(1))
        return sorted(set(versions))


def _get_hgnc_urls_for_version(version: str) -> dict[str, str]:
    """Build HGNC URLs for a specific version date.

    Args:
        version: Version string in YYYY-MM-DD format.

    Returns:
        Dictionary with 'withdrawn' and 'complete' URLs.
    """
    return {
        "withdrawn": f"{_ARCHIVE_BASE_URL}/withdrawn_{version}.txt",
        "complete": f"{_ARCHIVE_BASE_URL}/hgnc_complete_set_{version}.txt",
    }


def check_hgnc_release() -> ReleaseInfo:
    """Check for the latest HGNC release from the quarterly archive.

    Queries the Google Cloud Storage API to list files in the HGNC
    quarterly archive bucket and finds the latest release.

    Returns:
        ReleaseInfo with the latest HGNC release details.
    """
    logger.info("Querying HGNC files from GCS API: %s", _GCS_API_URL)

    with httpx.Client(follow_redirects=True, timeout=30.0) as client:
        response = client.get(_GCS_API_URL)
        response.raise_for_status()
        data = response.json()

    # Extract file names from the response
    items = data.get("items", [])
    file_names = [item.get("name", "") for item in items]

    # Find all complete set files: hgnc_complete_set_YYYY-MM-DD.txt
    complete_dates = []
    for name in file_names:
        match = re.search(r"hgnc_complete_set_(\d{4}-\d{2}-\d{2})\.txt$", name)
        if match:
            complete_dates.append(match.group(1))

    if not complete_dates:
        logger.warning("Could not find HGNC complete set files in GCS bucket")
        # Fall back to current release URLs
        from pysec2pri.constants import ALL_DATASOURCES

        last_modified = get_file_last_modified(ALL_DATASOURCES["hgnc"].download_urls["complete"])
        version = last_modified.strftime("%Y-%m-%d") if last_modified else None
        return ReleaseInfo(
            datasource="hgnc",
            version=version,
            release_date=last_modified,
            is_new=True,
            files=dict(ALL_DATASOURCES["hgnc"].download_urls),
        )

    # Get the latest date
    latest_date = max(complete_dates)
    logger.info("Latest HGNC release date: %s", latest_date)

    # Parse the date
    release_date = datetime.strptime(latest_date, "%Y-%m-%d")

    return ReleaseInfo(
        datasource="hgnc",
        version=latest_date,
        release_date=release_date,
        is_new=True,
        files=_get_hgnc_urls_for_version(latest_date),
    )


def _iso_to_datetime(value: str | None) -> datetime | None:
    """Parse an ``YYYY-MM-DD`` string into a datetime, or return None."""
    try:
        return datetime.strptime(str(value), "%Y-%m-%d")
    except (ValueError, TypeError):
        return None


def urls_and_date(
    version: str | None, config: DatasourceConfig, **kwargs: Any
) -> tuple[dict[str, str], datetime | None]:
    """Resolve HGNC download URLs and the source release date.

    Args:
        version: Specific archive date (``YYYY-MM-DD``), or ``None`` for latest.
        config: HGNC's datasource configuration (unused; HGNC's URLs are
            always derived from the version/archive directly).
        **kwargs: Unused; accepted for interface uniformity.

    Returns:
        Tuple of (file-key -> URL mapping, release date or None).
    """
    if version:
        urls = _get_hgnc_urls_for_version(version)
        # The HGNC archive version is itself the release date (YYYY-MM-DD).
        release_date = _iso_to_datetime(version)
        logger.info("HGNC version %s: %s", version, urls)
    else:
        release_info = check_hgnc_release()
        urls = release_info.files
        release_date = release_info.release_date
        logger.info("HGNC release %s: %s", release_info.version, urls)
    return urls, release_date
