"""Copy to ``src/pysec2pri/downloads/<source>.py``.

Only needed when the URL needs a version to have to find the (latest) release.
Otherwise the config's ``download_urls`` provide the whole download logic and
this file is not needed.

Example keeps every release at ``https://example.org/{version}/``, and lists
them at ``https://example.org/archive/``.
"""

from __future__ import annotations

import re
from datetime import datetime
from typing import TYPE_CHECKING, Any

import httpx
from mapkgsutils.download import ReleaseInfo

from pysec2pri.logging import logger
from pysec2pri.parsers.base import BaseDownloader

if TYPE_CHECKING:
    from pathlib import Path

    from pysec2pri.parsers.base import DatasourceConfig

__all__ = ["ExampleDownloader", "check_example_release", "urls_and_date"]

_ARCHIVE_URL = "https://example.org/archive/"


def _urls_for_version(version: str) -> dict[str, str]:
    """Return each file's download_urls key -> URL, for one release."""
    return {
        "withdrawn": f"https://example.org/{version}/withdrawn.tsv",
        "complete": f"https://example.org/{version}/complete.tsv",
    }


class ExampleDownloader(BaseDownloader):
    """Downloader for Example's per-release files."""

    datasource_name = "example"

    def get_download_urls(self, version: str | None = None, **kwargs: Any) -> dict[str, str]:
        """Return each file's download_urls key -> URL.

        Args:
            version: Release to build URLs for. Latest when ``None``.
            **kwargs: Unused; accepted for interface uniformity.
        """
        return _urls_for_version(version or check_example_release().version or "")

    def download(
        self,
        output_dir: Path,
        version: str | None = None,
        decompress: bool = True,
        **kwargs: Any,
    ) -> dict[str, Path]:
        """Download the files into *output_dir*, returning key -> local path."""
        return self._download_urls(self.get_download_urls(version), output_dir, decompress)

    def list_versions(self) -> list[str]:
        """Return every release Example keeps, oldest first.

        Only needed for a source with an archive: this is what
        ``pysec2pri list-versions`` and ``--consolidate`` walk.
        """
        with httpx.Client(follow_redirects=True, timeout=30.0) as client:
            response = client.get(_ARCHIVE_URL)
            response.raise_for_status()
        return sorted(set(re.findall(r'href="(\d{4}-\d{2}-\d{2})/"', response.text)))


def check_example_release() -> ReleaseInfo:
    """Return the newest release: its version, its date, and its files."""
    with httpx.Client(follow_redirects=True, timeout=30.0) as client:
        response = client.get(_ARCHIVE_URL)
        response.raise_for_status()

    versions = re.findall(r'href="(\d{4}-\d{2}-\d{2})/"', response.text)
    if not versions:
        raise ValueError(f"No Example releases found at {_ARCHIVE_URL}")

    latest = max(versions)
    logger.info("Latest Example release: %s", latest)
    return ReleaseInfo(
        datasource="example",
        version=latest,
        release_date=datetime.strptime(latest, "%Y-%m-%d"),
        is_new=True,
        files=_urls_for_version(latest),
    )


def urls_and_date(
    version: str | None, config: DatasourceConfig, **kwargs: Any
) -> tuple[dict[str, str], datetime | None]:
    """Return the files to download and the date that release came out.

    The date becomes each mapping's ``mapping_date``.

    Args:
        version: Release to download, or ``None`` for the latest.
        config: Example's config. Unused here: the URLs come from the version.
        **kwargs: Unused in example, use if needed.
    """
    if version is None:
        release = check_example_release()
        return release.files, release.release_date
    return _urls_for_version(version), datetime.strptime(version, "%Y-%m-%d")
