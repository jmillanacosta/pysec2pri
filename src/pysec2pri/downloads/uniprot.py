"""Downloading and release detection for UniProt."""

from __future__ import annotations

import re
import tarfile
from datetime import datetime
from typing import TYPE_CHECKING, Any

import httpx
from mapkgsutils.download import ReleaseInfo

from pysec2pri.logging import logger
from pysec2pri.parsers.base import BaseDownloader

if TYPE_CHECKING:
    from pathlib import Path

    from pysec2pri.parsers.base import DatasourceConfig

__all__ = ["UniProtDownloader", "check_uniprot_release", "extract_tar", "urls_and_date"]

METALINK_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/RELEASE.metalink"
_NS = {"m": "http://www.metalinker.org/"}


class UniProtDownloader(BaseDownloader):
    """Downloader for UniProt data files."""

    datasource_name = "uniprot"

    def get_download_urls(
        self,
        version: str | None = None,
        **kwargs: object,
    ) -> dict[str, str]:
        """Get UniProt download URLs for *version*, or latest."""
        if version:
            return _get_uniprot_urls_for_version(version)
        if self._config:
            return dict(self._config.download_urls)
        raise ValueError("UniProt config not loaded")

    def download(
        self,
        output_dir: Path,
        version: str | None = None,
        decompress: bool = True,
        **kwargs: object,
    ) -> dict[str, Path]:
        """Download UniProt files into *output_dir*."""
        urls = self.get_download_urls(version)
        return self._download_urls(urls, output_dir, decompress)

    def list_versions(self) -> list[str]:
        """List all available UniProt previous-release versions.

        Scrapes the UniProt FTP previous_releases directory for version
        strings.

        Returns:
            Sorted list of version strings
            (e.g. ``["2024_01", "2024_02", ...]``).

        Raises:
            ValueError: If the archive URL is not configured.
        """
        if not self._config or not self._config.archive_url:
            raise ValueError("UniProt archive URL not configured")
        with httpx.Client(follow_redirects=True, timeout=30.0) as client:
            response = client.get(self._config.archive_url)
            response.raise_for_status()
        # FTP HTML index: links like "release-2024_01/"
        matches = re.findall(r'href="release-(\d{4}_\d{2})/', response.text)
        return sorted(set(matches))


def _get_uniprot_urls_for_version(version: str) -> dict[str, str]:
    """Build UniProt URLs for a specific release version.

    Args:
        version: Release version (e.g., "2024_01").

    Returns:
        Dictionary with 'knowledgebase_docs' URL (tar.gz containing sec_ac.txt).
    """
    base = "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases"
    return {
        "knowledgebase_docs": (
            f"{base}/release-{version}/knowledgebase/knowledgebase-docs-only{version}.tar.gz"
        ),
    }


def _fetch_uniprot_metalink(url: str = METALINK_URL) -> str:
    import requests

    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.text


def _parse_uniprot_release_xml(xml_text: str) -> tuple[str, datetime | None]:
    from defusedxml.ElementTree import fromstring

    root = fromstring(xml_text)

    version_node = root.find("m:version", _NS)
    if version_node is None or not version_node.text:
        raise ValueError("Missing UniProt version in metalink")

    version = version_node.text.strip()
    return version, None


def check_uniprot_release() -> ReleaseInfo:
    """Check for the latest UniProt release.

    Returns:
        ReleaseInfo with the latest UniProt release details.
    """
    from pysec2pri.constants import ALL_DATASOURCES

    xml_text = _fetch_uniprot_metalink()
    version, release_date = _parse_uniprot_release_xml(xml_text)

    return ReleaseInfo(
        datasource="uniprot",
        version=version,  # e.g. "2026_01"
        release_date=release_date,  # likely None unless you derive it elsewhere
        is_new=True,
        files=dict(ALL_DATASOURCES["uniprot"].download_urls),
    )


def urls_and_date(
    version: str | None, config: DatasourceConfig, **kwargs: Any
) -> tuple[dict[str, str], datetime | None]:
    """Resolve UniProt download URLs and the source release date.

    Args:
        version: Specific release version (e.g. ``"2024_01"``), or ``None``
            for latest.
        config: UniProt's datasource configuration.
        **kwargs: Unused; accepted for interface uniformity.

    Returns:
        Tuple of (file-key -> URL mapping, release date or None).
    """
    if version:
        urls = _get_uniprot_urls_for_version(version)
        logger.info("UniProt version %s: %s", version, urls)
        return urls, None
    return dict(config.download_urls), None


def extract_tar(tar_path: Path, output_dir: Path) -> dict[str, Path]:
    """Extract UniProt tar.gz and return paths to sec_ac.txt and delac_sp.txt.

    Args:
        tar_path: Path to the downloaded tar.gz file.
        output_dir: Directory to extract to.

    Returns:
        Dictionary mapping file keys to extracted paths.
    """
    extracted = {}
    with tarfile.open(tar_path, "r:gz") as tar:
        for member in tar.getmembers():
            if member.name.endswith("sec_ac.txt"):
                tar.extract(member, output_dir)
                extracted["sec_ac"] = output_dir / member.name
            elif member.name.endswith("delac_sp.txt"):
                tar.extract(member, output_dir)
                extracted["delac_sp"] = output_dir / member.name
    return extracted
