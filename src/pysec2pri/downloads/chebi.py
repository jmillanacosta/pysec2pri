"""Downloading and release detection for ChEBI."""

from __future__ import annotations

import re
from datetime import datetime
from typing import TYPE_CHECKING, Any

import httpx
from mapkgsutils.download import ReleaseInfo, get_file_last_modified

from pysec2pri.constants import ALL_DATASOURCES
from pysec2pri.logging import logger
from pysec2pri.parsers.base import BaseDownloader

if TYPE_CHECKING:
    from pathlib import Path

    from pysec2pri.parsers.base import DatasourceConfig

__all__ = ["ChEBIDownloader", "check_chebi_release", "urls_and_date"]


class ChEBIDownloader(BaseDownloader):
    """Downloader for ChEBI data files.

    Supports both TSV flat files (>= release 245) and legacy SDF files.
    """

    datasource_name = "chebi"

    def __init__(
        self,
        version: str | None = None,
        show_progress: bool = True,
        subset: str = "3star",
        use_sdf: bool = False,
    ) -> None:
        """Initialize the ChEBI downloader.

        Args:
            version: Version/release identifier.
            show_progress: Whether to show progress bars.
            subset: "3star" or "complete" - which compound subset to use.
            use_sdf: Force SDF format even for releases >= 245.
        """
        super().__init__(version=version, show_progress=show_progress)
        self.subset = subset
        self.use_sdf = use_sdf

    def get_download_urls(
        self,
        version: str | None = None,
        **kwargs: Any,
    ) -> dict[str, str]:
        """Get ChEBI download URLs based on version.

        Resolves the distribution era for *version* via
        :meth:`~mapkgsutils.parsers.base.DatasourceConfig.era_for` (see
        ``chebi.yaml``'s ``distribution_eras``): the ``sdf`` era (releases
        <= 244) returns a single SDF URL, the ``tsv`` era (releases >= 245)
        returns the three flat-file URLs. When no era matches (e.g.
        ``version`` is ``None`` for "latest"), falls back to the top-level
        ``download_urls`` (TSV).

        Args:
            version: Specific release version.
            **kwargs: Optional 'subset' and 'use_sdf' overrides.

        Returns:
            Dictionary with file URLs keyed by type.
        """
        v = version or self.version
        subset = kwargs.get("subset", self.subset)
        use_sdf = kwargs.get("use_sdf", self.use_sdf)

        if not self._config:
            raise ValueError("ChEBI config not loaded")

        era = self._config.era_for(v)
        if use_sdf and (era is None or era.format != "sdf"):
            era = next((e for e in self._config.distribution_eras if e.format == "sdf"), era)

        if era is not None and era.format == "sdf":
            sdf_key = "sdf_3star" if subset == "3star" else "sdf_complete"
            url = era.download_urls[sdf_key].format(version=v)
            return {"sdf": url}

        download_urls = era.download_urls if era is not None else self._config.download_urls
        return {
            "secondary_ids": download_urls["secondary_ids"].format(version=v),
            "names": download_urls["names"].format(version=v),
            "compounds": download_urls["compounds"].format(version=v),
        }

    def download(
        self,
        output_dir: Path,
        version: str | None = None,
        decompress: bool = True,
        **kwargs: Any,
    ) -> dict[str, Path]:
        """Download ChEBI files.

        Args:
            output_dir: Directory to save files.
            version: Specific version to download.
            decompress: Whether to decompress .gz files.
            **kwargs: Optional 'subset' and 'use_sdf' overrides.

        Returns:
            Dictionary mapping file keys to downloaded paths.
        """
        v = version or self.version
        urls = self.get_download_urls(v, **kwargs)

        return self._download_urls(urls, output_dir, decompress)

    def get_format(self, version: str | None = None) -> str:
        """Get the format that will be used for a given version.

        Args:
            version: Version to check. If None, uses self.version.

        Returns:
            "tsv" or "sdf"
        """
        v = version or self.version
        if self.use_sdf:
            return "sdf"
        if self._config:
            era = self._config.era_for(v)
            if era is not None and era.format:
                return era.format
        return "tsv" if self.is_new_format(v) else "sdf"

    def list_versions(self) -> list[str]:
        """List all available ChEBI archive release numbers.

        Queries the base archive URL plus every distribution era's
        ``archive_url`` (e.g. the legacy SDF-only archive for releases <=
        244, see ``chebi.yaml``'s ``distribution_eras``) and unions all
        discovered release numbers.

        Returns:
            Sorted list of version strings (e.g. ``["100", "101", ..., "251"]``).

        Raises:
            ValueError: If the archive URL is not configured.
        """
        if not self._config or not self._config.archive_url:
            raise ValueError("ChEBI archive URL not configured")

        archive_urls = {self._config.archive_url}
        archive_urls.update(
            era.archive_url for era in self._config.distribution_eras if era.archive_url
        )

        versions: set[str] = set()
        with httpx.Client(follow_redirects=True, timeout=30.0) as client:
            for url in archive_urls:
                try:
                    response = client.get(url)
                    response.raise_for_status()
                    versions.update(re.findall(r"rel(\d+)", response.text))
                except httpx.HTTPError:
                    pass  # That archive location unavailable, still return what we have

        return sorted(versions, key=int)


def _get_chebi_urls_for_version(version: str, subset: str = "3star") -> dict[str, str]:
    """Build ChEBI URLs for a specific release version.

    Delegates to :meth:`ChEBIDownloader.get_download_urls`, which resolves
    the distribution era (TSV >= 245, legacy SDF <= 244) from ``chebi.yaml``'s
    ``distribution_eras``.

    Args:
        version: Release number (e.g., "232", "245").
        subset: Either "3star" or "complete". For new releases (>=245),
            this determines which compounds are included via the compounds.tsv
            filter. For legacy releases (<245), this determines the SDF file.

    Returns:
        Dictionary with file URLs keyed by type (sdf, secondary_ids, etc.).
    """
    return ChEBIDownloader(version=version, subset=subset).get_download_urls(version)


def check_chebi_release() -> ReleaseInfo:
    """Check for the latest ChEBI release.

    Returns:
        ReleaseInfo with the latest ChEBI release details.
    """
    archive_url = ALL_DATASOURCES["chebi"].archive_url
    with httpx.Client(follow_redirects=True) as client:
        response = client.get(archive_url)
        response.raise_for_status()

    # Parse the archive index to find the latest release number
    matches = re.findall(r'href="rel(\d+)/"', response.text)
    if not matches:
        raise ValueError("Could not find ChEBI releases in archive")

    latest_release = max(int(m) for m in matches)
    version = str(latest_release)

    # Return URLs based on release version
    urls = _get_chebi_urls_for_version(version, subset="3star")
    # Get release date from secondary_ids (new) or sdf (legacy)
    check_url = urls.get("secondary_ids") or urls.get("sdf")
    release_date = get_file_last_modified(check_url) if check_url else None

    return ReleaseInfo(
        datasource="chebi",
        version=version,
        release_date=release_date,
        is_new=True,  # Caller determines if it's new
        files=urls,
    )


def urls_and_date(
    version: str | None, config: DatasourceConfig, **kwargs: Any
) -> tuple[dict[str, str], datetime | None]:
    """Resolve ChEBI download URLs and the source release date.

    Args:
        version: Specific release version, or ``None`` for latest.
        config: ChEBI's datasource configuration.
        **kwargs: Optional ``subset`` override.

    Returns:
        Tuple of (file-key -> URL mapping, release date or None).
    """
    subset = kwargs.get("subset") or config.default_subset() or "3star"
    release_date: datetime | None = None
    if version:
        urls = _get_chebi_urls_for_version(version, subset=subset)
        logger.info("ChEBI version %s: %s", version, urls)
    else:
        release_info = check_chebi_release()
        urls = release_info.files
        release_date = release_info.release_date
        logger.info("ChEBI release %s: %s", release_info.version, urls)
    return urls, release_date
