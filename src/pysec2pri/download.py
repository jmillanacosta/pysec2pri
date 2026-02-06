"""Download and release detection for biological database sources."""

from __future__ import annotations

import gzip
import re
from collections.abc import Generator, Iterable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import httpx
from tqdm import tqdm

from pysec2pri.constants import ALL_DATASOURCES
from pysec2pri.logging import logger
from pysec2pri.parsers.base import DatasourceConfig

__all__ = [
    "ReleaseInfo",
    "check_release",
    "download_datasource",
    "download_file",
    "get_download_urls",
    "get_latest_release_info",
]


@dataclass
class ReleaseInfo:
    """Information about a datasource release."""

    datasource: str
    version: str | None
    release_date: datetime | None
    is_new: bool
    files: dict[str, str]  # key -> URL mapping


def download_file(
    url: str,
    output_path: Path,
    decompress_gz: bool = True,
    timeout: float | None = None,
    show_progress: bool = True,
    description: str | None = None,
) -> Path:
    """Download a file from URL to the specified path.

    Args:
        url: URL to download from.
        output_path: Where to save the file.
        decompress_gz: Whether to decompress .gz files automatically.
        timeout: Request timeout in seconds.
        show_progress: Whether to show a progress bar.
        description: Description for the progress bar.

    Returns:
        Path to the downloaded (and optionally decompressed) file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Get filename for progress bar description
    if description is None:
        description = output_path.name
    with httpx.Client(timeout=timeout, follow_redirects=True) as client:
        with client.stream("GET", url) as response:
            response.raise_for_status()

            # Get total size if available
            total_size = int(response.headers.get("content-length", 0))

            # Determine if we need to decompress
            is_gzip = url.endswith(".gz") and decompress_gz
            final_path = output_path

            if is_gzip:
                _download_gzip(
                    output_path,
                    show_progress,
                    response,
                    total_size,
                    final_path,
                    description,
                )
            else:
                _download_nogzip(
                    output_path,
                    total_size,
                    response,
                    show_progress,
                    description,
                )
    return final_path


def get_file_last_modified(url: str, timeout: float = 30.0) -> datetime | None:
    """Get the Last-Modified date from a URL via HEAD request.

    Args:
        url: URL to check.
        timeout: Request timeout in seconds.

    Returns:
        The Last-Modified datetime or None if unavailable.
    """
    try:
        with httpx.Client(timeout=timeout, follow_redirects=True) as client:
            response = client.head(url)
            if "last-modified" in response.headers:
                from email.utils import parsedate_to_datetime

                return parsedate_to_datetime(response.headers["last-modified"])
    except (httpx.HTTPError, ValueError):
        pass
    return None


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


def check_uniprot_release() -> ReleaseInfo:
    """Check for the latest UniProt release.

    Returns:
        ReleaseInfo with the latest UniProt release details.
    """
    sec_ac_url = ALL_DATASOURCES["uniprot"].download_urls["secondary"]
    last_modified = get_file_last_modified(sec_ac_url)
    version = last_modified.strftime("%Y-%m-%d") if last_modified else None

    return ReleaseInfo(
        datasource="uniprot",
        version=version,
        release_date=last_modified,
        is_new=True,
        files=dict(ALL_DATASOURCES["uniprot"].download_urls),
    )


def check_hmdb_release() -> ReleaseInfo:
    """Check for the latest HMDB release by downloading and checking XML.

    Returns:
        ReleaseInfo with the latest HMDB release details.
    """
    # HMDB requires downloading to check the version in XML
    metabolites_url = ALL_DATASOURCES["hmdb"].download_urls["metabolites"]
    last_modified = get_file_last_modified(metabolites_url)
    version = last_modified.strftime("%Y-%m-%d") if last_modified else None

    return ReleaseInfo(
        datasource="hmdb",
        version=version,
        release_date=last_modified,
        is_new=True,
        files=dict(ALL_DATASOURCES["hmdb"].download_urls),
    )


def check_hgnc_release() -> ReleaseInfo:
    """Check for the latest HGNC release from the quarterly archive.

    Queries the Google Cloud Storage API to list files in the HGNC
    quarterly archive bucket and finds the latest release.

    Returns:
        ReleaseInfo with the latest HGNC release details.
    """
    gcs_api_url = (
        "https://storage.googleapis.com/storage/v1/b/public-download-files/o"
        "?prefix=hgnc/archive/archive/quarterly/tsv/"
    )
    logger.info(f"Querying HGNC files from GCS API: {gcs_api_url}")

    with httpx.Client(follow_redirects=True, timeout=30.0) as client:
        response = client.get(gcs_api_url)
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
    logger.info(f"Latest HGNC release date: {latest_date}")

    # Build URLs for the quarterly archive files
    # Files are hosted on Google Cloud Storage
    base_url = (
        "https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv"
    )
    complete_url = f"{base_url}/hgnc_complete_set_{latest_date}.txt"
    withdrawn_url = f"{base_url}/withdrawn_{latest_date}.txt"

    # Parse the date
    release_date = datetime.strptime(latest_date, "%Y-%m-%d")

    return ReleaseInfo(
        datasource="hgnc",
        version=latest_date,
        release_date=release_date,
        is_new=True,
        files={
            "complete": complete_url,
            "withdrawn": withdrawn_url,
        },
    )


def get_latest_release_info(datasource: str) -> ReleaseInfo:
    """Get release information for a datasource.

    Args:
        datasource: Name of the datasource (chebi, hmdb, hgnc, ncbi, uniprot).

    Returns:
        ReleaseInfo with the latest release details.

    Raises:
        ValueError: If the datasource is not supported.
    """
    checkers = {
        "chebi": check_chebi_release,
        "hmdb": check_hmdb_release,
        "hgnc": check_hgnc_release,
        "ncbi": check_ncbi_release,
        "uniprot": check_uniprot_release,
    }

    if datasource.lower() not in checkers:
        raise ValueError(f"Unknown datasource: {datasource}. Supported: {list(checkers.keys())}")

    return checkers[datasource.lower()]()


def _get_hgnc_urls_for_version(version: str) -> dict[str, str]:
    """Build HGNC URLs for a specific version date.

    Args:
        version: Version string in YYYY-MM-DD format.

    Returns:
        Dictionary with 'withdrawn' and 'complete' URLs.
    """
    base_url = (
        "https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv"
    )
    return {
        "withdrawn": f"{base_url}/withdrawn_{version}.txt",
        "complete": f"{base_url}/hgnc_complete_set_{version}.txt",
    }


def _get_chebi_urls_for_version(
    version: str,
    subset: str = "3star",
) -> dict[str, str]:
    """Build ChEBI URLs for a specific release version.

    For releases >= 245, returns TSV flat file URLs.
    For releases < 245, returns SDF file URL.

    URLs are loaded from chebi.yaml config file.

    Args:
        version: Release number (e.g., "232", "245").
        subset: Either "3star" or "complete". For new releases (>=245),
            this determines which compounds are included via the compounds.tsv
            filter. For legacy releases (<245), this determines the SDF file.

    Returns:
        Dictionary with file URLs keyed by type (sdf, secondary_ids, etc.).
    """
    config = ALL_DATASOURCES.get("chebi")
    if not config:
        raise ValueError("ChEBI config not loaded")

    new_format_version = config.new_format_version or 245
    download_urls = config.download_urls
    release_num = int(version)

    if release_num >= new_format_version:
        # New TSV format (>= 245)
        new_urls = download_urls.get("new", {})
        return {
            "secondary_ids": new_urls["secondary_ids"].format(version=version),
            "names": new_urls["names"].format(version=version),
            "compounds": new_urls["compounds"].format(version=version),
        }
    else:
        # Legacy SDF format (< 245)
        legacy_urls = download_urls.get("legacy", {})
        sdf_key = "sdf_3star" if subset == "3star" else "sdf_complete"
        url = legacy_urls[sdf_key].format(version=version)
        return {"sdf": url}


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


def get_download_urls(
    datasource: str,
    version: str | None = None,
    subset: str = "3star",
) -> dict[str, str]:
    """Get download URLs for a datasource.

    Args:
        datasource: Name of the datasource.
        version: Specific version to get URLs for.
        subset: For ChEBI: "3star" or "complete" (only affects legacy SDF).

    Returns:
        Dictionary mapping file keys to URLs.
    """
    datasource_lower = datasource.lower()
    config = ALL_DATASOURCES.get(datasource_lower)
    if not config:
        raise ValueError(f"Unknown datasource: {datasource}")

    if datasource_lower == "hgnc":
        if version:
            return _get_hgnc_urls_for_version(version)
        return check_hgnc_release().files
    elif datasource_lower == "chebi":
        if version:
            return _get_chebi_urls_for_version(version, subset=subset)
        return check_chebi_release().files
    elif datasource_lower == "uniprot":
        if version:
            return _get_uniprot_urls_for_version(version)
        return dict(config.download_urls)
    else:
        return dict(config.download_urls)


def _get_datasource_urls(
    datasource_lower: str,
    config: DatasourceConfig,
    version: str | None = None,
    subset: str = "3star",
) -> dict[str, str]:
    """Get download URLs for a datasource.

    Args:
        datasource_lower: Lowercase datasource name.
        config: Datasource configuration.
        version: Specific version to get URLs for.
        subset: For ChEBI: "3star" or "complete" (only affects legacy SDF).

    Returns:
        Dictionary mapping file keys to URLs.
    """
    if datasource_lower == "hgnc":
        if version:
            urls = _get_hgnc_urls_for_version(version)
            logger.info(f"HGNC version {version}: {urls}")
        else:
            release_info = check_hgnc_release()
            urls = release_info.files
            logger.info(f"HGNC release {release_info.version}: {urls}")
        return urls

    if datasource_lower == "chebi":
        if version:
            urls = _get_chebi_urls_for_version(version, subset=subset)
            logger.info(f"ChEBI version {version}: {urls}")
        else:
            release_info = check_chebi_release()
            urls = release_info.files
            logger.info(f"ChEBI release {release_info.version}: {urls}")
        return urls

    if datasource_lower == "uniprot":
        if version:
            urls = _get_uniprot_urls_for_version(version)
            logger.info(f"UniProt version {version}: {urls}")
            return urls
        return dict(config.download_urls)

    # NCBI, HMDB - no versioned archives available
    if version:
        logger.warning(
            f"{datasource_lower} does not have versioned archives. "
            f"Downloading latest version instead."
        )
    return dict(config.download_urls)


def download_datasource(
    datasource: str,
    output_dir: Path,
    decompress: bool = True,
    version: str | None = None,
    subset: str = "3star",
) -> dict[str, Path]:
    """Download all files for a datasource.

    For datasources with dynamic URLs (like HGNC quarterly archive),
    this function first checks for the latest release and uses those URLs.
    If a version is specified, it downloads that specific version.

    Args:
        datasource: Name of the datasource.
        output_dir: Directory to save files.
        decompress: Whether to decompress .gz files.
        version: Specific version to download. Format depends on datasource:
            - HGNC: date string "YYYY-MM-DD" (e.g., "2026-01-06")
            - ChEBI: release number (e.g., "232")
        subset: For ChEBI: "3star" or "complete" (only affects legacy SDF).

    Returns:
        Dictionary mapping file keys to downloaded paths.
    """
    datasource_lower = datasource.lower()
    config = ALL_DATASOURCES.get(datasource_lower)
    if not config:
        raise ValueError(f"Unknown datasource: {datasource}")

    output_dir.mkdir(parents=True, exist_ok=True)
    urls = _get_datasource_urls(datasource_lower, config, version, subset)

    return _download_urls(urls, output_dir, decompress)


def _download_urls(
    urls: dict[str, str],
    output_dir: Path,
    decompress: bool = True,
) -> dict[str, Path]:
    """Download files from URLs to output directory.

    Args:
        urls: Dictionary mapping file keys to URLs.
        output_dir: Directory to save files.
        decompress: Whether to decompress .gz files.

    Returns:
        Dictionary mapping file keys to downloaded paths.
    """
    downloaded = {}

    for key, url in urls.items():
        filename = url.split("/")[-1]

        # Handle tar.gz files (UniProt archive)
        if filename.endswith(".tar.gz"):
            output_path = output_dir / filename
            logger.info(f"Downloading {key}: {url}")
            download_file(url, output_path, decompress_gz=False)
            extracted = _extract_uniprot_tar(output_path, output_dir)
            downloaded.update(extracted)
            logger.info(f"Extracted: {list(extracted.keys())}")
            continue

        if decompress and filename.endswith(".gz"):
            filename = filename[:-3]

        output_path = output_dir / filename
        logger.info(f"Downloading {key}: {url}")
        download_file(url, output_path, decompress_gz=decompress)
        downloaded[key] = output_path
        logger.info(f"Saved to: {output_path}")

    return downloaded


def _extract_uniprot_tar(tar_path: Path, output_dir: Path) -> dict[str, Path]:
    """Extract UniProt tar.gz and return paths to sec_ac.txt and delac_sp.txt.

    Args:
        tar_path: Path to the downloaded tar.gz file.
        output_dir: Directory to extract to.

    Returns:
        Dictionary mapping file keys to extracted paths.
    """
    import tarfile

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
    info = get_latest_release_info(datasource)

    # Determine if this is a new release
    if current_version and info.version:
        info = ReleaseInfo(
            datasource=info.datasource,
            version=info.version,
            release_date=info.release_date,
            is_new=info.version != current_version,
            files=info.files,
        )
    elif current_date and info.release_date:
        info = ReleaseInfo(
            datasource=info.datasource,
            version=info.version,
            release_date=info.release_date,
            is_new=info.release_date > current_date,
            files=info.files,
        )

    return info


@contextmanager
def iter_with_progress(
    iterator: Iterable[bytes],
    *,
    enabled: bool,
    total: int | None,
    description: str,
) -> Iterator[Iterable[bytes]]:
    """Wrap an iterator with optional progress bar.

    Args:
        iterator: The byte iterator to wrap.
        enabled: Whether to show progress.
        total: Total size in bytes.
        description: Description for progress bar.

    Yields:
        The wrapped iterator.
    """
    if enabled and total and total > 0:
        with tqdm(
            total=total,
            unit="B",
            unit_scale=True,
            desc=description,
        ) as pbar:

            def gen() -> Generator[bytes, None, None]:
                """Yield chunks."""
                for chunk in iterator:
                    pbar.update(len(chunk))
                    yield chunk

            yield gen()
    else:
        yield iterator


def _download_gzip(
    output_path: Path,
    show_progress: bool,
    response: httpx.Response,
    total_size: int,
    final_path: Path,
    description: str | None = None,
) -> None:
    """Help download gzipped files.

    Args:
        output_path: Path for the output file.
        show_progress: Whether to show progress bar.
        response: HTTP response object.
        total_size: Total size in bytes.
        final_path: Final destination path.
        description: Description for progress bar.
    """
    temp_path = output_path.with_suffix(output_path.suffix + ".gz")

    with temp_path.open("wb") as f:
        with iter_with_progress(
            response.iter_bytes(chunk_size=8192),
            enabled=show_progress,
            total=total_size,
            description=f"Downloading {description}",
        ) as chunks:
            for chunk in chunks:
                f.write(chunk)

    if show_progress:
        compressed_size = temp_path.stat().st_size
    else:
        compressed_size = None

    with gzip.open(temp_path, "rb") as f_in, final_path.open("wb") as f_out:
        with iter_with_progress(
            iter(lambda: f_in.read(8192), b""),
            enabled=show_progress,
            total=compressed_size,
            description=f"Decompressing {description}",
        ) as chunks:
            for chunk in chunks:
                f_out.write(chunk)

    temp_path.unlink()


def _download_nogzip(
    output_path: Path,
    total_size: int,
    response: httpx.Response,
    show_progress: bool,
    description: str | None,
) -> None:
    """Help download non-gzipped files.

    Args:
        output_path: Path for the output file.
        total_size: Total size in bytes.
        response: HTTP response object.
        show_progress: Whether to show progress bar.
        description: Description for progress bar.
    """
    with output_path.open("wb") as f:
        with iter_with_progress(
            response.iter_bytes(chunk_size=8192),
            enabled=show_progress,
            total=total_size,
            description=f"Downloading {description}",
        ) as chunks:
            for chunk in chunks:
                f.write(chunk)
