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

from pysec2pri.constants import (
    ALL_DATASOURCES,
    CHEBI,
    HGNC,
    HMDB,
    NCBI,
    UNIPROT,
)
from pysec2pri.logging import logger

__all__ = [
    "ReleaseInfo",
    "check_release",
    "download_datasource",
    "download_file",
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
    archive_url = CHEBI.archive_url
    with httpx.Client(follow_redirects=True) as client:
        response = client.get(archive_url)
        response.raise_for_status()

    # Parse the archive index to find the latest release number
    matches = re.findall(r'href="rel(\d+)/"', response.text)
    if not matches:
        raise ValueError("Could not find ChEBI releases in archive")

    latest_release = max(int(m) for m in matches)
    version = str(latest_release)

    sdf_url = (
        f"http://ftp.ebi.ac.uk/pub/databases/chebi/archive/"
        f"rel{latest_release}/SDF/chebi_3_stars.sdf.gz"
    )

    return ReleaseInfo(
        datasource="chebi",
        version=version,
        release_date=get_file_last_modified(sdf_url),
        is_new=True,  # Caller determines if it's new
        files={"sdf": sdf_url},
    )


def check_ncbi_release() -> ReleaseInfo:
    """Check for the latest NCBI Gene release.

    Returns:
        ReleaseInfo with the latest NCBI release details.
    """
    history_url = NCBI.download_urls["gene_history"]
    last_modified = get_file_last_modified(history_url)
    version = last_modified.strftime("%Y-%m-%d") if last_modified else None

    return ReleaseInfo(
        datasource="ncbi",
        version=version,
        release_date=last_modified,
        is_new=True,
        files=dict(NCBI.download_urls),
    )


def check_uniprot_release() -> ReleaseInfo:
    """Check for the latest UniProt release.

    Returns:
        ReleaseInfo with the latest UniProt release details.
    """
    sec_ac_url = UNIPROT.download_urls["secondary"]
    last_modified = get_file_last_modified(sec_ac_url)
    version = last_modified.strftime("%Y-%m-%d") if last_modified else None

    return ReleaseInfo(
        datasource="uniprot",
        version=version,
        release_date=last_modified,
        is_new=True,
        files=dict(UNIPROT.download_urls),
    )


def check_hmdb_release() -> ReleaseInfo:
    """Check for the latest HMDB release by downloading and checking XML.

    Returns:
        ReleaseInfo with the latest HMDB release details.
    """
    # HMDB requires downloading to check the version in XML
    metabolites_url = HMDB.download_urls["metabolites"]
    last_modified = get_file_last_modified(metabolites_url)
    version = last_modified.strftime("%Y-%m-%d") if last_modified else None

    return ReleaseInfo(
        datasource="hmdb",
        version=version,
        release_date=last_modified,
        is_new=True,
        files=dict(HMDB.download_urls),
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
        last_modified = get_file_last_modified(HGNC.download_urls["complete"])
        version = last_modified.strftime("%Y-%m-%d") if last_modified else None
        return ReleaseInfo(
            datasource="hgnc",
            version=version,
            release_date=last_modified,
            is_new=True,
            files=dict(HGNC.download_urls),
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


def download_datasource(
    datasource: str,
    output_dir: Path,
    decompress: bool = True,
) -> dict[str, Path]:
    """Download all files for a datasource.

    For datasources with dynamic URLs (like HGNC quarterly archive),
    this function first checks for the latest release and uses those URLs.

    Args:
        datasource: Name of the datasource.
        output_dir: Directory to save files.
        decompress: Whether to decompress .gz files.

    Returns:
        Dictionary mapping file keys to downloaded paths.
    """
    datasource_lower = datasource.lower()
    config = ALL_DATASOURCES.get(datasource_lower)
    if not config:
        raise ValueError(f"Unknown datasource: {datasource}")

    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded = {}

    # For HGNC, get dynamic URLs from the archive
    if datasource_lower == "hgnc":
        release_info = check_hgnc_release()
        urls = release_info.files
        logger.info(f"HGNC release {release_info.version}: {urls}")
    else:
        urls = config.download_urls

    for key, url in urls.items():
        # Determine output filename
        filename = url.split("/")[-1]
        if decompress and filename.endswith(".gz"):
            filename = filename[:-3]

        output_path = output_dir / filename
        logger.info(f"Downloading {key}: {url}")
        download_file(url, output_path, decompress_gz=decompress)
        downloaded[key] = output_path
        logger.info(f"Saved to: {output_path}")

    return downloaded


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
