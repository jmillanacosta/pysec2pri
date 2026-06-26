"""Downloading and release detection for HMDB."""

from __future__ import annotations

from datetime import datetime
from typing import TYPE_CHECKING, Any

from mapkgsutils.download import ReleaseInfo, get_file_last_modified

from pysec2pri.constants import ALL_DATASOURCES
from pysec2pri.logging import logger

if TYPE_CHECKING:
    from pysec2pri.parsers.base import DatasourceConfig

__all__ = ["check_hmdb_release", "urls_and_date"]

# HMDB's live site sits behind Cloudflare bot protection, which blocks even
# the lightweight HEAD request the generic Last-Modified fallback relies on.
# HMDB also has not published a new release since Version 5.0, so a static
# lookup is both necessary and safe to keep current. Dates are the "Released
# on" values published on HMDB's own Downloads page
# (https://hmdb.ca/downloads, Version 5.0 / "Current Version" tab), keyed by
# each config's ``config_id`` ("hmdb_metabolites"/"hmdb_proteins").
_RELEASE_DATES: dict[str, datetime] = {
    "hmdb_metabolites": datetime(2021, 11, 17),
    "hmdb_proteins": datetime(2021, 11, 9),
}


def check_hmdb_release() -> ReleaseInfo:
    """Check for the latest HMDB release by downloading and checking XML.

    Returns:
        ReleaseInfo with the latest HMDB release details.
    """
    # HMDB requires downloading to check the version in XML
    metabolites_url = ALL_DATASOURCES["hmdb_metabolites"].download_urls["metabolites"]
    last_modified = get_file_last_modified(metabolites_url)
    version = last_modified.strftime("%Y-%m-%d") if last_modified else None

    return ReleaseInfo(
        datasource="hmdb",
        version=version,
        release_date=last_modified,
        is_new=True,
        files=dict(ALL_DATASOURCES["hmdb_metabolites"].download_urls),
    )


def urls_and_date(
    version: str | None, config: DatasourceConfig, **kwargs: Any
) -> tuple[dict[str, str], datetime | None]:
    """Resolve HMDB download URLs and the source release date.

    HMDB sits behind Cloudflare bot protection (blocks even the generic
    Last-Modified fallback) and has not published a new release since
    Version 5.0, so the release date comes from the static
    :data:`_RELEASE_DATES` lookup rather than being resolved live.
    """
    if version:
        logger.warning(
            "%s does not have versioned archives. Downloading latest version instead.",
            config.name,
        )
    return dict(config.download_urls), _RELEASE_DATES[config.config_id]
