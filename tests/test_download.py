"""Tests for release-date resolution in the download layer."""

from __future__ import annotations

from datetime import datetime

import pytest

from pysec2pri.download import CloudflareBlockedError, _get_datasource_urls, resolve_release_date


class TestHmdbStaticReleaseDates:
    """HMDB sits behind Cloudflare and rarely releases, so its dates are static.

    See ``_HMDB_RELEASE_DATES`` in ``download.py`` for the source page.
    """

    def test_hmdb_metabolites_uses_static_date(self) -> None:
        """``hmdb_metabolites`` resolves to its known Version 5.0 release date."""
        date = resolve_release_date("hmdb_metabolites")
        assert date == datetime(2021, 11, 17)

    def test_hmdb_proteins_uses_static_date(self) -> None:
        """``hmdb_proteins`` resolves to its known Version 5.0 release date."""
        date = resolve_release_date("hmdb_proteins")
        assert date == datetime(2021, 11, 9)


class TestGenericFallbackResilience:
    """The generic Last-Modified fallback must never crash the download."""

    def test_cloudflare_block_degrades_to_none(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """A Cloudflare-blocked Last-Modified lookup yields ``None``, not a raise."""
        from pysec2pri import download as download_module
        from pysec2pri.constants import ALL_DATASOURCES

        def _raise_blocked(url: str, timeout: float = 30.0) -> datetime | None:
            raise CloudflareBlockedError(f"blocked: {url}")

        monkeypatch.setattr(download_module, "get_file_last_modified", _raise_blocked)
        urls, release_date = _get_datasource_urls("ncbi", ALL_DATASOURCES["ncbi"])
        assert release_date is None
        assert urls
