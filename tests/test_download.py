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


class _FakeHttpResponse:
    """Stand-in for an httpx.Response carrying a directory-listing body."""

    def __init__(self, text: str = "") -> None:
        """Store the fake response body."""
        self.text = text

    def raise_for_status(self) -> None:
        """No-op: the fake response is always considered successful."""


class _FakeHttpClient:
    """Stand-in for httpx.Client returning a fixed directory listing."""

    def __init__(self, *args: object, **kwargs: object) -> None:
        """Ignore all init args; this fake never opens a real connection."""

    def __enter__(self) -> _FakeHttpClient:
        """Support use as a context manager, like the real httpx.Client."""
        return self

    def __exit__(self, *args: object) -> None:
        """Support use as a context manager, like the real httpx.Client."""
        return None

    def get(self, url: str) -> _FakeHttpResponse:
        """Return a fixed Ensembl-style directory listing, regardless of *url*."""
        return _FakeHttpResponse(
            '<a href="homo_sapiens_core_115_38/">homo_sapiens_core_115_38/</a>'
        )


class _FailingHttpClient(_FakeHttpClient):
    """Stand-in for httpx.Client whose every request fails."""

    def get(self, url: str) -> _FakeHttpResponse:
        """Simulate a network failure, regardless of *url*."""
        import httpx

        raise httpx.ConnectError("simulated network failure")


class TestEnsemblDownloader:
    """Tests for EnsemblDownloader URL templating and assembly auto-discovery."""

    def test_get_download_urls_templates_version_species_assembly(self) -> None:
        """Explicit assembly skips discovery and fills the URL template directly."""
        from pysec2pri.parsers.ensembl import EnsemblDownloader

        downloader = EnsemblDownloader(
            version="115", species=9606, assembly="38", show_progress=False
        )
        urls = downloader.get_download_urls("115")
        assert urls["stable_id_event"] == (
            "https://ftp.ensembl.org/pub/release-115/mysql/"
            "homo_sapiens_core_115_38/stable_id_event.txt.gz"
        )
        assert urls["gene"].endswith("homo_sapiens_core_115_38/gene.txt.gz")

    def test_get_download_urls_resolves_mouse_species_token(self) -> None:
        """A different taxon ID resolves to its own Ensembl species token."""
        from pysec2pri.parsers.ensembl import EnsemblDownloader

        downloader = EnsemblDownloader(
            version="115", species=10090, assembly="39", show_progress=False
        )
        urls = downloader.get_download_urls("115")
        assert "mus_musculus_core_115_39" in urls["gene"]

    def test_get_download_urls_unknown_species_raises(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A taxon ID unknown both statically and via the live fallback raises."""
        from pysec2pri.parsers import ensembl as ensembl_module

        monkeypatch.setattr("pysec2pri.parsers.ensembl.httpx.Client", _FailingHttpClient)
        downloader = ensembl_module.EnsemblDownloader(
            version="115", species=99999, assembly="1", show_progress=False
        )
        with pytest.raises(ValueError, match="Unknown species"):
            downloader.get_download_urls("115")

    def test_get_download_urls_falls_back_to_live_species_lookup(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A taxon ID absent from config/ensembl.yaml resolves via the live REST fallback."""
        from pysec2pri.parsers import ensembl as ensembl_module

        class _FakeSpeciesResponse(_FakeHttpResponse):
            def json(self) -> dict[str, object]:
                """Return a fake Ensembl REST species list with two dog breed entries."""
                return {
                    "species": [
                        {"name": "canis_lupus_familiaris", "taxon_id": "9615"},
                        {"name": "canis_lupus_familiarisboxer", "taxon_id": "9615"},
                    ]
                }

        class _FakeSpeciesClient(_FakeHttpClient):
            def get(self, url: str) -> _FakeHttpResponse:
                """Return the fake species list for the REST URL, else the base fake response."""
                if "rest.ensembl.org" in url:
                    return _FakeSpeciesResponse()
                return super().get(url)

        monkeypatch.setattr("pysec2pri.parsers.ensembl.httpx.Client", _FakeSpeciesClient)
        downloader = ensembl_module.EnsemblDownloader(
            version="115", species=9615, assembly="1", show_progress=False
        )
        urls = downloader.get_download_urls("115")
        assert "canis_lupus_familiaris_core_115_1" in urls["gene"]
        # The breed-specific entry must not win over the shorter canonical name.
        assert "canis_lupus_familiarisboxer" not in urls["gene"]

    def test_get_download_urls_without_explicit_assembly_calls_discovery(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Omitting assembly triggers _resolve_core_dir (monkeypatched HTTP)."""
        from pysec2pri.parsers import ensembl as ensembl_module

        monkeypatch.setattr("pysec2pri.parsers.ensembl.httpx.Client", _FakeHttpClient)
        downloader = ensembl_module.EnsemblDownloader(
            version="115", species=9606, show_progress=False
        )
        urls = downloader.get_download_urls("115")
        assert "homo_sapiens_core_115_38" in urls["gene"]

    def test_resolve_core_dir_discovers_assembly_from_listing(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """_resolve_core_dir parses the assembly suffix out of the directory listing."""
        from pysec2pri.parsers import ensembl as ensembl_module

        monkeypatch.setattr("pysec2pri.parsers.ensembl.httpx.Client", _FakeHttpClient)
        downloader = ensembl_module.EnsemblDownloader(show_progress=False)
        assert downloader._resolve_core_dir("115", "homo_sapiens") == "38"

    def test_resolve_core_dir_falls_back_to_default_on_http_error(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A failed listing request falls back to parse_options.default_assembly."""
        from pysec2pri.parsers import ensembl as ensembl_module

        monkeypatch.setattr("pysec2pri.parsers.ensembl.httpx.Client", _FailingHttpClient)
        downloader = ensembl_module.EnsemblDownloader(show_progress=False)
        assert downloader._resolve_core_dir("115", "homo_sapiens") == "38"

    def test_resolve_core_dir_falls_back_when_token_not_found(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A listing that doesn't mention the requested token falls back to the default."""
        from pysec2pri.parsers import ensembl as ensembl_module

        monkeypatch.setattr("pysec2pri.parsers.ensembl.httpx.Client", _FakeHttpClient)
        downloader = ensembl_module.EnsemblDownloader(show_progress=False)
        assert downloader._resolve_core_dir("115", "mus_musculus") == "38"
