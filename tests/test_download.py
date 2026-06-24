"""Tests for release-date resolution in the download layer."""

from __future__ import annotations

from datetime import datetime
from types import SimpleNamespace
from typing import Any

import pytest
from sssom_schema import Mapping

from pysec2pri import api as api_module
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


_DOG_LISTING = (
    '<a href="homo_sapiens_core_115_38/">homo_sapiens_core_115_38/</a>\n'
    '<a href="canis_lupus_familiaris_core_115_1/">canis_lupus_familiaris_core_115_1/</a>\n'
    '<a href="canis_lupus_familiarisboxer_core_115_1/">'
    "canis_lupus_familiarisboxer_core_115_1/</a>\n"
)
_DOG_REST_SPECIES = {
    "species": [
        {"name": "canis_lupus_familiaris", "taxon_id": "9615"},
        {"name": "canis_lupus_familiarisboxer", "taxon_id": "9615"},
    ]
}


class TestDiscoverEnsemblSpecies:
    """Tests for the ``species="all"`` bulk-discovery primitive, fully mocked."""

    def test_dedupes_breed_variants_to_canonical_token(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A breed-specific assembly sharing its species' taxon ID is deduped away."""
        from pysec2pri.parsers import ensembl as ensembl_module

        class _RestSpeciesResponse(_FakeHttpResponse):
            def json(self) -> dict[str, list[dict[str, str]]] | dict[str, object]:
                """Return a fake REST species list tagging both dog tokens with one taxon ID."""
                return _DOG_REST_SPECIES

        class _RestClient(_FakeHttpClient):
            def get(self, url: str) -> _FakeHttpResponse:
                """Serve the REST species list, or the FTP listing for any other URL."""
                if "rest.ensembl.org" in url:
                    return _RestSpeciesResponse()
                return _FakeHttpResponse(_DOG_LISTING)

        monkeypatch.setattr("pysec2pri.parsers.ensembl.httpx.Client", _RestClient)

        species = ensembl_module.discover_ensembl_species("115")

        assert set(species) == {("homo_sapiens", "38"), ("canis_lupus_familiaris", "1")}

    def test_falls_back_to_unfiltered_listing_when_rest_api_unreachable(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """If the REST species list can't be fetched, every FTP entry is kept (no dedup)."""
        from pysec2pri.parsers import ensembl as ensembl_module

        class _NoRestClient(_FailingHttpClient):
            def get(self, url: str) -> _FakeHttpResponse:
                """Fail the REST call, but still serve the FTP listing."""
                if "rest.ensembl.org" in url:
                    return super().get(url)
                return _FakeHttpResponse(_DOG_LISTING)

        monkeypatch.setattr("pysec2pri.parsers.ensembl.httpx.Client", _NoRestClient)

        species = ensembl_module.discover_ensembl_species("115")

        assert set(species) == {
            ("homo_sapiens", "38"),
            ("canis_lupus_familiaris", "1"),
            ("canis_lupus_familiarisboxer", "1"),
        }


## Ensembl species orchestration


class _FakeMappingSet:
    """Stand-in for a Sec2PriMappingSet exposing just what the orchestration needs."""

    def __init__(
        self,
        mappings: list[Mapping] | None = None,
        primary_ids: set[str] | None = None,
        primary_labels: dict[str, set[str]] | None = None,
    ) -> None:
        """Store mappings and optional primary-ID/label stores."""
        self.mappings = mappings or []
        self._primary_ids = primary_ids or set()
        self._primary_labels = primary_labels or {}


def _fake_mapping(token: str) -> Mapping:
    """Build a minimal valid Mapping tagged with *token*, for orchestration tests."""
    return Mapping(
        subject_id=f"ENSEMBL:{token}_OLD",
        object_id=f"ENSEMBL:{token}_NEW",
        predicate_id="IAO:0100001",
        mapping_justification="semapv:BackgroundKnowledgeBasedMatching",
    )


def _patch_species_discovery(monkeypatch: pytest.MonkeyPatch, tokens: list[str]) -> None:
    """Stub discovery/release-resolution so every test exercises the real version=None path.

    Every real caller invokes the bulk path with ``version=None``, relying
    on ``check_ensembl_release()`` to resolve "latest" -- rather than
    hardcoding an arbitrary version literal per test, every test here goes
    through that same resolution, just against a fake "latest".
    """
    monkeypatch.setattr(
        "pysec2pri.download.check_ensembl_release", lambda: SimpleNamespace(version="999-test")
    )
    monkeypatch.setattr(
        "pysec2pri.parsers.ensembl.discover_ensembl_species",
        lambda version: [(t, "1") for t in tokens],
    )
    monkeypatch.setattr("pysec2pri.download.resolve_release_date", lambda *a, **kw: None)


class TestGenerateEnsemblAllSpecies:
    """Tests for ``_generate_ensembl_all_species``, fully mocked (no network)."""

    def test_combines_successful_species_and_skips_failures(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Successful species combine into one mapping set; a failed one is skipped."""
        _patch_species_discovery(monkeypatch, ["species_a", "species_b", "species_fail"])

        def fake_process(
            kind: str, version: str, species: str, assembly: str, release_date: object
        ) -> _FakeMappingSet:
            """Succeed for every species except "species_fail", which raises."""
            if species == "species_fail":
                raise RuntimeError("simulated failure")
            return _FakeMappingSet([_fake_mapping(species)], primary_ids={f"ENSEMBL:{species}_NEW"})

        monkeypatch.setattr(api_module, "_process_one_ensembl_species", fake_process)

        result = api_module._generate_ensembl_all_species("ids", None, show_progress=False)

        assert {m.subject_id for m in result.mappings} == {
            "ENSEMBL:species_a_OLD",
            "ENSEMBL:species_b_OLD",
        }
        assert result._primary_ids == {"ENSEMBL:species_a_NEW", "ENSEMBL:species_b_NEW"}
        assert result.mapping_set_id.endswith("/999-test/all")

    def test_merges_primary_labels_across_species(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Per-species ``_primary_labels`` dicts are unioned by label across species."""
        _patch_species_discovery(monkeypatch, ["species_a", "species_b"])
        monkeypatch.setattr(
            api_module,
            "_process_one_ensembl_species",
            lambda kind, version, token, assembly, release_date: _FakeMappingSet(
                primary_labels={"SHAREDSYM": {f"ENSEMBL:{token}_NEW"}}
            ),
        )

        result = api_module._generate_ensembl_all_species("labels", None, show_progress=False)

        assert result._primary_labels["SHAREDSYM"] == {
            "ENSEMBL:species_a_NEW",
            "ENSEMBL:species_b_NEW",
        }

    def test_raises_when_every_species_fails(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """A run where no species could be processed raises, rather than returning empty."""
        _patch_species_discovery(monkeypatch, ["species_a"])

        def fake_process(*args: Any, **kwargs: Any) -> _FakeMappingSet:
            """Fail regardless of which species is requested."""
            raise RuntimeError("simulated failure")

        monkeypatch.setattr(api_module, "_process_one_ensembl_species", fake_process)

        with pytest.raises(ValueError, match="No Ensembl species could be processed"):
            api_module._generate_ensembl_all_species("ids", None, show_progress=False)
