"""Downloading and release detection for Ensembl."""

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

__all__ = [
    "EnsemblDownloader",
    "check_ensembl_release",
    "discover_ensembl_species",
    "urls_and_date",
]

_ENSEMBL_REST_SPECIES_URL = "https://rest.ensembl.org/info/species?content-type=application/json"


def _lookup_ensembl_species_token(taxon_id: str | int) -> str | None:
    """Resolve a taxon ID to its Ensembl species token via Ensembl's own REST API.

    ``config/ensembl.yaml``'s ``species.available`` only curates a handful
    of species for CLI defaults/help; Ensembl itself publishes 300+ (e.g.
    dog, pig, chicken). This is the live fallback for anything not in that
    static list, queried on demand rather than hardcoded, since the list
    changes every release.

    Several taxon IDs map to more than one entry (per-breed/per-strain
    assemblies, e.g. multiple dog breeds): the *shortest* matching name is
    returned, which is reliably the canonical species-level entry (breed/
    strain variants always have additional characters appended).

    Args:
        taxon_id: Canonical NCBI taxon ID.

    Returns:
        The Ensembl species token (e.g. ``"canis_lupus_familiaris"``), or
        ``None`` if the REST API is unreachable or has no matching entry.
    """
    try:
        with httpx.Client(follow_redirects=True, timeout=30.0) as client:
            response = client.get(_ENSEMBL_REST_SPECIES_URL)
            response.raise_for_status()
            data = response.json()
    except httpx.HTTPError:
        return None

    matches = [
        str(entry["name"])
        for entry in data.get("species", [])
        if entry.get("name") and str(entry.get("taxon_id")) == str(taxon_id)
    ]
    return min(matches, key=len) if matches else None


def _ensembl_taxon_by_token() -> dict[str, str]:
    """Return ``{species_token: taxon_id}`` from Ensembl's REST species list.

    Returns ``{}`` if the REST API is unreachable.
    """
    try:
        with httpx.Client(follow_redirects=True, timeout=30.0) as client:
            response = client.get(_ENSEMBL_REST_SPECIES_URL)
            response.raise_for_status()
            data = response.json()
    except httpx.HTTPError:
        return {}
    return {
        str(entry["name"]): str(entry["taxon_id"])
        for entry in data.get("species", [])
        if entry.get("name") and entry.get("taxon_id")
    }


def discover_ensembl_species(version: str) -> list[tuple[str, str]]:
    """Return ``(species_token, assembly)`` for every species published at *version*.

    Scrapes the release's ``mysql/`` directory listing exactly once -- a
    single HTTP request covers every species at once, unlike per-species
    assembly auto-discovery (:meth:`EnsemblDownloader._resolve_core_dir`),
    which would otherwise need one request per species to do the same.
    Used by the ``species="all"`` bulk path (see
    :func:`pysec2pri.api._generate_ensembl_all_species`).

    The raw directory listing includes one entry per *assembly*, not per
    species -- several breed/strain assemblies can share one species (e.g.
    half a dozen dog breeds). Those are deduped to a single canonical
    (shortest-named) token per taxon ID, cross-referenced via Ensembl's
    REST species list, so the bulk run produces one mapping set per
    species rather than one per assembly.

    Args:
        version: Ensembl release number.

    Returns:
        Sorted list of ``(species_token, assembly)`` tuples, e.g.
        ``("homo_sapiens", "38")``.

    Raises:
        httpx.HTTPError: If the release's directory listing can't be
            fetched at all (unlike per-species discovery, there is no
            sensible fallback when the *entire* species list is needed).
    """
    index_url = f"https://ftp.ensembl.org/pub/release-{version}/mysql/"
    with httpx.Client(follow_redirects=True, timeout=60.0) as client:
        response = client.get(index_url)
        response.raise_for_status()

    pattern = re.compile(rf'href="([a-z0-9_]+)_core_{re.escape(str(version))}_(\d+)/"')
    published = dict(sorted(set(pattern.findall(response.text))))

    taxon_by_token = _ensembl_taxon_by_token()
    if not taxon_by_token:
        return sorted(published.items())

    tokens_by_taxon: dict[str, list[str]] = {}
    for token in published:
        key = taxon_by_token.get(token, f"_untagged:{token}")
        tokens_by_taxon.setdefault(key, []).append(token)

    canonical_tokens = {min(tokens, key=len) for tokens in tokens_by_taxon.values()}
    return sorted((tok, published[tok]) for tok in canonical_tokens)


class EnsemblDownloader(BaseDownloader):
    """Downloader for Ensembl core flat-file dumps.

    Download paths are templated by release, species, and assembly suffix
    (e.g. ``homo_sapiens_core_115_38``); the assembly suffix is
    auto-discovered per release via :meth:`_resolve_core_dir`.
    """

    datasource_name = "ensembl"

    def __init__(
        self,
        version: str | None = None,
        species: str | int = 9606,
        assembly: str | None = None,
        show_progress: bool = True,
    ) -> None:
        """Initialize the Ensembl downloader.

        Args:
            version: Ensembl release number (e.g. ``"115"``).
            species: Canonical NCBI taxon ID, resolved to Ensembl's own
                species token via ``config.species_token()``.
            assembly: Force a specific assembly suffix (e.g. ``"38"``),
                skipping auto-discovery.
            show_progress: Whether to show progress bars.
        """
        super().__init__(version=version, show_progress=show_progress)
        self.species = species
        self.assembly = assembly

    def _resolve_core_dir(self, version: str, token: str) -> str:
        r"""Discover the per-species assembly suffix for *version*.

        Lists ``.../release-{version}/mysql/`` and matches
        ``{token}_core_{version}_(\\d+)``; falls back to
        ``parse_options.default_assembly`` when the listing can't be
        fetched or no match is found.

        Args:
            version: Ensembl release number.
            token: Species token (e.g. ``"homo_sapiens"``).

        Returns:
            Assembly suffix string (e.g. ``"38"``).
        """
        default_assembly = "38"
        if self._config:
            default_assembly = str(
                self._config.parse_options.get("default_assembly", default_assembly)
            )

        index_url = f"https://ftp.ensembl.org/pub/release-{version}/mysql/"
        try:
            with httpx.Client(follow_redirects=True, timeout=30.0) as client:
                response = client.get(index_url)
                response.raise_for_status()
        except httpx.HTTPError:
            return default_assembly

        match = re.search(
            rf"{re.escape(token)}_core_{re.escape(str(version))}_(\d+)", response.text
        )
        return match.group(1) if match else default_assembly

    def _resolve_species_token(self, taxon_id: str | int) -> str:
        """Resolve *taxon_id* to its Ensembl species token.

        Tries ``config/ensembl.yaml``'s static ``species.available`` map
        first (no network); falls back to a live lookup against Ensembl's
        REST species list (see :func:`_lookup_ensembl_species_token`) for
        taxon IDs not curated there -- Ensembl publishes 300+ species, of
        which ``available`` only lists a handful for CLI defaults/help.

        Args:
            taxon_id: Canonical NCBI taxon ID.

        Returns:
            The Ensembl species token.

        Raises:
            ValueError: If *taxon_id* is unknown both statically and live.
        """
        if self._config:
            try:
                return self._config.species_token(taxon_id)
            except ValueError:
                pass

        token = _lookup_ensembl_species_token(taxon_id)
        if token is not None:
            return token

        known = ""
        if self._config:
            available = (self._config.species or {}).get("available") or {}
            known = ", ".join(sorted(str(k) for k in available))
        raise ValueError(
            f"Unknown species taxon ID {taxon_id!r} for Ensembl: not in config/ensembl.yaml's "
            f"static list ({known or '(none configured)'}) and not found via Ensembl's live "
            "species list either."
        )

    def get_download_urls(
        self,
        version: str | None = None,
        **kwargs: Any,
    ) -> dict[str, str]:
        """Get Ensembl download URLs for a release/species.

        Args:
            version: Ensembl release number.
            **kwargs: Optional ``species``/``assembly`` overrides.

        Returns:
            Dictionary with file URLs keyed by table name.
        """
        v = version or self.version
        if v is None:
            raise ValueError("Ensembl requires an explicit release version.")
        if not self._config:
            raise ValueError("Ensembl config not loaded")

        species = kwargs.get("species", self.species)
        token = self._resolve_species_token(species)
        assembly = kwargs.get("assembly") or self.assembly or self._resolve_core_dir(str(v), token)

        return {
            key: url.format(version=v, species=token, assembly=assembly)
            for key, url in self._config.download_urls.items()
        }

    def download(
        self,
        output_dir: Path,
        version: str | None = None,
        decompress: bool = True,
        keys: list[str] | None = None,
        **kwargs: Any,
    ) -> dict[str, Path]:
        """Download Ensembl core flat files.

        Args:
            output_dir: Directory to save files.
            version: Ensembl release number.
            decompress: Whether to decompress ``.gz`` files.
            keys: Optional list of file-key names to download (defaults to
                all keys in ``ensembl.yaml``'s ``download_urls``).
            **kwargs: Forwarded to :meth:`get_download_urls`.

        Returns:
            Dictionary mapping file keys to downloaded paths.
        """
        urls = self.get_download_urls(version, **kwargs)
        if keys is not None:
            urls = {k: v for k, v in urls.items() if k in keys}
        return self._download_urls(urls, output_dir, decompress)

    def list_versions(self) -> list[str]:
        """List all Ensembl release numbers published on the FTP server.

        Returns:
            Sorted list of release-number strings, e.g. ``["100", ..., "115"]``.
        """
        with httpx.Client(follow_redirects=True, timeout=30.0) as client:
            response = client.get("https://ftp.ensembl.org/pub/")
            response.raise_for_status()
        versions = sorted({int(m) for m in re.findall(r"release-(\d+)/", response.text)})
        return [str(v) for v in versions]


def _get_ensembl_urls_for_version(
    version: str,
    species: str | int = 9606,
    assembly: str | None = None,
) -> dict[str, str]:
    """Build Ensembl URLs for a specific release/species.

    Delegates to :meth:`EnsemblDownloader.get_download_urls`, which resolves
    *species* to Ensembl's own species token and auto-discovers the assembly
    suffix when *assembly* is not given.

    Args:
        version: Ensembl release number (e.g. ``"115"``).
        species: Canonical NCBI taxon ID.
        assembly: Force a specific assembly suffix, skipping auto-discovery.

    Returns:
        Dictionary with file URLs keyed by table name.
    """
    return EnsemblDownloader(version=version, species=species, assembly=assembly).get_download_urls(
        version
    )


def check_ensembl_release() -> ReleaseInfo:
    """Check for the latest Ensembl release.

    The release number is global across every species (Ensembl cuts all
    species' core databases under the same release number), so this always
    checks against the default species (human).

    Returns:
        ReleaseInfo with the latest Ensembl release details.
    """
    versions = EnsemblDownloader(show_progress=False).list_versions()
    if not versions:
        raise ValueError("Could not find Ensembl releases on the FTP server")
    version = versions[-1]

    urls = _get_ensembl_urls_for_version(version)
    check_url = urls.get("stable_id_event")
    release_date = get_file_last_modified(check_url) if check_url else None

    return ReleaseInfo(
        datasource="ensembl",
        version=version,
        release_date=release_date,
        is_new=True,
        files=urls,
    )


def urls_and_date(
    version: str | None, config: DatasourceConfig, **kwargs: Any
) -> tuple[dict[str, str], datetime | None]:
    """Resolve Ensembl download URLs and the source release date.

    Args:
        version: Ensembl release number, or ``None`` for latest.
        config: Ensembl's datasource configuration.
        **kwargs: Optional ``species`` override.

    Returns:
        Tuple of (file-key -> URL mapping, release date or None).
    """
    species = kwargs.get("species", config.default_species())
    release_date: datetime | None = None
    if version:
        urls = _get_ensembl_urls_for_version(version, species=species)
        logger.info("Ensembl version %s (species %s): %s", version, species, urls)
    else:
        # The release number is global across species, but the URLs (and the
        # resolved release date) are species-specific, so resolve the latest
        # version and build URLs for the requested species rather than reusing
        # check_ensembl_release()'s human-only URLs.
        versions = EnsemblDownloader(show_progress=False).list_versions()
        if not versions:
            raise ValueError("Could not find Ensembl releases on the FTP server")
        version = versions[-1]
        urls = _get_ensembl_urls_for_version(version, species=species)
        check_url = urls.get("stable_id_event")
        release_date = get_file_last_modified(check_url) if check_url else None
        logger.info("Ensembl release %s (species %s): %s", version, species, urls)
    return urls, release_date
