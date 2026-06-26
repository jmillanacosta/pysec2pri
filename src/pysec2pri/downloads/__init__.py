"""Per-datasource downloader classes and release/URL resolution.

Each submodule holds one datasource's concrete download knowledge (its
:class:`~mapkgsutils.parsers.base.BaseDownloader` subclass, its
``check_<source>_release`` function, and its URL/release-date resolver).
The registries below feed that knowledge into :mod:`mapkgsutils.download`'s
datasource-agnostic dispatchers.
"""

from __future__ import annotations

from pysec2pri.downloads import chebi, ensembl, hgnc, hmdb, ncbi, uniprot
from pysec2pri.downloads.chebi import ChEBIDownloader
from pysec2pri.downloads.ensembl import EnsemblDownloader
from pysec2pri.downloads.hgnc import HGNCDownloader
from pysec2pri.downloads.uniprot import UniProtDownloader
from pysec2pri.parsers.base import BaseDownloader

#: ``{datasource_name: check_release_fn}``, for
#: :func:`mapkgsutils.download.get_latest_release_info`.
CHECK_RELEASE = {
    "chebi": chebi.check_chebi_release,
    "ensembl": ensembl.check_ensembl_release,
    "hgnc": hgnc.check_hgnc_release,
    "hmdb": hmdb.check_hmdb_release,
    "ncbi": ncbi.check_ncbi_release,
    "uniprot": uniprot.check_uniprot_release,
}

#: ``{datasource_name: downloader_class}``, for :func:`mapkgsutils.download.list_versions`.
#: Datasources with no versioned archive (NCBI, VGNC, Wikidata, HMDB) are absent.
DOWNLOADERS: dict[str, type[BaseDownloader]] = {
    "chebi": chebi.ChEBIDownloader,
    "ensembl": ensembl.EnsemblDownloader,
    "hgnc": hgnc.HGNCDownloader,
    "uniprot": uniprot.UniProtDownloader,
}

#: ``{datasource_name: urls_and_date_fn}``, for
#: :func:`mapkgsutils.download.resolve_datasource_urls`.
#: Datasources absent here (NCBI, VGNC, Wikidata) use the generic
#: "latest config URLs, no release date" fallback.
URLS_AND_DATE = {
    "chebi": chebi.urls_and_date,
    "ensembl": ensembl.urls_and_date,
    "hgnc": hgnc.urls_and_date,
    "uniprot": uniprot.urls_and_date,
    "hmdb_metabolites": hmdb.urls_and_date,
    "hmdb_proteins": hmdb.urls_and_date,
}

#: ``{datasource_name: extractor}``, for
#: :func:`mapkgsutils.download.download_datasource_with_release`.
#: Only UniProt publishes a ``.tar.gz`` archive needing member-level extraction.
TAR_EXTRACTORS = {
    "uniprot": uniprot.extract_tar,
}

__all__ = [
    "CHECK_RELEASE",
    "DOWNLOADERS",
    "TAR_EXTRACTORS",
    "URLS_AND_DATE",
    "ChEBIDownloader",
    "EnsemblDownloader",
    "HGNCDownloader",
    "UniProtDownloader",
]
