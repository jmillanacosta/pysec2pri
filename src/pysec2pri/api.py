"""Main functions for pysec2pri.

This module provides functions for parsing biological database
secondary-to-primary mapping files and generating and using the standardized Mapping sets.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

from pysec2pri.context import ContextSpec, load_xref_mapping
from pysec2pri.exports import (
    write_json,
    write_name2synonym,
    write_output,
    write_owl,
    write_pri_ids,
    write_rdf,
    write_sec2pri,
    write_sssom,
)
from pysec2pri.exports import (
    write_label2prev as write_label_sec2pri,
)

if TYPE_CHECKING:
    from datetime import datetime

    import pandas as pd

    from pysec2pri.context import DecisionRecord, XrefMapping
    from pysec2pri.diff import MappingDiff
    from pysec2pri.parsers.base import Sec2PriMappingSet

from pysec2pri.parsers.base import AmbiguousMappingSet, IdMappingSet, LabelMappingSet

__all__ = [
    # Utilities
    "ContextSpec",
    "combine_mapping_sets",
    "crosswalk",
    "find_ambiguous",
    # Need to remove at some point the old functions
    # (see aliased functions at the bottom)
    "generate_chebi",
    # Core
    "generate_chebi_ids",
    "generate_chebi_labels",
    "generate_chebi_primary_ids",
    "generate_chebi_primary_labels",
    "generate_chebi_synonyms",
    "generate_ensembl",
    "generate_ensembl_ids",
    "generate_ensembl_label_history",
    "generate_ensembl_labels",
    "generate_ensembl_primary_ids",
    "generate_ensembl_primary_labels",
    "generate_hgnc",
    "generate_hgnc_ids",
    "generate_hgnc_labels",
    "generate_hgnc_labels",
    "generate_hgnc_primary_ids",
    "generate_hmdb",
    "generate_hmdb_ids",
    "generate_hmdb_primary_ids",
    "generate_hmdb_proteins",
    "generate_hmdb_proteins_ids",
    "generate_ncbi",
    "generate_ncbi_ids",
    "generate_ncbi_labels",
    "generate_ncbi_labels",
    "generate_ncbi_primary_ids",
    "generate_ncbi_primary_labels",
    "generate_uniprot",
    "generate_uniprot_ids",
    "generate_uniprot_primary_ids",
    "generate_vgnc",
    "generate_vgnc_ids",
    "generate_vgnc_labels",
    "generate_vgnc_primary_ids",
    "generate_vgnc_primary_labels",
    "generate_wikidata",
    "generate_wikidata_ids",
    "generate_wikidata_labels",
    "generate_wikidata_labels",
    "list_versions",
    "load_label_mapping",
    "load_mapping",
    "load_xref_mapping",
    "resolve_ids",
    "resolve_labels",
    "save",
    "write_all_formats",
    "write_diff_output",
    "write_json",
    "write_label_sec2pri",
    "write_name2synonym",
    "write_output",
    "write_owl",
    "write_rdf",
    "write_sec2pri",
    "write_sssom",
]


def _auto_download(
    datasource: str,
    version: str | None = None,
    keys: list[str] | None = None,
) -> tuple[dict[str, Path], datetime | None]:
    """Download files for *datasource* into a temp dir.

    Args:
        datasource: Datasource name (e.g. ``"hgnc"``).
        version: Optional specific version to download.
        keys: Optional list of file-key names to download. When given, only
            those keys are fetched (e.g. ``["complete"]``).

    Returns:
        Tuple of (file-key -> downloaded path mapping, source release date or
        None). The release date is set on the parser so the generated mapping
        set's ``mapping_date`` reflects the upstream release.
    """
    import tempfile

    from pysec2pri.download import download_datasource_with_release

    tmpdir = Path(tempfile.mkdtemp(prefix=f"pysec2pri_{datasource}_"))
    return download_datasource_with_release(datasource, tmpdir, version=version, keys=keys)


def generate_chebi(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
    subset: str = "3star",
    mapping_sets: str = "ids",
) -> Sec2PriMappingSet:
    """Return ChEBI mappings (IDs, synonyms, or both).

    Downloads the latest release automatically when ``input_path`` is omitted.
    Pass an SDF file (releases < 245) or a directory of TSV flat files
    (releases >= 245) to use a local copy.

    Args:
        input_path: Local SDF file or TSV directory. Auto-downloaded if ``None``.
        version: Release number (e.g. ``"245"``).
        show_progress: Whether to show progress bars.
        subset: ``"3star"`` (default) or ``"complete"``.
        mapping_sets: ``"ids"`` (default), ``"synonyms"``, or ``"all"``.
    """
    import tempfile

    from pysec2pri.download import check_chebi_release, resolve_release_date
    from pysec2pri.parsers import ChEBIParser
    from pysec2pri.parsers.chebi import ChEBIDownloader

    release_date = None
    if input_path is None:
        if version is None:
            version = check_chebi_release().version or "245"
        release_date = resolve_release_date("chebi", version, subset=subset)
        downloader = ChEBIDownloader(version=version, subset=subset)
        tmpdir = Path(tempfile.mkdtemp(prefix=f"pysec2pri_chebi_{version}_"))
        downloaded = downloader.download(tmpdir, version=version)
        # SDF releases (< 245) download a single file; TSV releases (>= 245)
        # download a directory of flat files that `parse()` auto-discovers.
        input_path = downloaded.get("sdf", tmpdir)

    parser = ChEBIParser(version=version, show_progress=show_progress, subset=subset)
    parser.release_date = release_date

    if mapping_sets == "synonyms":
        return parser.parse_synonyms(Path(input_path))
    if mapping_sets == "all":
        ids = parser.parse(Path(input_path))
        syns = parser.parse_synonyms(Path(input_path))
        return combine_mapping_sets(ids, syns)
    return parser.parse(Path(input_path))


def generate_chebi_synonyms(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
    subset: str = "3star",
) -> Sec2PriMappingSet:
    """Return ChEBI synonym (name) mappings."""
    return generate_chebi(
        input_path=input_path,
        version=version,
        show_progress=show_progress,
        subset=subset,
        mapping_sets="synonyms",
    )


def generate_hgnc(
    input_path: Path | str | None = None,
    complete_set_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return HGNC secondary to primary ID mappings.

    Downloads the withdrawn and complete set files automatically when
    ``input_path`` / ``complete_set_path`` are omitted.  The complete set is
    used to populate the full list of current primary IDs so that
    :meth:`~pysec2pri.parsers.base.Sec2PriMappingSet.to_pri_ids` returns the
    authoritative list (~45 k IDs) rather than just the ~5 k primaries that
    happen to have a secondary.

    Args:
        input_path: Local HGNC withdrawn TSV. Auto-downloaded if ``None``.
        complete_set_path: Local HGNC complete set TSV. Auto-downloaded if
            ``None``.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import HGNCParser

    release_date = None
    if input_path is None or complete_set_path is None:
        files, release_date = _auto_download("hgnc", version)
        if input_path is None:
            input_path = files["withdrawn"]
        if complete_set_path is None:
            complete_set_path = files["complete"]

    parser = HGNCParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse(Path(input_path), complete_set_path=Path(complete_set_path))


def generate_hgnc_primary_ids(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return a mapping set containing the full list of current HGNC primary IDs.

    Only the HGNC complete set file is downloaded/read.  The returned mapping
    set has an empty ``mappings`` list; its ``_primary_ids`` store is
    populated with every current HGNC ID so that ``to_pri_ids()`` produces
    the authoritative complete list, not just the subset of primaries that
    happen to have an associated secondary.

    Args:
        input_path: Local HGNC complete set TSV. Auto-downloaded if ``None``.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import HGNCParser

    release_date = None
    if input_path is None:
        files, release_date = _auto_download("hgnc", version, keys=["complete"])
        input_path = files["complete"]

    parser = HGNCParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse_primary_ids(Path(input_path))


def generate_hgnc_primary_labels(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return a mapping set containing the full list of current HGNC primary Symbols.

    Only the HGNC complete set file is downloaded/read.  The returned mapping
    set has an empty ``mappings`` list; its ``_primary_labels`` store is
    populated with every current HGNC Symbol so that ``to_pri_labels()`` produces
    the authoritative complete list, not just the subset of primaries that
    happen to have an associated secondary.

    Args:
        input_path: Local HGNC complete set TSV. Auto-downloaded if ``None``.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import HGNCParser

    release_date = None
    if input_path is None:
        files, release_date = _auto_download("hgnc", version, keys=["complete"])
        input_path = files["complete"]

    parser = HGNCParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse_primary_labels(Path(input_path))


def generate_hgnc_labels(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
    statuses: list[str] | None = None,
) -> Sec2PriMappingSet:
    """Return HGNC label to previous-label mappings.

    Downloads the complete set file automatically when ``input_path`` is omitted.

    Args:
        input_path: Local HGNC complete set TSV. Auto-downloaded if ``None``.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
        statuses: Entry statuses to include (e.g. ``["Approved"]``).
    """
    from pysec2pri.parsers import HGNCParser

    release_date = None
    if input_path is None:
        files, release_date = _auto_download("hgnc", version)
        input_path = files["complete"]

    parser = HGNCParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse_labels(Path(input_path), statuses=statuses)


def generate_vgnc(
    input_path: Path | str | None = None,
    complete_set_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return VGNC secondary to primary ID mappings.

    Downloads the withdrawn and gene-set files automatically when
    ``input_path`` / ``complete_set_path`` are omitted. Never filtered by
    species: VGNC IDs are unique across all species (like HGNC IDs), and
    the withdrawn file's ``taxon_id`` column is not populated upstream (see
    :mod:`pysec2pri.parsers.vgnc`). The gene-set file is used to populate the
    full list of current primary IDs (across every species) so that
    :meth:`~pysec2pri.parsers.base.Sec2PriMappingSet.to_pri_ids` returns the
    authoritative list rather than just the primaries that happen to have a
    secondary.

    Args:
        input_path: Local VGNC withdrawn TSV. Auto-downloaded if ``None``.
        complete_set_path: Local VGNC gene-set TSV. Auto-downloaded if
            ``None``.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import VGNCParser

    release_date = None
    if input_path is None or complete_set_path is None:
        files, release_date = _auto_download("vgnc", version)
        if input_path is None:
            input_path = files["withdrawn"]
        if complete_set_path is None:
            complete_set_path = files["complete"]

    parser = VGNCParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse(Path(input_path), complete_set_path=Path(complete_set_path))


def generate_vgnc_primary_ids(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return a mapping set containing the full list of current VGNC primary IDs.

    Only the VGNC gene-set file is downloaded/read. The returned mapping set
    has an empty ``mappings`` list; its ``_primary_ids`` store is populated
    with every current VGNC ID across all species (not filtered by species,
    see :mod:`pysec2pri.parsers.vgnc`).

    Args:
        input_path: Local VGNC gene-set TSV. Auto-downloaded if ``None``.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import VGNCParser

    release_date = None
    if input_path is None:
        files, release_date = _auto_download("vgnc", version, keys=["complete"])
        input_path = files["complete"]

    parser = VGNCParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse_primary_ids(Path(input_path))


def generate_vgnc_primary_labels(
    input_path: Path | str | None = None,
    species: str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return a mapping set with the full list of current VGNC primary symbols for one species.

    Only the VGNC gene-set file is downloaded/read. The returned mapping set
    has an empty ``mappings`` list; its ``_primary_labels`` store is
    populated with every current approved symbol for ``species``.

    Args:
        input_path: Local VGNC gene-set TSV. Auto-downloaded if ``None``.
        species: NCBI taxon ID to filter by. Defaults to ``config/vgnc.yaml``'s
            ``species.default`` (chimpanzee) when omitted.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.constants import VGNC
    from pysec2pri.parsers import VGNCParser

    species = species if species is not None else str(VGNC.default_species())

    release_date = None
    if input_path is None:
        files, release_date = _auto_download("vgnc", version, keys=["complete"])
        input_path = files["complete"]

    parser = VGNCParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse_primary_labels(Path(input_path), species=species)


def generate_vgnc_labels(
    input_path: Path | str | None = None,
    species: str | None = None,
    version: str | None = None,
    show_progress: bool = True,
    statuses: list[str] | None = None,
) -> Sec2PriMappingSet:
    """Return VGNC label to previous-label mappings for one species.

    Downloads the gene-set file automatically when ``input_path`` is
    omitted. Symbols are not unique across species (the same approved
    symbol can legitimately name orthologous genes in different species),
    so ``species`` scopes ambiguity detection and ``to_pri_labels()`` to one
    species' own namespace (see :mod:`pysec2pri.parsers.vgnc`).

    Args:
        input_path: Local VGNC gene-set TSV. Auto-downloaded if ``None``.
        species: NCBI taxon ID to filter by. Defaults to ``config/vgnc.yaml``'s
            ``species.default`` (chimpanzee) when omitted.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
        statuses: Entry statuses to include (e.g. ``["Approved"]``).
    """
    from pysec2pri.constants import VGNC
    from pysec2pri.parsers import VGNCParser

    species = species if species is not None else str(VGNC.default_species())

    release_date = None
    if input_path is None:
        files, release_date = _auto_download("vgnc", version)
        input_path = files["complete"]

    parser = VGNCParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse_labels(Path(input_path), species=species, statuses=statuses)


def generate_chebi_primary_ids(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
    subset: str = "3star",
) -> Sec2PriMappingSet:
    """Return a mapping set containing the full list of current ChEBI primary IDs.

    Reads ``compounds.tsv`` to extract every current ChEBI compound ID.
    The returned mapping set has an empty ``mappings`` list; ``_primary_ids``
    is populated with every current ``CHEBI:<n>`` CURIE.

    Args:
        input_path: Local ``compounds.tsv`` file or directory containing it.
            Auto-downloaded if ``None``.
        version: Release number (e.g. ``"245"``).
        show_progress: Whether to show progress bars.
        subset: ``"3star"`` (default) or ``"complete"``.
    """
    import tempfile

    from pysec2pri.download import check_chebi_release, resolve_release_date
    from pysec2pri.parsers.chebi import ChEBIDownloader, ChEBIParser

    release_date = None
    if input_path is None:
        if version is None:
            version = check_chebi_release().version or "245"
        release_date = resolve_release_date("chebi", version, subset=subset)
        downloader = ChEBIDownloader(version=version, subset=subset)
        tmpdir = Path(tempfile.mkdtemp(prefix=f"pysec2pri_chebi_{version}_"))
        downloader.download(tmpdir, version=version, keys=["compounds"])
        input_path = tmpdir

    parser = ChEBIParser(version=version, show_progress=show_progress, subset=subset)
    parser.release_date = release_date
    return parser.parse_primary_ids(Path(input_path))


def generate_chebi_primary_labels(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
    subset: str = "3star",
) -> Sec2PriMappingSet:
    """Return a mapping set containing the full list of current ChEBI compound names.

    Reads ``compounds.tsv`` to extract every current compound's canonical name.
    The returned mapping set has an empty ``mappings`` list; ``_primary_labels``
    is populated.

    Args:
        input_path: Local ``compounds.tsv`` file or directory containing it.
            Auto-downloaded if ``None``.
        version: Release number (e.g. ``"245"``).
        show_progress: Whether to show progress bars.
        subset: ``"3star"`` (default) or ``"complete"``.
    """
    import tempfile

    from pysec2pri.download import check_chebi_release, resolve_release_date
    from pysec2pri.parsers.chebi import ChEBIDownloader, ChEBIParser

    release_date = None
    if input_path is None:
        if version is None:
            version = check_chebi_release().version or "245"
        release_date = resolve_release_date("chebi", version, subset=subset)
        downloader = ChEBIDownloader(version=version, subset=subset)
        tmpdir = Path(tempfile.mkdtemp(prefix=f"pysec2pri_chebi_{version}_"))
        downloader.download(tmpdir, version=version, keys=["compounds"])
        input_path = tmpdir

    parser = ChEBIParser(version=version, show_progress=show_progress, subset=subset)
    parser.release_date = release_date
    return parser.parse_primary_labels(Path(input_path))


def generate_ncbi_primary_ids(
    input_path: Path | str | None = None,
    species: str = "9606",
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return a mapping set containing the full list of current NCBI Gene primary IDs.

    Reads ``gene_info`` to extract every current Gene ID for the given taxonomy.
    The returned mapping set has an empty ``mappings`` list; ``_primary_ids``
    is populated with every current ``NCBIGene:<id>`` CURIE.

    Args:
        input_path: Local gene_info file. Auto-downloaded if ``None``.
        species: NCBI taxon ID to filter by (default: ``"9606"`` for human).
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import NCBIParser

    release_date = None
    if input_path is None:
        files, release_date = _auto_download("ncbi", version, keys=["gene_info"])
        input_path = files["gene_info"]

    parser = NCBIParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse_primary_ids(Path(input_path), species=species)


def generate_ncbi_primary_labels(
    input_path: Path | str | None = None,
    species: str = "9606",
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return a mapping set containing the full list of current NCBI Gene labels.

    Reads ``gene_info`` to extract every current gene label for the given
    taxonomy.  The returned mapping set has an empty ``mappings`` list;
    ``_primary_labels`` is populated.

    Args:
        input_path: Local gene_info file. Auto-downloaded if ``None``.
        species: NCBI taxon ID to filter by (default: ``"9606"`` for human).
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import NCBIParser

    release_date = None
    if input_path is None:
        files, release_date = _auto_download("ncbi", version, keys=["gene_info"])
        input_path = files["gene_info"]

    parser = NCBIParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse_primary_labels(Path(input_path), species=species)


def generate_hmdb_primary_ids(
    metabolites_path: Path | str | None = None,
    proteins_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return a mapping set containing the full list of current HMDB primary IDs.

    Reads one or both of ``hmdb_metabolites.xml`` and ``hmdb_proteins.xml``
    and collects all primary accession numbers.  The returned mapping set has
    an empty ``mappings`` list; ``_primary_ids`` is populated with every
    current ``HMDB:<acc>`` CURIE.

    Args:
        metabolites_path: Local metabolites XML file. Auto-downloaded if both
            paths are ``None``.
        proteins_path: Local proteins XML file (optional).
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers.hmdb import HMDBMetaboliteParser

    release_date = None
    if metabolites_path is None and proteins_path is None:
        files, release_date = _auto_download("hmdb_metabolites", version, keys=["metabolites"])
        metabolites_path = files["metabolites"]  # Metabolites default
    parser = HMDBMetaboliteParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse_primary_ids(
        metabolites_path=metabolites_path,
        proteins_path=proteins_path,
    )


def generate_uniprot_primary_ids(
    acindex_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return a mapping set containing the full list of current UniProt primary ACs.

    Parses ``acindex.txt`` to extract every accession number that currently
    appears in UniProtKB/Swiss-Prot.  The returned mapping set has an empty
    ``mappings`` list; ``_primary_ids`` is populated with every current
    ``UniProtKB:<AC>`` CURIE.

    For versioned (legacy) releases the file is available at::

        https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/
        release-{version}/knowledgebase/docs/acindex.txt.gz

    Args:
        acindex_path: Local ``acindex.txt`` (plain or ``.gz``).
            Auto-downloaded from the current release when ``None``.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers.uniprot import UniProtParser

    release_date = None
    if acindex_path is None:
        files, release_date = _auto_download("uniprot", version, keys=["acindex"])
        acindex_path = files["acindex"]

    parser = UniProtParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse_primary_ids(Path(acindex_path))


def generate_ncbi(
    input_path: Path | str | None = None,
    gene_info_path: Path | str | None = None,
    species: str = "9606",
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return NCBI Gene secondary to primary ID mappings.

    Downloads the gene_history file automatically when ``input_path`` is
    omitted.  When ``gene_info_path`` is supplied (or auto-downloaded), the
    full list of current primary IDs is read from ``gene_info`` and stored in
    ``_primary_ids``, so that :meth:`~pysec2pri.parsers.base.Sec2PriMappingSet.to_pri_ids`
    returns the authoritative complete set rather than only the subset of
    primaries that happen to appear in ``gene_history``.

    Args:
        input_path: Local gene_history file. Auto-downloaded if ``None``.
        gene_info_path: Local gene_info file used to populate the full primary
            ID list. Auto-downloaded together with ``input_path`` when both
            are ``None``.
        species: NCBI taxon ID to filter (default: ``"9606"`` for human).
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import NCBIParser

    release_date = None
    if input_path is None or gene_info_path is None:
        files, release_date = _auto_download("ncbi", version)
        if input_path is None:
            input_path = files["gene_history"]
        if gene_info_path is None:
            gene_info_path = files["gene_info"]

    parser = NCBIParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse(Path(input_path), species=species, gene_info_path=Path(gene_info_path))


def generate_ncbi_labels(
    input_path: Path | str | None = None,
    species: str = "9606",
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return NCBI Gene label to previous-label mappings.

    Downloads the gene_info file automatically when ``input_path`` is omitted.

    Args:
        input_path: Local gene_info file. Auto-downloaded if ``None``.
        species: NCBI taxon ID to filter (default: ``"9606"`` for human).
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import NCBIParser

    release_date = None
    if input_path is None:
        files, release_date = _auto_download("ncbi", version)
        input_path = files["gene_info"]

    parser = NCBIParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse_labels(Path(input_path), species=species)


def _auto_download_ensembl(
    version: str | None,
    species: str | int,
    keys: list[str],
    show_progress: bool,
) -> tuple[dict[str, Path], str, datetime | None]:
    """Download Ensembl files into a temp dir; resolve version/release date.

    Downloads are parameterized by species/assembly beyond just ``version``
    (unlike most datasources), so this goes through
    :class:`~pysec2pri.parsers.ensembl.EnsemblDownloader` directly rather
    than the generic :func:`_auto_download` (mirrors the ChEBI pattern).
    """
    import tempfile

    from pysec2pri.download import check_ensembl_release, resolve_release_date
    from pysec2pri.parsers.ensembl import EnsemblDownloader

    if version is None:
        version = check_ensembl_release().version
        if version is None:
            raise ValueError("Could not determine the latest Ensembl release.")
    release_date = resolve_release_date("ensembl", version, species=species)
    downloader = EnsemblDownloader(version=version, species=species, show_progress=show_progress)
    tmpdir = Path(tempfile.mkdtemp(prefix=f"pysec2pri_ensembl_{version}_"))
    files = downloader.download(tmpdir, version=version, keys=keys)
    return files, version, release_date


def generate_ensembl(
    input_path: Path | str | None = None,
    mapping_session_path: Path | str | None = None,
    gene_path: Path | str | None = None,
    species: str | int = 9606,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return Ensembl secondary to primary gene ID mappings.

    Downloads ``stable_id_event``/``mapping_session``/``gene`` automatically
    when their paths are omitted. Each release's ``stable_id_event`` table is
    cumulative, so this describes the whole state of Ensembl gene IDs at
    ``version`` (no cross-release chain resolution; see
    :mod:`pysec2pri.consolidate` for that).

    Args:
        input_path: Local ``stable_id_event`` file (the ``ids`` mapping
            set's ``primary_input``, see ``ensembl.yaml``). Auto-downloaded
            if ``None``.
        mapping_session_path: Local ``mapping_session`` file, used to resolve
            each row's ``mapping_date``. Auto-downloaded together with
            ``input_path`` when both are ``None``.
        gene_path: Local ``gene`` file, used to populate the full primary ID
            list. Auto-downloaded together with the others when ``None``.
        species: Canonical NCBI taxon ID (default ``9606`` for human).
        version: Ensembl release number. Latest release is used when ``None``.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers.ensembl import EnsemblParser

    release_date = None
    if input_path is None or mapping_session_path is None or gene_path is None:
        files, version, release_date = _auto_download_ensembl(
            version, species, ["stable_id_event", "mapping_session", "gene"], show_progress
        )
        if input_path is None:
            input_path = files["stable_id_event"]
        if mapping_session_path is None:
            mapping_session_path = files["mapping_session"]
        if gene_path is None:
            gene_path = files["gene"]

    parser = EnsemblParser(version=version, show_progress=show_progress, species=species)
    parser.release_date = release_date
    return parser.parse(
        Path(input_path),
        mapping_session_path=Path(mapping_session_path),
        gene_path=Path(gene_path),
    )


def generate_ensembl_labels(
    input_path: Path | str | None = None,
    gene_path: Path | str | None = None,
    xref_path: Path | str | None = None,
    species: str | int = 9606,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return Ensembl gene external-synonym to current-label mappings.

    Downloads ``external_synonym``/``gene``/``xref`` automatically when
    their paths are omitted.

    Args:
        input_path: Local ``external_synonym`` file (the ``labels`` mapping
            set's ``primary_input``, see ``ensembl.yaml``). Auto-downloaded
            if ``None``.
        gene_path: Local ``gene`` file. Auto-downloaded if ``None``.
        xref_path: Local ``xref`` file. Auto-downloaded if ``None``.
        species: Canonical NCBI taxon ID (default ``9606`` for human).
        version: Ensembl release number. Latest release is used when ``None``.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers.ensembl import EnsemblParser

    release_date = None
    if input_path is None or gene_path is None or xref_path is None:
        files, version, release_date = _auto_download_ensembl(
            version, species, ["gene", "xref", "external_synonym"], show_progress
        )
        if input_path is None:
            input_path = files["external_synonym"]
        if gene_path is None:
            gene_path = files["gene"]
        if xref_path is None:
            xref_path = files["xref"]

    parser = EnsemblParser(version=version, show_progress=show_progress, species=species)
    parser.release_date = release_date
    return parser.parse_labels(Path(gene_path), Path(xref_path), Path(input_path))


def generate_ensembl_primary_ids(
    gene_path: Path | str | None = None,
    species: str | int = 9606,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return a mapping set containing the full list of current Ensembl gene IDs.

    Only the ``gene`` file is downloaded/read. The returned mapping set has
    an empty ``mappings`` list; ``_primary_ids`` is populated.

    Args:
        gene_path: Local ``gene`` file. Auto-downloaded if ``None``.
        species: Canonical NCBI taxon ID (default ``9606`` for human).
        version: Ensembl release number. Latest release is used when ``None``.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers.ensembl import EnsemblParser

    release_date = None
    if gene_path is None:
        files, version, release_date = _auto_download_ensembl(
            version, species, ["gene"], show_progress
        )
        gene_path = files["gene"]

    parser = EnsemblParser(version=version, show_progress=show_progress, species=species)
    parser.release_date = release_date
    return parser.parse_primary_ids(Path(gene_path))


def generate_ensembl_primary_labels(
    gene_path: Path | str | None = None,
    xref_path: Path | str | None = None,
    species: str | int = 9606,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return a mapping set containing the full list of current Ensembl gene labels.

    Only the ``gene``/``xref`` files are downloaded/read. The returned
    mapping set has an empty ``mappings`` list; ``_primary_labels`` is
    populated.

    Args:
        gene_path: Local ``gene`` file. Auto-downloaded if ``None``.
        xref_path: Local ``xref`` file. Auto-downloaded if ``None``.
        species: Canonical NCBI taxon ID (default ``9606`` for human).
        version: Ensembl release number. Latest release is used when ``None``.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers.ensembl import EnsemblParser

    release_date = None
    if gene_path is None or xref_path is None:
        files, version, release_date = _auto_download_ensembl(
            version, species, ["gene", "xref"], show_progress
        )
        if gene_path is None:
            gene_path = files["gene"]
        if xref_path is None:
            xref_path = files["xref"]

    parser = EnsemblParser(version=version, show_progress=show_progress, species=species)
    parser.release_date = release_date
    return parser.parse_primary_labels(Path(gene_path), Path(xref_path))


def generate_ensembl_label_history(
    species: str | int = 9606,
    from_version: str | None = None,
    to_version: str | None = None,
    cache_dir: Path | str | None = None,
    show_progress: bool = True,
    force: bool = False,
) -> Sec2PriMappingSet:
    """Return cross-release previous-symbol -> current-symbol Ensembl mappings.

    Ensembl's core schema has no previous-gene-symbol table, so genuine
    previous-symbol -> current-symbol transitions are derived by diffing
    each release's primary-label snapshot per stable ID. Delegates to
    :func:`pysec2pri.consolidate.build_label_history`, which walks every
    historical release (or a bounded range) and resumes on re-run. This is
    a network-heavy, run-on-demand operation, not part of normal
    single-release generation.

    Args:
        species: Canonical NCBI taxon ID (default ``9606`` for human).
        from_version: Optional lower bound (inclusive) on the release walk.
        to_version: Optional upper bound (inclusive) on the release walk.
        cache_dir: Directory for the resumable sidecar/cache. Defaults to
            :func:`pysec2pri.consolidate.default_cache_dir`.
        show_progress: Whether to show progress bars.
        force: Re-walk every release, ignoring resume state.

    Returns:
        LabelMappingSet of previous -> current symbol transitions.
    """
    from pysec2pri.consolidate import build_label_history

    return build_label_history(
        species=species,
        from_version=from_version,
        to_version=to_version,
        cache_dir=Path(cache_dir) if cache_dir else None,
        show_progress=show_progress,
        force=force,
    )


def generate_uniprot(
    input_path: Path | str | None = None,
    delac_file: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return UniProt secondary to primary accession mappings.

    Downloads sec_ac.txt and delac_sp.txt automatically when ``input_path``
    is omitted.

    Args:
        input_path: Local sec_ac.txt. Auto-downloaded if ``None``.
        delac_file: Local delac_sp.txt.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import UniProtParser

    release_date = None
    if input_path is None:
        files, release_date = _auto_download("uniprot", version)
        input_path = files.get("sec_ac") or next(iter(files.values()))
        if delac_file is None:
            delac_file = files.get("delac_sp")

    parser = UniProtParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse(
        Path(input_path),
        delac_path=Path(delac_file) if delac_file else None,
    )


def generate_hmdb(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return HMDB metabolite secondary to primary accession mappings.

    Downloads hmdb_metabolites.xml automatically when ``input_path`` is omitted.

    Args:
        input_path: Local hmdb_metabolites.xml (or .zip/.gz). Auto-downloaded if ``None``.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import HMDBMetaboliteParser

    release_date = None
    if input_path is None:
        files, release_date = _auto_download("hmdb_metabolites", version)
        input_path = files["metabolites"]

    parser = HMDBMetaboliteParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse(Path(input_path))


def generate_hmdb_proteins(
    input_path: Path | str | None = None,
    version: str | None = None,
    show_progress: bool = True,
) -> Sec2PriMappingSet:
    """Return HMDB protein secondary to primary accession mappings.

    Downloads hmdb_proteins.xml automatically when ``input_path`` is omitted.

    Args:
        input_path: Local hmdb_proteins.xml (or .zip/.gz). Auto-downloaded if ``None``.
        version: Version string for metadata.
        show_progress: Whether to show progress bars.
    """
    from pysec2pri.parsers import HMDBProteinParser

    release_date = None
    if input_path is None:
        files, release_date = _auto_download("hmdb_proteins", version)
        input_path = files["proteins"]

    parser = HMDBProteinParser(version=version, show_progress=show_progress)
    parser.release_date = release_date
    return parser.parse(Path(input_path))


def generate_wikidata(
    input_path: Path | str | None = None,
    entity_type: str | None = None,
    version: str | None = None,
    endpoint: str | None = None,
    show_progress: bool = True,
    test_subset: bool = False,
) -> Sec2PriMappingSet:
    """Return Wikidata redirect mappings via SPARQL (or a pre-downloaded TSV).

    Queries the QLever Wikidata endpoint when ``input_path`` is omitted.
    If ``entity_type`` is ``None``, all entity types (metabolites, genes,
    proteins) are queried and combined.

    Args:
        input_path: Pre-downloaded TSV file. Queries SPARQL if ``None``.
        entity_type: ``"metabolites"``, ``"chemicals"``, ``"genes"``, or
        `"proteins"``. Queries all types when ``None``.
        version: Version string for metadata (defaults to today's date).
        endpoint: Custom SPARQL endpoint URL.
        show_progress: Whether to show progress bars.
        test_subset: Use test queries limited to 10 results.
    """
    from pysec2pri.parsers import WikidataParser

    parser = WikidataParser(
        version=version,
        show_progress=show_progress,
        entity_type=entity_type or "metabolites",
        endpoint=endpoint,
        test_subset=test_subset,
    )

    if input_path is not None:
        return parser.parse_from_file(Path(input_path))

    if entity_type is None:
        return parser.parse_all()

    return parser.parse()


def generate_wikidata_labels(
    input_path: Path | str | None = None,
    entity_type: str | None = None,
    version: str | None = None,
    endpoint: str | None = None,
    show_progress: bool = True,
    test_subset: bool = False,
) -> LabelMappingSet:
    """Return Wikidata label mappings (previous label  to current label).

    Queries the QLever Wikidata endpoint when ``input_path`` is omitted.
    If ``entity_type`` is ``None``, all entity types are queried and
    their label mappings combined.

    Args:
        input_path: Pre-downloaded TSV file. Queries SPARQL if ``None``.
        entity_type: ``"metabolites"``, ``"chemicals"``, ``"genes"``, or
            ``"proteins"``. Queries all types when ``None``.
        version: Version string for metadata.
        endpoint: Custom SPARQL endpoint URL.
        show_progress: Whether to show progress bars.
        test_subset: Use test queries limited to 10 results.

    Returns:
        :class:`~pysec2pri.parsers.base.LabelMappingSet` with label mappings.
    """
    from pysec2pri.parsers import WikidataParser

    entity_types = ["metabolites", "genes", "proteins"] if entity_type is None else [entity_type]

    sets = []
    for etype in entity_types:
        parser = WikidataParser(
            version=version,
            show_progress=show_progress,
            entity_type=etype,
            endpoint=endpoint,
            test_subset=test_subset,
        )
        sets.append(parser.parse_labels(Path(input_path) if input_path else None))

    if len(sets) == 1:
        return sets[0]

    # Combine multiple LabelMappingSets into one
    all_mappings = [m for ms in sets for m in (ms.mappings or [])]
    combined = sets[0]
    combined.mappings = all_mappings
    combined.compute_cardinalities()
    return combined


def combine_mapping_sets(
    id_mappings: Sec2PriMappingSet | None,
    synonym_mappings: Sec2PriMappingSet | None,
) -> Sec2PriMappingSet:
    """Combine two mapping sets into one.

    Args:
        id_mappings: First mapping set (e.g. ID mappings).
        synonym_mappings: Second mapping set (e.g. synonym mappings).

    Returns:
        Combined mapping set.

    Raises:
        ValueError: If both mapping sets are ``None``.
    """
    if id_mappings is None and synonym_mappings is None:
        msg = "At least one mapping set must be provided"
        raise ValueError(msg)
    if id_mappings is None:
        return synonym_mappings  # type: ignore[return-value]
    if synonym_mappings is None:
        return id_mappings
    combined = list(id_mappings.mappings or [])
    combined.extend(synonym_mappings.mappings or [])
    id_mappings.mappings = combined
    return id_mappings


# Output helpers

_FORMAT_EXTENSIONS: dict[str, str] = {
    "rdf": ".ttl",
    "owl": "_owl.ttl",
    "json": ".json",
}


def _output_filename(base_name: str, fmt: str) -> Path:
    """Return the default filename for *base_name* in format *fmt*."""
    ext = _FORMAT_EXTENSIONS.get(fmt, ".tsv")
    if fmt == "owl":
        return Path(f"{base_name}{ext}")
    return Path(f"{base_name}_{fmt}{ext}")


def save(
    mapping_set: Sec2PriMappingSet,
    output_format: str,
    output: Path | str | None = None,
    *,
    base_name: str,
) -> Path:
    """Write *mapping_set* and return the path that was written.

    Delegates to :meth:`~pysec2pri.parsers.base.Sec2PriMappingSet.save` for
    single formats and :func:`write_all_formats` for ``"all"``.

    Args:
        mapping_set: The mapping set to write.
        output_format: One of ``sssom``, ``sec2pri``, ``pri_ids``,
            ``name2synonym``, ``label_sec2pri``, ``pri_labels``,
            ``rdf``, ``json``, ``owl``, or ``all``.
        output: Explicit output path or directory.  When ``None``, a
            default name derived from *base_name* is used.
        base_name: Stem used to derive file names, e.g. ``"hgnc_2026-04-07"``.

    Returns:
        The directory (for ``"all"``) or file path that was written.
    """
    out = Path(output) if output else None

    if output_format == "all":
        if out is None:
            out_dir = Path(base_name)
        elif out.suffix:
            out_dir = out.parent / base_name
        else:
            out_dir = out
        write_all_formats(mapping_set, out_dir, base_name)
        return out_dir

    # Resolve output path, then delegate to the mapping-set method
    if out is None:
        out_path = _output_filename(base_name, output_format)
    elif out.is_dir():
        out_path = out / _output_filename(base_name, output_format).name
    else:
        out_path = out

    return mapping_set.save(output_format, out_path)


def write_all_formats(
    mapping_set: Sec2PriMappingSet,
    output_dir: Path,
    base_name: str,
    include_name2synonym: bool = True,
) -> None:
    """Write mapping set in all output formats to a directory.

    Args:
        mapping_set: The mapping set to write.
        output_dir: Directory to write files to.
        base_name: Base name for output files (e.g., "chebi_3star_245").
        include_name2synonym: Whether to include name2synonym format.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    write_sssom(mapping_set, output_dir / f"{base_name}_sssom.tsv")

    if isinstance(mapping_set, IdMappingSet):
        write_sec2pri(mapping_set, output_dir / f"{base_name}_sec2pri.tsv")
        write_pri_ids(mapping_set, output_dir / f"{base_name}_pri_ids.txt")

    if isinstance(mapping_set, LabelMappingSet):
        write_label_sec2pri(mapping_set, output_dir / f"{base_name}_label_sec2pri.tsv")
        if include_name2synonym:
            write_name2synonym(mapping_set, output_dir / f"{base_name}_name2synonym.tsv")

    write_rdf(mapping_set, output_dir / f"{base_name}.ttl")
    write_json(mapping_set, output_dir / f"{base_name}.json")
    write_owl(mapping_set, output_dir / f"{base_name}_owl.ttl")


def write_diff_output(
    result: MappingDiff,
    output_path: Path,
) -> None:
    """Write diff results to a TSV file.

    Args:
        result: MappingDiff object with added/removed/changed mappings.
        output_path: Path to write the TSV file.
    """
    import polars as pl

    dfs = []

    if result.added_count > 0:
        added_df = result.added.with_columns(
            pl.lit("added").alias("change_type"),
            pl.lit(None).alias("old_subject_id"),
        ).select(
            [
                "change_type",
                "object_id",
                pl.col("subject_id").alias("new_subject_id"),
                "old_subject_id",
            ]
        )
        dfs.append(added_df)

    if result.removed_count > 0:
        removed_df = result.removed.with_columns(
            pl.lit("removed").alias("change_type"),
            pl.lit(None).alias("new_subject_id"),
        ).select(
            [
                "change_type",
                "object_id",
                "new_subject_id",
                pl.col("subject_id").alias("old_subject_id"),
            ]
        )
        dfs.append(removed_df)

    if result.changed_count > 0:
        changed_df = result.changed.with_columns(
            pl.lit("changed").alias("change_type"),
        ).select(
            [
                "change_type",
                "object_id",
                "new_subject_id",
                "old_subject_id",
            ]
        )
        dfs.append(changed_df)

    if dfs:
        combined = pl.concat(dfs)
        combined.write_csv(output_path, separator="\t")


def load_mapping(path: Path | str) -> IdMappingSet:
    """Load an ID mapping set from a pysec2pri TSV file.

    Accepts the ``sec2pri`` TSV format (columns ``subject_id``, ``object_id``,
    ``predicate_id``, ``mapping_cardinality``) and the full SSSOM TSV format
    (comment-prefixed metadata lines are skipped automatically).

    Args:
        path: Path to the TSV file to load.

    Returns:
        An :class:`~pysec2pri.parsers.base.IdMappingSet` populated from the
        file, ready to pass to :func:`resolve_ids`.
    """
    import pandas as pd
    from sssom_schema import Mapping

    path = Path(path)
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    ms = IdMappingSet(
        mapping_set_id=str(path),
        license="https://creativecommons.org/licenses/by/4.0/",
    )
    mappings: list[Mapping] = []
    for _, row in df.iterrows():
        m = Mapping(
            subject_id=row.get("subject_id") or "",
            object_id=row.get("object_id") or "",
            predicate_id=row.get("predicate_id") or "",
            mapping_justification=row.get("mapping_justification")
            or "semapv:BackgroundKnowledgeBasedMatching",
            mapping_cardinality=row.get("mapping_cardinality") or None,
        )
        mappings.append(m)
    ms.mappings = mappings
    return ms


def load_label_mapping(path: Path | str) -> LabelMappingSet:
    """Load a label/label mapping set from a pysec2pri TSV file.

    Accepts two column-name conventions:

    - **New** (``label_sec2pri`` tabular output): ``secondary_id``,
      ``secondary_label``, ``primary_id``, ``primary_label``,
      ``predicate_id``, ``mapping_cardinality``.
    - **Legacy** (SSSOM or old tabular output): ``subject_id``,
      ``subject_label``, ``object_id``, ``object_label``, ``predicate_id``.

    Full SSSOM TSV (comment-prefixed metadata lines) is also accepted.

    Args:
        path: Path to the TSV file to load.

    Returns:
        A :class:`~pysec2pri.parsers.base.LabelMappingSet` populated from
        the file, ready to pass to :func:`resolve_labels`.
    """
    import pandas as pd
    from sssom_schema import Mapping

    path = Path(path)
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    ms = LabelMappingSet(
        mapping_set_id=str(path),
        license="https://creativecommons.org/licenses/by/4.0/",
    )

    def _col(row: pd.Series, *names: str) -> str | None:
        for name in names:
            val = row.get(name)
            if val and str(val).strip():
                return str(val).strip()
        return None

    mappings: list[Mapping] = []
    for _, row in df.iterrows():
        m = Mapping(
            subject_id=_col(row, "secondary_id", "subject_id") or "",
            subject_label=_col(row, "secondary_label", "subject_label"),
            object_id=_col(row, "primary_id", "object_id") or "",
            object_label=_col(row, "primary_label", "object_label"),
            predicate_id=_col(row, "predicate_id") or "",
            mapping_justification=_col(row, "mapping_justification")
            or "semapv:BackgroundKnowledgeBasedMatching",
            mapping_cardinality=_col(row, "mapping_cardinality"),
        )
        mappings.append(m)
    ms.mappings = mappings
    return ms


def resolve_ids(
    input_path: Path | str | list[str],
    mapping_set: Sec2PriMappingSet,
    at: str | list[str] | None = None,
    *,
    output_path: Path | str | None = None,
    suffix: str = "_primary",
    sep: str | None = None,
    synonyms: str | None = None,
    label_mapping_set: Sec2PriMappingSet | None = None,
    xref: str | None = None,
    xref_mapping: XrefMapping | None = None,
    report_path: Path | str | None = None,
) -> pd.DataFrame | str | list[str]:
    r"""Resolve secondary IDs to primary IDs.

    Direct lookup: when *input_path* is a plain identifier string or a list
    of identifier strings (i.e. not a path to an existing file), the function
    returns the resolved primary ID(s).  *at*, *output_path*, *suffix*, and
    *sep* are ignored in this mode::

        resolve_ids("HMDB00001", hmdb_ms)  # -> "HMDB:HMDB0000001"
        resolve_ids(["HMDB00001", "HMDB00002"], hmdb_ms)  # -> ["...", "..."]

    DataFrame mode: when *input_path* points to an existing TSV/CSV
    file, *at* is required.  The file is read with
    ``pandas.read_csv`` and for each column named in *at* a new column
    ``<col><suffix>`` is appended containing the resolved primary IDs.
    Identifiers not present in *mapping_set* are kept unchanged.

    Args:
        input_path: An identifier string, a list of identifier strings, or
            the path to a TSV/CSV file.
        mapping_set: A :class:`~pysec2pri.parsers.base.Sec2PriMappingSet`
            (e.g. the result of ``generate_hgnc()``).
        at: Column name(s) to resolve.  Required in DataFrame mode;
            ignored in direct-lookup mode.
        output_path: If given, the resulting DataFrame is written to this
            path (DataFrame mode only).
        suffix: Suffix appended to each resolved column name
            (default ``"_primary"``).
        sep: Delimiter for reading the file.  Inferred from the extension
            when ``None`` (``"\\t"`` for ``.tsv``, ``","`` otherwise).
        xref: *DataFrame mode only.* Column with a per-row cross-reference
            token, passed through to :func:`~pysec2pri.update_ids.update_ids`.
        xref_mapping: The :class:`~pysec2pri.context.XrefMapping` crosswalk
            table to resolve *xref* tokens against. Required when *xref* is
            given.
        report_path: When given, every disambiguation attempt (from
            *synonyms* and/or *xref*) is logged to this TSV.

    Returns:
        A resolved identifier string, a list of resolved strings (direct-lookup
        mode), or a :class:`pandas.DataFrame` with one additional column per
        entry in *at* (DataFrame mode).
    """
    import pandas as pd

    from pysec2pri.update_ids import build_lookup, update_ids

    # list direct-lookup mode
    if isinstance(input_path, list):
        lkp = build_lookup(mapping_set)
        return [lkp.get(v, v) for v in input_path]

    # single-value lookup mode
    input_path_obj = Path(input_path)
    if not input_path_obj.exists():
        lkp = build_lookup(mapping_set)
        return lkp.get(str(input_path), str(input_path))

    # DataFrame mode
    if at is None:
        raise TypeError("resolve_ids() requires 'at' when input_path is a file")

    if sep is None:
        sep = "\t" if input_path_obj.suffix.lower() == ".tsv" else ","

    df = pd.read_csv(input_path_obj, sep=sep, dtype=str)
    lkp = build_lookup(mapping_set)
    result = update_ids(
        df,
        mapping_set,
        at=at,
        suffix=suffix,
        lookup=lkp,
        synonyms=synonyms,
        label_mapping_set=label_mapping_set,
        xref=xref,
        xref_mapping=xref_mapping,
        report_path=report_path,
    )

    if output_path is not None:
        output_path = Path(output_path)
        out_sep = "\t" if output_path.suffix.lower() == ".tsv" else ","
        result.to_csv(output_path, sep=out_sep, index=False)

    return result


def resolve_labels(
    input_path: Path | str | list[str],
    mapping_set: Sec2PriMappingSet,
    at: str | list[str] | None = None,
    *,
    output_path: Path | str | None = None,
    suffix: str = "_current",
    sep: str | None = None,
    synonyms: str | None = None,
    xref: str | None = None,
    xref_mapping: XrefMapping | None = None,
    report_path: Path | str | None = None,
) -> pd.DataFrame | str | list[str]:
    r"""Resolve previous/alias labels to current labels.

    Direct lookup: when *input_path* is a plain label string or a list
    of label strings (i.e. not a path to an existing file), the function
    returns the resolved current label(s).  *at*, *output_path*, *suffix*,
    and *sep* are ignored in this mode::

        resolve_labels("Ibuprofen", chebi_ms)  # -> "ibuprofen"
        resolve_labels(["Ibuprofen", "Glucose"], chebi_ms)  # -> ["...", "..."]

    DataFrame mode: when *input_path* points to an existing TSV/CSV
    file, *at* is required.  For each column named in *at* a new
    column ``<col><suffix>`` is appended containing the resolved current
    labels.  Symbols not present in *mapping_set* are kept unchanged.

    Args:
        input_path: A label string, a list of label strings, or the path
            to a TSV/CSV file.
        mapping_set: A :class:`~pysec2pri.parsers.base.LabelMappingSet`
            (e.g. the result of ``generate_hgnc_labels()``).
        at: Column name(s) to resolve.  Required in DataFrame mode;
            ignored in direct-lookup mode.
        output_path: If given, the resulting DataFrame is written to this
            path (DataFrame mode only).
        suffix: Suffix appended to each resolved column name
            (default ``"_current"``).
        sep: Delimiter for reading the file.  Inferred from the extension
            when ``None`` (``"\\t"`` for ``.tsv``, ``","`` otherwise).
        xref: *DataFrame mode only.* Column with a per-row cross-reference
            token, passed through to
            :func:`~pysec2pri.update_ids.update_labels`.
        xref_mapping: The :class:`~pysec2pri.context.XrefMapping` crosswalk
            table to resolve *xref* tokens against. Required when *xref* is
            given.
        report_path: When given, every disambiguation attempt (from
            *synonyms* and/or *xref*) is logged to this TSV.

    Returns:
        A resolved label string, a list of resolved strings (direct-lookup
        mode), or a :class:`pandas.DataFrame` with one additional column per
        entry in *at* (DataFrame mode).
    """
    import pandas as pd

    from pysec2pri.update_ids import build_label_lookup, update_labels

    # list direct-lookup mode
    if isinstance(input_path, list):
        lkp = build_label_lookup(mapping_set)
        return [lkp.get(v, v) for v in input_path]

    # single-value direct-lookup mode
    input_path_obj = Path(input_path)
    if not input_path_obj.exists():
        lkp = build_label_lookup(mapping_set)
        return lkp.get(str(input_path), str(input_path))

    # DataFrame mode
    if at is None:
        raise TypeError("resolve_labels() requires 'at' when input_path is a file")

    if sep is None:
        sep = "\t" if input_path_obj.suffix.lower() == ".tsv" else ","

    df = pd.read_csv(input_path_obj, sep=sep, dtype=str)
    lkp = build_label_lookup(mapping_set)
    result = pd.DataFrame(
        update_labels(
            df,
            mapping_set,
            at=at,
            suffix=suffix,
            lookup=lkp,
            synonyms=synonyms,
            xref=xref,
            xref_mapping=xref_mapping,
            report_path=report_path,
        )
    )

    if output_path is not None:
        output_path = Path(output_path)
        out_sep = "\t" if output_path.suffix.lower() == ".tsv" else ","
        result.to_csv(output_path, sep=out_sep, index=False)

    return result


def list_versions(datasource: str) -> Any:
    """List all available archive versions for a datasource.

    For datasources that publish versioned archives (ChEBI, HGNC, UniProt),
    this queries the remote archive index and returns all available version
    strings sorted in ascending order.

    NCBI and HMDB do not maintain versioned archives; calling this function
    for those datasources raises :class:`ValueError`.

    Args:
        datasource: Datasource name, one of ``"chebi"``, ``"hgnc"``, or
            ``"uniprot"``.

    Returns:
        Sorted list of version strings.  Format depends on the datasource:

        - **chebi**: integer release numbers, e.g. ``["200", ..., "245"]``
        - **hgnc**: ISO dates, e.g. ``["2023-01-01", ..., "2026-04-07"]``
        - **uniprot**: release IDs, e.g. ``["2024_01", "2024_02", ...]``

    Raises:
        ValueError: If *datasource* is unknown or has no versioned archive.
    """
    from pysec2pri.download import list_versions as _list_versions

    return _list_versions(datasource)


def _crosswalk_via_xref(
    tokens: list[str],
    xref_mapping: XrefMapping,
    to: str,
    decisions: list[DecisionRecord],
) -> dict[str, str]:
    """Resolve each of *tokens* through a flat ``subject_id -> object_*`` crosswalk.

    Unlike the secondary/primary ambiguity handled by
    :mod:`pysec2pri.update_ids`, a crosswalk table can itself carry more than
    one distinct target for the same token; that case is also left
    unresolved (and logged) rather than guessed.
    """
    from pysec2pri.context import DecisionRecord

    index = xref_mapping.by_subject()
    result: dict[str, str] = {}
    for tok in tokens:
        if tok in result:
            continue
        candidates = index.get(tok) or []
        if not candidates:
            result[tok] = ""
            decisions.append(
                DecisionRecord("xref_crosswalk", tok, None, None, False, "no crossreference entry")
            )
            continue
        targets = {(c.object_id, c.object_label or "") for c in candidates}
        if len(targets) > 1:
            result[tok] = ""
            decisions.append(
                DecisionRecord(
                    "xref_crosswalk",
                    tok,
                    None,
                    None,
                    False,
                    f"ambiguous crosswalk: {len(targets)} distinct targets",
                )
            )
            continue
        record = candidates[0]
        value = record.object_id if to == "hgnc_id" else (record.object_label or "")
        result[tok] = value
        decisions.append(
            DecisionRecord(
                "xref_crosswalk", tok, record.predicate_id, value, True, "unique crosswalk match"
            )
        )
    return result


def _crosswalk_symbol(
    input_data: str | list[str] | pd.DataFrame,
    *,
    to: str,
    at: str | None,
    version: str | None,
    label_mapping_set: Sec2PriMappingSet | None,
    report_path: Path | str | None,
    show_progress: bool,
) -> dict[str, str] | pd.DataFrame:
    """``frm == "symbol"`` branch of :func:`crosswalk`."""
    from pysec2pri.context import DecisionRecord, write_decision_log
    from pysec2pri.update_ids import (
        build_ambiguous_labels_set,
        build_primary_token_to_id,
        update_labels,
    )

    label_ms = label_mapping_set or generate_hgnc_labels(
        version=version, show_progress=show_progress
    )
    token_to_id = build_primary_token_to_id(label_ms)
    ambiguous_symbols = build_ambiguous_labels_set(label_ms)

    def _log_ambiguous(values: set[str]) -> None:
        if report_path is None:
            return
        decisions = [
            DecisionRecord("symbol", v, None, None, False, "ambiguous symbol, no context provided")
            for v in sorted(values & ambiguous_symbols)
        ]
        write_decision_log(decisions, report_path)

    if isinstance(input_data, (str, list)):
        resolved = update_labels(input_data, label_ms)
        _log_ambiguous(set(resolved))
        if to == "symbol":
            return resolved
        return {k: (token_to_id.get(v, "") if v else "") for k, v in resolved.items()}

    if at is None:
        raise TypeError("crosswalk() requires 'at' when input_data is a DataFrame")
    result = update_labels(input_data, label_ms, at=at)
    result[f"{at}_symbol"] = result[f"{at}_current"]
    result[f"{at}_hgnc_id"] = result[f"{at}_symbol"].map(
        lambda v: token_to_id.get(v, "") if v else ""
    )
    _log_ambiguous(set(input_data[at].astype(str)))
    return result


def _resolve_crosswalk_xref_mapping(
    frm: str,
    xref_mapping: XrefMapping | None,
    xref_source: str,
    show_progress: bool,
) -> XrefMapping:
    """Return *xref_mapping*, or download the configured *xref_source* for it."""
    from pysec2pri.context import download_xref_source
    from pysec2pri.parsers.base import get_datasource_config

    if xref_mapping is not None:
        return xref_mapping

    cfg = get_datasource_config("hgnc")
    src = cfg.xref_source(xref_source)
    if src is None:
        known = ", ".join(s.id for s in cfg.xref_sources) or "(none configured)"
        raise ValueError(f"Unknown xref_source {xref_source!r}. Known: {known}")
    subject_col = src.subject_id_cols.get(frm)
    if subject_col is None:
        known_keys = ", ".join(sorted(src.subject_id_cols)) or "(none)"
        raise ValueError(
            f"Unsupported 'frm' vocabulary {frm!r} for xref_source {xref_source!r}. "
            f"Known: {known_keys}"
        )
    return download_xref_source(src, subject_col, show_progress=show_progress)


def _crosswalk_xref(
    input_data: str | list[str] | pd.DataFrame,
    *,
    frm: str,
    to: str,
    at: str | None,
    xref_mapping: XrefMapping | None,
    xref_source: str,
    report_path: Path | str | None,
    show_progress: bool,
) -> dict[str, str] | pd.DataFrame:
    """Non-``"symbol"`` ``frm`` branch of :func:`crosswalk`: a flat crosswalk-table lookup."""
    from pysec2pri.context import write_decision_log

    mapping = _resolve_crosswalk_xref_mapping(frm, xref_mapping, xref_source, show_progress)
    decisions: list[DecisionRecord] = []

    if isinstance(input_data, str):
        result_dict = _crosswalk_via_xref([input_data], mapping, to, decisions)
    elif isinstance(input_data, list):
        result_dict = _crosswalk_via_xref(list(dict.fromkeys(input_data)), mapping, to, decisions)
    else:
        if at is None:
            raise TypeError("crosswalk() requires 'at' when input_data is a DataFrame")
        if at not in input_data.columns:
            raise KeyError(f"Column {at!r} not found in DataFrame.")
        unique_tokens = [str(v) for v in input_data[at].dropna().unique()]
        lookup = _crosswalk_via_xref(unique_tokens, mapping, to, decisions)
        result = input_data.copy()
        result[f"{at}_{to}"] = result[at].astype(str).map(lambda v: lookup.get(v, ""))
        if report_path is not None:
            write_decision_log(decisions, report_path)
        return result

    if report_path is not None:
        write_decision_log(decisions, report_path)
    return result_dict


def crosswalk(
    input_data: str | list[str] | pd.DataFrame,
    *,
    frm: str,
    to: str = "hgnc_id",
    at: str | None = None,
    version: str | None = None,
    label_mapping_set: Sec2PriMappingSet | None = None,
    xref_mapping: XrefMapping | None = None,
    xref_source: str = "hgnc_custom",
    report_path: Path | str | None = None,
    show_progress: bool = True,
) -> dict[str, str] | pd.DataFrame:
    r"""Map a gene identifier from one vocabulary to another, via HGNC.

    A thin convenience wrapper, not a separate resolver: ``frm="symbol"``
    resolves through :func:`generate_hgnc_labels` +
    :func:`~pysec2pri.update_ids.update_labels`, so a *previous* HGNC symbol
    still resolves to its current identity (the temporal aspect) and a
    genuinely ambiguous symbol is left blank and reported rather than
    guessed. ``frm`` in ``"ensembl"``, ``"entrez"``, ``"refseq"``, or
    ``"uniprot"`` resolves through HGNC's own cross-reference crosswalk
    table, downloaded from the ``xref_sources`` declared in
    ``config/hgnc.yaml`` unless *xref_mapping* is supplied explicitly.

    Args:
        input_data: A single identifier string, a list of identifier
            strings, or a :class:`pandas.DataFrame` (requires *at*).
        frm: Source vocabulary: ``"symbol"``, ``"ensembl"``, ``"entrez"``,
            ``"refseq"``, or ``"uniprot"``.
        to: Target vocabulary: ``"hgnc_id"`` (default) or ``"symbol"``.
        at: *DataFrame mode only.* Column containing the *frm* values.
        version: HGNC release version to resolve against (``None`` = latest).
            Ignored when *label_mapping_set* is given.
        label_mapping_set: Pre-built :class:`~pysec2pri.parsers.base.LabelMappingSet`
            for ``frm == "symbol"`` (e.g. the result of
            ``generate_hgnc_labels()``), to avoid re-downloading. When
            ``None``, it is generated automatically.
        xref_mapping: Pre-built crosswalk table for non-``"symbol"`` *frm*
            values. When ``None``, it is downloaded automatically (see
            *xref_source*); ignored when ``frm == "symbol"``.
        xref_source: Which ``xref_sources`` entry declared in
            ``config/hgnc.yaml`` to download when *xref_mapping* is not
            given (default ``"hgnc_custom"``).
        report_path: When given, every disambiguation attempt is logged as a
            TSV (see :func:`pysec2pri.context.write_decision_log`).
        show_progress: Forwarded to the HGNC generator / downloader.

    Returns:
        dict[str, str]: when *input_data* is a ``str``/``list[str]``, mapping
        each unique input value to its resolved target (``""`` when
        unresolved/ambiguous).

        pandas.DataFrame: when *input_data* is a ``DataFrame``, a copy with
        the columns :func:`~pysec2pri.update_ids.update_labels` adds
        (``<at>_current``/``<at>_current_id``) for ``frm="symbol"``, plus a
        uniform ``<at>_symbol``/``<at>_hgnc_id`` pair; for other *frm*
        values, a single ``<at>_<to>`` column.
    """
    if to not in ("hgnc_id", "symbol"):
        raise ValueError(f"Unsupported 'to' vocabulary: {to!r}. Choose 'hgnc_id' or 'symbol'.")

    if frm == "symbol":
        return _crosswalk_symbol(
            input_data,
            to=to,
            at=at,
            version=version,
            label_mapping_set=label_mapping_set,
            report_path=report_path,
            show_progress=show_progress,
        )

    return _crosswalk_xref(
        input_data,
        frm=frm,
        to=to,
        at=at,
        xref_mapping=xref_mapping,
        xref_source=xref_source,
        report_path=report_path,
        show_progress=show_progress,
    )


def find_ambiguous(
    mapping_set: Sec2PriMappingSet,
) -> AmbiguousMappingSet:
    """Find identifiers that are ambiguous in *mapping_set*.

    An identifier is ambiguous when it appears both as a ``subject_id`` (i.e. a
    secondary/previous term) and as a current primary identifier.  Such
    entries cannot be automatically resolved without risk of corrupting
    references that are already current.

    This is a convenience wrapper around
    :meth:`~pysec2pri.parsers.base.Sec2PriMappingSet.find_ambiguous`.

    Args:
        mapping_set: A :class:`~pysec2pri.parsers.base.Sec2PriMappingSet`
            (e.g. the result of ``generate_hgnc()``).

    Returns:
        An :class:`~pysec2pri.parsers.base.AmbiguousMappingSet` whose
        ``mappings`` list contains one entry for each conflicting subject, with a
        ``comment`` explaining the conflict.
    """
    return mapping_set.find_ambiguous()


# The functions here use old names, remove at some point

generate_chebi_ids = generate_chebi
generate_chebi_labels = generate_chebi_synonyms
generate_ensembl_ids = generate_ensembl
generate_hgnc_ids = generate_hgnc
generate_hgnc_labels = generate_hgnc_labels
generate_ncbi_ids = generate_ncbi
generate_ncbi_labels = generate_ncbi_labels
generate_hmdb_ids = generate_hmdb
generate_hmdb_proteins_ids = generate_hmdb_proteins
generate_uniprot_ids = generate_uniprot
generate_vgnc_ids = generate_vgnc
generate_wikidata_ids = generate_wikidata
generate_wikidata_labels = generate_wikidata_labels
