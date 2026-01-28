#######
 Usage
#######

Quick Start
===========

The easiest way to generate mappings is with auto-download:

.. code-block:: bash

    # Download latest ChEBI and generate SSSOM mappings
    pysec2pri chebi -o chebi_sec2pri.sssom.tsv

    # Same for other databases
    pysec2pri hmdb -o hmdb_sec2pri.sssom.tsv
    pysec2pri hgnc -o hgnc_sec2pri.sssom.tsv
    pysec2pri ncbi -o ncbi_sec2pri.sssom.tsv
    pysec2pri uniprot -o uniprot_sec2pri.sssom.tsv

Python API
==========

Download and Parse
------------------

The recommended workflow is to download files first, then parse:

.. code-block:: python

    from pysec2pri import parse_chebi, write_sssom
    from pysec2pri.download import download_datasource

    # Download latest files
    files = download_datasource("chebi", output_dir="./data")

    # Parse and write output
    mapping_set = parse_chebi(files["sdf"])
    write_sssom(mapping_set, "chebi_sec2pri.sssom.tsv")

Parse Local Files
-----------------

If you already have local files:

.. code-block:: python

    from pysec2pri import (
        parse_chebi,
        parse_hmdb,
        parse_hgnc,
        parse_ncbi,
        parse_uniprot,
        write_sssom,
    )

    # Parse ChEBI SDF file
    chebi_mappings = parse_chebi("ChEBI_complete_3star.sdf")
    write_sssom(chebi_mappings, "chebi_sec2pri.sssom.tsv")

    # Parse HMDB XML (ZIP or directory)
    hmdb_mappings = parse_hmdb("hmdb_metabolites.zip")
    write_sssom(hmdb_mappings, "hmdb_sec2pri.sssom.tsv")

    # Parse HGNC files
    hgnc_mappings = parse_hgnc(
        "withdrawn.txt",
        complete_set_file="hgnc_complete_set.txt"
    )
    write_sssom(hgnc_mappings, "hgnc_sec2pri.sssom.tsv")

    # Parse NCBI Gene files
    ncbi_mappings = parse_ncbi("gene_history.gz", "gene_info.gz")
    write_sssom(ncbi_mappings, "ncbi_sec2pri.sssom.tsv")

    # Parse UniProt files
    uniprot_mappings = parse_uniprot("sec_ac.txt", "delac_sp.txt")
    write_sssom(uniprot_mappings, "uniprot_sec2pri.sssom.tsv")

Compare Releases
----------------

Use the diff module to compare mapping sets between releases:

.. code-block:: python

    from pysec2pri import parse_chebi
    from pysec2pri.diff import diff_mapping_sets, summarize_diff

    old_set = parse_chebi("chebi_v220.sdf")
    new_set = parse_chebi("chebi_v221.sdf")

    diff = diff_mapping_sets(old_set, new_set)
    print(summarize_diff(diff))

    # Access detailed changes
    print(f"Added: {diff.added_count}")
    print(f"Removed: {diff.removed_count}")

Check for Updates
-----------------

Check if a new release is available:

.. code-block:: python

    from pysec2pri.download import check_release, get_latest_release_info

    # Get latest release info
    info = get_latest_release_info("chebi")
    print(f"Latest version: {info.version}")
    print(f"Release date: {info.release_date}")

    # Check against current version
    info = check_release("chebi", current_version="220")
    if info.is_new:
        print("New release available!")

Working with MappingSet
-----------------------

The ``MappingSet`` object provides access to individual mappings:

.. code-block:: python

    mapping_set = parse_chebi("ChEBI_complete_3star.sdf")

    print(f"Source: {mapping_set.datasource_name}")
    print(f"Total mappings: {len(mapping_set)}")

    for mapping in mapping_set.iter_mappings():
        print(f"{mapping.secondary_id} -> {mapping.primary_id}")

Command Line Interface
======================

See the :doc:`cli` documentation for command-line usage.
