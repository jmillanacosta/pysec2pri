########################
 Command Line Interface
########################

pysec2pri automatically installs the command ``pysec2pri``. See ``pysec2pri --help`` for
usage details.

Overview
========

All commands **auto-download source files by default** from upstream databases.
Use input options to specify local files instead.

Quick Examples
--------------

.. code-block:: bash

    # Auto-download and generate ChEBI mappings
    pysec2pri chebi -o chebi_sec2pri.sssom.tsv

    # Use a local file instead
    pysec2pri chebi --input ChEBI_complete_3star.sdf -o chebi.sssom.tsv

    # Download files only (without processing)
    pysec2pri download chebi -o ./data

    # Check for new releases
    pysec2pri check-release all

    # Compare two SSSOM files
    pysec2pri diff old.sssom.tsv new.sssom.tsv

Common Options
--------------

All datasource commands support these common options:

``-o, --output``
    Output SSSOM TSV file path.

``--output-dir``
    Directory for downloads and output (default: current directory).

``-v, --version``
    Version/release string to include in output metadata.

``-d, --date``
    Mapping date (YYYY-MM-DD format).

``--keep-download``
    Keep downloaded files after processing (by default, they are cleaned up).

Full Command Reference
======================

.. click:: pysec2pri.cli:main
    :prog: pysec2pri
    :show-nested:
