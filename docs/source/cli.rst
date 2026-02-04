######################
Command Line Interface
######################

The ``pysec2pri`` CLI provides commands for each supported database.

All commands output SSSOM-compliant TSV files by default.

Usage
=====

.. code-block:: bash

    pysec2pri [OPTIONS] COMMAND [ARGS]...

Commands
========

.. click:: pysec2pri.cli:main
    :prog: pysec2pri
    :nested: full

Output Formats
==============

- ``sssom`` (default): Full SSSOM TSV with metadata header
- ``sec2pri``: Simple TSV with subject_id, object_id, predicate_id
- ``name2synonym``: TSV with subject_label, object_label (ChEBI only)
- ``symbol2prev``: TSV with gene symbol mappings (HGNC/NCBI only)
