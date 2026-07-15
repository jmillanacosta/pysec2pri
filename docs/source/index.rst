###################################
 pySec2Pri |release| Documentation
###################################

**pysec2pri** maps secondary (retired/withdrawn) biological database identifiers
and labels to primary (current) ones, with
`SSSOM <https://mapping-commons.github.io/sssom/>`_ output by default.

Supported databases
===================

.. list-table::
   :header-rows: 1

   * - Database
     - ``source``
     - Mapping sets
   * - ChEBI
     - ``chebi``
     - ``ids``, ``labels``
   * - Ensembl
     - ``ensembl``
     - ``ids``, ``labels``
   * - HGNC
     - ``hgnc``
     - ``ids``, ``labels``
   * - HMDB Metabolites
     - ``hmdb_metabolites``
     - ``ids``
   * - HMDB Proteins
     - ``hmdb_proteins``
     - ``ids``
   * - NCBI Gene
     - ``ncbi``
     - ``ids``, ``labels``
   * - UniProt
     - ``uniprot``
     - ``ids``
   * - VGNC
     - ``vgnc``
     - ``ids``, ``labels``
   * - Wikidata
     - ``wikidata``
     - ``ids``, ``labels``

:func:`~pysec2pri.api.sources` returns this list at run time. To add another,
see :doc:`adding_a_source`: it takes a config and a parser.

Quick Start
===========

**Generate a mapping set (CLI):**

.. code-block:: bash

    pysec2pri hgnc ids
    pysec2pri chebi labels

**Update IDs or labels in a file (CLI):**

.. code-block:: bash

    pysec2pri update-ids data.tsv hgnc --at gene_id -o data_primary.tsv
    pysec2pri update-labels data.tsv hgnc --at label
    # reuse a saved mapping file
    pysec2pri update-ids data.tsv hgnc --at gene_id --mapping hgnc_sssom.tsv

**Python API:**

.. code-block:: python

    from pysec2pri import generate_ids, generate_labels, resolve_ids, resolve_labels
    from pysec2pri import load_mapping, load_label_mapping, sources

    sources()                                # every source

    ms = generate_ids("hgnc")
    resolve_ids("HGNC:131", ms)              # : "HGNC:145"
    resolve_ids(["HGNC:131", "HGNC:2"], ms)  # : ["HGNC:145", ...]

    lms = generate_labels("hgnc")
    resolve_labels("BRCA1_OLD", lms)         # : "BRCA1"

    # read back a saved SSSOM file
    ms = load_mapping("hgnc_sssom.tsv")
    lms = load_label_mapping("hgnc_labels_sssom.tsv")

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   installation
   cli

.. toctree::
   :maxdepth: 1
   :caption: Features

   generate
   update_ids
   exports
   gtex

.. toctree::
   :maxdepth: 1
   :caption: Extending

   adding_a_source

.. toctree::
   :maxdepth: 1
   :caption: API Reference

   api
   parsers
   models
   download
