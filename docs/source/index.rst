###################################
 pySec2Pri |release| Documentation
###################################

**pysec2pri** maps secondary (retired/withdrawn) biological database identifiers
and symbols to primary (current) ones, with
`SSSOM <https://mapping-commons.github.io/sssom/>`_ output by default.

Supported databases: ChEBI, HMDB, HGNC, NCBI Gene, UniProt, Wikidata.

Quick Start
===========

**Generate a mapping set (CLI):**

.. code-block:: bash

    pysec2pri hgnc ids
    pysec2pri chebi synonyms

**Update IDs or symbols in a file (CLI):**

.. code-block:: bash

    pysec2pri update-ids data.tsv hgnc --at gene_id -o data_primary.tsv
    pysec2pri update-symbols data.tsv hgnc --at symbol
    # reuse a saved mapping file
    pysec2pri update-ids data.tsv hgnc --at gene_id --mapping hgnc_sssom.tsv

**Python API:**

.. code-block:: python

    from pysec2pri import generate_hgnc, resolve_ids
    from pysec2pri import generate_hgnc_symbols, resolve_symbols
    from pysec2pri import load_mapping, load_label_mapping

    ms = generate_hgnc()
    resolve_ids("HGNC:131", ms)              # → "HGNC:145"
    resolve_ids(["HGNC:131", "HGNC:2"], ms)  # → ["HGNC:145", ...]

    lms = generate_hgnc_symbols()
    resolve_symbols("BRCA1_OLD", lms)        # → "BRCA1"

    # load from a saved sec2pri / SSSOM file
    ms = load_mapping("hgnc_sec2pri.tsv")
    lms = load_label_mapping("hgnc_symbol2prev.tsv")

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   installation
   cli

.. toctree::
   :maxdepth: 1
   :caption: Features

   download
   exports
   diff

.. toctree::
   :maxdepth: 1
   :caption: API Reference

   api
   parsers
   models
   update_ids
