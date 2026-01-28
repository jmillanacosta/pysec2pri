###################################
 pySec2Pri |release| Documentation
###################################

**pysec2pri** converts secondary (retired/withdrawn) biological database identifiers
to primary (current) identifiers. It generates mappings in SSSOM format for use with
BridgeDb and other identifier mapping tools.

Supported databases:

- **ChEBI** - Chemical Entities of Biological Interest
- **HMDB** - Human Metabolome Database
- **HGNC** - HUGO Gene Nomenclature Committee
- **NCBI Gene** - Entrez Gene
- **UniProt** - Protein sequence database
- **Wikidata** - Redirect mappings from Wikidata knowledge graph

Quick Start
===========

.. code-block:: python

    from pysec2pri import parse_chebi, write_sssom

    mapping_set = parse_chebi("ChEBI_complete_3star.sdf")
    write_sssom(mapping_set, "chebi_sec2pri.sssom.tsv")

Or via CLI:

.. code-block:: bash

    pysec2pri chebi ChEBI_complete_3star.sdf -o chebi_sec2pri.sssom.tsv

.. toctree::
    :maxdepth: 2
    :caption: Getting Started
    :name: start

    installation
    usage
    cli

.. toctree::
    :maxdepth: 2
    :caption: Features
    :name: features

    download
    diff
    wikidata

.. toctree::
    :maxdepth: 2
    :caption: API Reference
    :name: api

    api
    models
    constants

********************
 Indices and Tables
********************

- :ref:`genindex`
- :ref:`modindex`
- :ref:`search`
