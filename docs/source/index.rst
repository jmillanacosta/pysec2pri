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



.. code-block:: bash

    pysec2pri chebi

This will automatically download and process the latest ChEBI release, generating a default SSSOM mapping file in the current directory.

You can also process a local file:

.. code-block:: bash

    pysec2pri chebi ChEBI_complete_3star.sdf

This will generate a default SSSOM mapping file from your local SDF file.

Other possibilities:

- **Specify output file:**

    .. code-block:: bash

            pysec2pri chebi ChEBI_complete_3star.sdf --output my_mappings.sssom.tsv

    Use this to choose a custom output filename.

- **Specifying outputs:**

    .. code-block:: bash

.. toctree::
    :maxdepth: 2
    :caption: Getting Started
    :name: start

    installation
    cli

.. toctree::
    :maxdepth: 2
    :caption: Helpers
    :name: features

    download
    diff

.. toctree::
    :maxdepth: 2
    :caption: API Reference
    :name: api

    api
    models
    constants
    parsers
