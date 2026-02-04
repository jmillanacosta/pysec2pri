###################################
 pySec2Pri |release| Documentation
###################################

**pysec2pri** maps secondary (retired/withdrawn) biological database identifiers
to primary (current) identifiers using the `SSSOM <https://mapping-commons.github.io/sssom/>`_ standard.

SSSOM Interoperability
======================

SSSOM outputs use **sssom-schema**:

- Mappings use ``sssom_schema.Mapping`` objects
- Mapping sets extend ``sssom_schema.MappingSet``
- Exports to SSSOM TSV via the ``sssom`` Python library

Supported Databases
===================

+--------------+------------------------------------------+
| Database     | Mapping Types                            |
+==============+==========================================+
| **ChEBI**    | Secondary→Primary IDs, Name→Synonyms     |
+--------------+------------------------------------------+
| **HMDB**     | Secondary→Primary IDs                    |
+--------------+------------------------------------------+
| **HGNC**     | Withdrawn→Current IDs, Symbol→Previous   |
+--------------+------------------------------------------+
| **NCBI Gene**| Discontinued→Current IDs, Symbol→Aliases |
+--------------+------------------------------------------+
| **UniProt**  | Secondary→Primary accessions             |
+--------------+------------------------------------------+
| **Wikidata** | Redirect mappings (SPARQL)               |
+--------------+------------------------------------------+

Quick Start
===========

**CLI:**

.. code-block:: bash

    # Parse ChEBI and output SSSOM
    pysec2pri chebi chebi_3star.sdf -o chebi.sssom.tsv

    # Parse HGNC with custom format
    pysec2pri hgnc hgnc_withdrawn.tsv --format sec2pri

**Python API:**

.. code-block:: python

    from pysec2pri.api import parse_chebi, write_sssom

    mapping_set = parse_chebi("chebi_3star.sdf")
    write_sssom(mapping_set, "chebi.sssom.tsv")

    # Access sssom_schema.Mapping objects directly
    for mapping in mapping_set.mappings:
        print(f"{mapping.object_id} → {mapping.subject_id}")

.. toctree::
    :maxdepth: 2
    :caption: Getting Started

    installation
    cli

.. toctree::
    :maxdepth: 2
    :caption: Features

    download
    diff
    exports

.. toctree::
    :maxdepth: 2
    :caption: API Reference

    api
    parsers
    models
