Wikidata Integration
====================

The pysec2pri package includes support for querying Wikidata to retrieve
redirect mappings for metabolites, genes, and proteins. This uses SPARQL
queries against the QLever endpoint.

Parsing Wikidata
----------------

Query Wikidata for redirect mappings:

.. code-block:: python

    from pysec2pri import parse_wikidata

    # Query metabolite redirects (ChEBI, HMDB, PubChem, etc.)
    metabolites = parse_wikidata("metabolites")

    # Query gene redirects
    genes = parse_wikidata("genes")

    # Query protein redirects
    proteins = parse_wikidata("proteins")

Using the WikidataParser
------------------------

For more control, use the WikidataParser class directly:

.. code-block:: python

    from pysec2pri import WikidataParser

    parser = WikidataParser()

    # Execute custom SPARQL query
    results = parser.query_sparql('''
        SELECT ?item ?redirect WHERE {
            ?redirect owl:sameAs ?item .
            ?item wdt:P31 wd:Q11173 .  # chemical compounds
        }
        LIMIT 100
    ''')

    # Parse with predefined query type
    mapping_set = parser.parse(query_type="metabolites")

CLI Usage
---------

Query Wikidata from the command line:

.. code-block:: bash

    # Query metabolite redirects
    pysec2pri wikidata metabolites --output metabolite_redirects.sssom.tsv

    # Query gene redirects
    pysec2pri wikidata genes --output gene_redirects.sssom.tsv

    # Query protein redirects
    pysec2pri wikidata proteins --output protein_redirects.sssom.tsv

Available Query Types
---------------------

- ``metabolites``: Wikidata items with ChEBI, HMDB, or PubChem IDs
- ``genes``: Wikidata items with NCBI Gene or HGNC IDs
- ``proteins``: Wikidata items with UniProt IDs

SPARQL Endpoint
---------------

By default, pysec2pri uses the QLever endpoint for Wikidata:
``https://qlever.dev/api/wikidata``

This provides fast query performance for large result sets.

API Reference
-------------

.. automodule:: pysec2pri.parsers.wikidata
   :members:
   :undoc-members:
   :show-inheritance:
