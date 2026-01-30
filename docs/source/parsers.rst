Parsers
=======

This section documents the parsers used in pysec2pri for each supported database.

ChEBIParser
-----------
Parses ChEBI SDF files to extract secondary-to-primary identifier mappings and name-to-synonym relationships. Accepts SDF files as input.

HMDBParser
----------
Parses HMDB XML files (or ZIP containing XML) to extract secondary-to-primary identifier mappings and synonym relationships. Uses streaming XML parsing for large files.

HGNCParser
----------
Parses HGNC withdrawn and complete set TSV files to extract secondary-to-primary gene identifier mappings and symbol relationships. Handles legacy and current file formats.

NCBIParser
----------
Parses NCBI gene_history and gene_info TSV files to extract secondary-to-primary gene identifier mappings and gene symbols. Supports compressed files and taxon filtering.

UniProtParser
-------------
Parses UniProt sec_ac.txt (secondary accessions) and delac_sp.txt (deleted accessions) files to extract secondary-to-primary protein accession mappings. Can also process FASTA files for primary IDs.

WikidataParser
--------------
Queries Wikidata using SPARQL (via QLever endpoint) to extract redirect mappings for supported entity types. Returns results as Polars DataFrames and MappingSets.

Usage
-----

Each parser is available in `pysec2pri.parsers` and is used internally by the main API functions. You can also use them directly for advanced workflows.

Example:

.. code-block:: python

    from pysec2pri.parsers import ChEBIParser
    parser = ChEBIParser()
    mapping_set = parser.parse("ChEBI_complete_3star.sdf")

See the API documentation for details on each parser's methods and options.
