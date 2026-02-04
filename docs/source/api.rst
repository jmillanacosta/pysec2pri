###
API
###

Functions for parsing biological databases and generating
SSSOM-compliant mappings.

All functions return SSSOM-schema-compliant ``MappingSet`` objects
for integration with the SSSOM ecosystem.

Parsing Functions
=================

.. autofunction:: pysec2pri.api.parse_chebi
.. autofunction:: pysec2pri.api.parse_chebi_synonyms
.. autofunction:: pysec2pri.api.parse_hmdb
.. autofunction:: pysec2pri.api.parse_hgnc
.. autofunction:: pysec2pri.api.parse_hgnc_symbols
.. autofunction:: pysec2pri.api.parse_ncbi
.. autofunction:: pysec2pri.api.parse_ncbi_symbols
.. autofunction:: pysec2pri.api.parse_uniprot

Export Functions
================

.. autofunction:: pysec2pri.exports.write_sssom
.. autofunction:: pysec2pri.exports.write_sec2pri
.. autofunction:: pysec2pri.exports.write_name2synonym
.. autofunction:: pysec2pri.exports.write_symbol2prev
