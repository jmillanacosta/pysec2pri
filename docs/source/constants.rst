###########
 Constants
###########

Configuration and constants for datasource handling.

Datasource Configuration
========================

.. autoclass:: pysec2pri.constants.DatasourceConfig
    :members:

Pre-defined datasource configurations:

.. autodata:: pysec2pri.constants.CHEBI
.. autodata:: pysec2pri.constants.HMDB
.. autodata:: pysec2pri.constants.HGNC
.. autodata:: pysec2pri.constants.NCBI
.. autodata:: pysec2pri.constants.UNIPROT

Prefix Maps
===========

.. autodata:: pysec2pri.constants.STANDARD_PREFIX_MAP

.. autodata:: pysec2pri.constants.MAPPING_JUSTIFICATION

Helper Functions
================

.. autofunction:: pysec2pri.constants.get_datasource_config
.. autofunction:: pysec2pri.constants.get_prefix_map_for_datasource
