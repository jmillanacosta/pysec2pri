#######
Parsers
#######

Database-specific parsers that produce SSSOM-compliant MappingSets.

All parsers inherit from :class:`~pysec2pri.parsers.base.BaseParser`
and return :class:`~pysec2pri.parsers.base.Sec2PriMappingSet` objects.

Anybody can add a new database by adding the necessary config yaml,
parser class, and a DataSourceConfig entry in constants.py.
If possible, you can also add automated download in ``src/pysec2pri/download.py``.

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Database
     - Input
     - Methods
   * - **ChEBI**
     - TSV directory (≥ release 245) or SDF file (< 245)
     - ``parse()``, ``parse_synonyms()``
   * - **Ensembl**
     - ``stable_id_event``, ``mapping_session``, ``gene``, ``xref``,
       ``external_synonym``
     - ``parse()``, ``parse_labels()``, ``parse_all()``
   * - **HMDB**
     - ``hmdb_metabolites.xml`` or ``hmdb_proteins.xml``
     - ``parse()``
   * - **HGNC**
     - ``hgnc_complete_set.txt``, ``withdrawn.txt``
     - ``parse()``, ``parse_labels()``, ``parse_all()``
   * - **NCBI Gene**
     - ``gene_history``, ``gene_info``
     - ``parse()``, ``parse_labels()``, ``parse_all()``
   * - **UniProt**
     - ``sec_ac.txt``, ``delac_sp.txt``
     - ``parse()``
   * - **Wikidata**
     - SPARQL endpoint (live) or pre-fetched JSON
     - ``parse()``, ``parse_all()``, ``parse_from_file()``

Module Reference
================

.. autoclass:: pysec2pri.parsers.chebi.ChEBIParser
    :members:

.. autoclass:: pysec2pri.parsers.ensembl.EnsemblParser
    :members:

.. autoclass:: pysec2pri.parsers.hmdb.HMDBParser
    :members:

.. autoclass:: pysec2pri.parsers.hgnc.HGNCParser
    :members:

.. autoclass:: pysec2pri.parsers.ncbi.NCBIParser
    :members:

.. autoclass:: pysec2pri.parsers.uniprot.UniProtParser
    :members:

.. autoclass:: pysec2pri.parsers.wikidata.WikidataParser
    :members:

Adding a New Parser
===================

1. **Create config YAML** (``config/mydb.yaml``) with ``mappingset``,
   ``mapping``, and ``download_urls`` sections - see existing configs for
   reference.

2. **Create parser class** (``src/pysec2pri/parsers/mydb.py``):

   .. code-block:: python

      from pysec2pri.parsers.base import BaseParser

      class MyDBParser(BaseParser):
          datasource_name = "mydb"

          def parse(self, input_path):
              raw = self._load(input_path)
              mappings = self._build_id_mappings(raw)
              return self._create_mapping_set(mappings, mapping_type="id")

3. **Register in constants** (``src/pysec2pri/constants.py``):

   .. code-block:: python

      MYDB = get_datasource_config("mydb")
      ALL_DATASOURCES = [..., MYDB]

4. **Expose in API and CLI** - add a ``parse_mydb()`` function in
   ``src/pysec2pri/api.py`` and a command in ``src/pysec2pri/cli.py``.
