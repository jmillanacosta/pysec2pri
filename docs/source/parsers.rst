#######
Parsers
#######

Database-specific parsers that produce SSSOM-compliant MappingSets.

All parsers inherit from :class:`~pysec2pri.parsers.base.BaseParser`
and return ``sssom_schema.MappingSet`` objects. Anybody can add a new database
by adding the necessary config yaml, parser class, and a DataSourceConfig entry
in constants.py. If possible, you can also add automated download in ``src/pysec2pri/download.py``.

Supported Databases
===================

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Database
     - Input Format
     - Description
   * - **ChEBI**
     - SDF files
     - Parses ChEBI SDF files for secondary→primary ID mappings.
   * - **HMDB**
     - XML files
     - Parses HMDB XML files for metabolite ID mappings.
   * - **HGNC**
     - TSV files
     - Parses HGNC TSV files for gene ID and symbol mappings.
   * - **NCBI Gene**
     - gene_history files
     - Parses NCBI gene_history files for gene ID mappings.
   * - **UniProt**
     - sec_ac.txt files
     - Parses UniProt sec_ac.txt files for protein accession mappings.
   * - **Wikidata**
     - SPARQL queries
     - Queries Wikidata SPARQL for redirect mappings.

Usage
=====

.. code-block:: python

    from pysec2pri.parsers import ChEBIParser

    parser = ChEBIParser()
    mapping_set = parser.parse("ChEBI_complete_3star.sdf")
    # mapping_set is an sssom_schema.MappingSet

Adding a New Parser
===================

To add support for a new database:

1. **Create config YAML** (``config/mydb.yaml``):

   .. code-block:: yaml

      # Metadata for the mapping set
      mappingset:
        curie_map:
          MYDB: https://example.org/mydb/
          # ... standard prefixes
        mapping_set_id: https://example.org/mydb
        mapping_set_title: MyDB Secondary to Primary Mapping
        license: https://example.org/mydb/license
        # ... other SSSOM metadata

      # Metadata applied to each mapping
      mapping:
        mapping_justification: semapv:BackgroundKnowledgeBasedMatching

      # Download URLs
      download_urls:
        main: https://example.org/mydb/data.txt

2. **Create parser class** (``src/pysec2pri/parsers/mydb.py``):

   .. code-block:: python

      from pysec2pri.parsers.base import BaseParser, Sec2PriMappingSet

      class MyDBParser(BaseParser):
          datasource_name = "mydb"

          def parse(self, input_path):
              # Parse your file format
              mappings = self._extract_mappings(input_path)
              return self.create_mapping_set(mappings, mapping_type="id")

3. **Register in constants** (``src/pysec2pri/constants.py``):

   .. code-block:: python

      MYDB = get_datasource_config("mydb")
      ALL_DATASOURCES = [..., MYDB]

4. **Add CLI command and API function** (``src/pysec2pri/cli.py``):

   .. code-block:: python

      @app.command()
      def mydb(...):
          ...

Module Reference
================

.. automodule:: pysec2pri.parsers
    :members:
