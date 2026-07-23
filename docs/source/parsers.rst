#######
Parsers
#######

Database-specific parsers that produce SSSOM-compliant MappingSets.

All parsers inherit from :class:`~pysec2pri.parsers.base.BaseParser`
and return :class:`~pysec2pri.parsers.base.BaseMappingSet` objects.

Anybody can add a new database by adding the necessary config yaml,
parser class, and a DataSourceConfig entry in constants.py.
If possible, you can also add automated download in ``src/pysec2pri/download.py``.

.. list-table::
   :header-rows: 1
   :widths: 15 30 15 15 20 20

   * - Database
     - Input
     - ``ids``
     - ``labels``
     - ``primary_ids``
     - ``primary_labels``
   * - **ChEBI**
     - TSV directory (≥ release 245) or SDF file (< 245)
     - ``parse()``
     - ``parse_synonyms()``
     - ``parse_primary_ids()``
     - ``parse_primary_labels()``
   * - **Ensembl**
     - ``stable_id_event``, ``mapping_session``, ``gene``, ``xref``,
       ``external_synonym``
     - ``parse()``
     - ``parse_labels()``
     - ``parse_primary_ids()``
     - ``parse_primary_labels()``
   * - **HMDB**
     - ``hmdb_metabolites.xml`` or ``hmdb_proteins.xml``
     - ``parse()``
     - --
     - ``parse_primary_ids()``
     - --
   * - **HGNC**
     - ``hgnc_complete_set.txt``, ``withdrawn.txt``
     - ``parse()``
     - ``parse_labels()``
     - ``parse_primary_ids()``
     - ``parse_primary_labels()``
   * - **NCBI Gene**
     - ``gene_history``, ``gene_info``
     - ``parse()``
     - ``parse_labels()``
     - ``parse_primary_ids()``
     - ``parse_primary_labels()``
   * - **UniProt**
     - ``sec_ac.txt``, ``delac_sp.txt``
     - ``parse()``
     - --
     - ``parse_primary_ids()``
     - --
   * - **VGNC**
     - ``all_vgnc_gene_set_All.tsv``, ``all_vgnc_withdrawn.tsv``
     - ``parse()``
     - ``parse_labels()``
     - ``parse_primary_ids()``
     - ``parse_primary_labels()``
   * - **Wikidata**
     - SPARQL endpoint (live) or pre-fetched JSON
     - ``parse_all()``
     - ``parse_labels()``
     - --
     - --

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

.. autoclass:: pysec2pri.parsers.vgnc.VGNCParser
    :members:

.. autoclass:: pysec2pri.parsers.wikidata.WikidataParser
    :members:

Adding a New Parser
===================

See :doc:`adding_a_source`.
