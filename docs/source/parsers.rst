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
   * - **VGNC**
     - ``all_vgnc_gene_set_All.tsv``, ``all_vgnc_withdrawn.tsv``
     - ``parse()``, ``parse_labels()``, ``parse_all()``
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

.. autoclass:: pysec2pri.parsers.vgnc.VGNCParser
    :members:

.. autoclass:: pysec2pri.parsers.wikidata.WikidataParser
    :members:

Adding a New Parser
===================

See :doc:`adding_a_source`.