########
 Models
########

Data models for representing secondary-to-primary identifier mappings.

Enums
=====

.. autoclass:: pysec2pri.models.MappingCardinality
    :members:
    :undoc-members:

.. autoclass:: pysec2pri.models.PredicateID
    :members:
    :undoc-members:

Base Classes
============

.. autoclass:: pysec2pri.models.BaseMapping
    :members:
    :show-inheritance:

.. autoclass:: pysec2pri.models.IdMapping
    :members:
    :show-inheritance:

.. autoclass:: pysec2pri.models.SymbolMapping
    :members:
    :show-inheritance:

Datasource Classes
==================

.. autoclass:: pysec2pri.models.ChEBIMapping
    :members:
    :show-inheritance:

.. autoclass:: pysec2pri.models.HMDBMapping
    :members:
    :show-inheritance:

.. autoclass:: pysec2pri.models.HGNCMapping
    :members:
    :show-inheritance:

.. autoclass:: pysec2pri.models.NCBIGeneMapping
    :members:
    :show-inheritance:

.. autoclass:: pysec2pri.models.UniProtMapping
    :members:
    :show-inheritance:

Container
=========

.. autoclass:: pysec2pri.models.MappingSet
    :members:
    :show-inheritance:

Utilities
=========

.. autofunction:: pysec2pri.models.compute_cardinality
