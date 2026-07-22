################
Adding a source
################

To add an **Example** source you need a config YAML and a parser module.
:mod:`pysec2pri.api` and :mod:`pysec2pri.cli` will pick up the new source
automatically once those are in place.

The goal of the parser is to extract all primary IDs, all primary labels, and
all available secondary mappings for them and generate an ``IdMappingSet``\s
and/or a ``LabelMappingSet``\s:

``docs/source/example/`` contains a working example to copy from.



.. list-table::
   :header-rows: 1

   * - type
     - maps
     - example
   * - ``ids``
     - old ID → current ID
     - ``HGNC:31139`` → ``HGNC:11979``
   * - ``labels``
     - old symbol → current symbol
     - ``TMEM8`` → ``TMEM8A``


************
0. Checklist
************

* ``src/pysec2pri/config/example.yaml``: ``parser_class``, ``download_urls``,
  ``mapping_sets``
* ``src/pysec2pri/parsers/example.py``: ``ExampleParser`` with every method
  the config names
* ``src/pysec2pri/parsers/__init__.py``: parser class in ``_LAZY_EXPORTS``,
  the ``TYPE_CHECKING`` block, and ``__all__``
* ``src/pysec2pri/downloads/example.py``: only if URLs are not fixed
* ``src/pysec2pri/downloads/__init__.py``: downloader in ``CHECK_RELEASE``,
  ``DOWNLOADERS``, ``URLS_AND_DATE``, and ``__all__`` (only if a downloader was
  written)
* ``tests/data/`` + ``tests/test_example.py``: sample file and a parsing test


****************************
1. ``config/example.yaml``
****************************

Copy ``docs/source/example/example.yaml`` to
``src/pysec2pri/config/<source>.yaml`` and edit it.

.. literalinclude:: example/example.yaml
   :language: yaml

``inputs`` (under ``mapping_sets``) specifies which files to download and
which argument each one becomes in :func:`~pysec2pri.generate_ids` and
:func:`~pysec2pri.generate_labels`.

Datasets within a release
=========================

Most sources publish one dataset per release. HGNC is one, so the source and
the version name it completely::

    https://sec2pri.bigcat-bioinformatics.org/hgnc/2024-01-01

Some publish several. ChEBI release 245 is both a ``3star`` subset and a
``complete`` one. Ensembl release 115 is a different set of mappings for every
species. Those are separate mapping sets with separate IRIs.

``products`` names the parser arguments that also become slugs in the generated
mapping document IRI and the set of ``record_id``\s it contains.

.. code-block:: yaml

    products: ["subset"]

Means that each subset of release 245 gets its own IRI pattern:

.. list-table::
   :header-rows: 1

   * - ``--subset``
     - ``mapping_set_id``
     - ``record_id``
   * - ``3star``
     - ``https://sec2pri.bigcat-bioinformatics.org/chebi/245/3star``
     - ``sec2pri:chebi/245/3star/<hash>``
   * - ``complete``
     - ``https://sec2pri.bigcat-bioinformatics.org/chebi/245/complete``
     - ``sec2pri:chebi/245/complete/<hash>``

Each name in ``products`` also becomes a CLI option (e.g. ``--subset``).
``record_id`` is built automatically from each row's ``subject_id``/
``object_id``; a row can include its own ``subset`` (or ``species``, etc.).


**************************
2. ``parsers/example.py``
**************************

Each source publishes its deprecation and synonym release files differently.
In this simplified example there are two files:

* **the withdrawn file**: the secondary mappings.
* **the complete set**: every current entry ID and primary label.

But in your case the whole information could be extracted from one or several files.

Copy ``docs/source/example/example.py`` to
``src/pysec2pri/parsers/<source>.py``. Name the class to match the
``parser_class`` value in the config YAML.

:func:`~pysec2pri.Example.ExampleParser.parse` must return the
``IdMappingSet`` and
:func:`~pysec2pri.Example.ExampleParser.parse_labels` the
``LabelMappingSet``.

.. literalinclude:: example/example.py
   :language: python


*****************************
3. ``parsers/__init__.py``
*****************************

The parser package uses **lazy imports** (PEP 562): the heavy dependencies
(``polars``, ``sssom``, ``rdflib``) are only loaded the first time a parser
class is actually used. For this to work, every parser class must be registered
in three places in ``src/pysec2pri/parsers/__init__.py``.

``_LAZY_EXPORTS``
=================

This dict maps each public name to the submodule that defines it. The
module-level ``__getattr__`` hook reads it to import and cache the class on
first access, so ``pysec2pri.parsers.ExampleParser`` resolves correctly without
eagerly importing the module at startup.

.. code-block:: python

    _LAZY_EXPORTS: dict[str, str] = {
        ...
        "ExampleParser": "pysec2pri.parsers.example",   # ← add this
    }

``TYPE_CHECKING`` block
=======================

Static analysis tools (mypy, Pyright) and IDEs never execute ``__getattr__``,
so they cannot resolve lazily-exported names. The ``if TYPE_CHECKING:`` block
provides the explicit import they need. It is never executed at runtime, so it
costs nothing.

.. code-block:: python

    if TYPE_CHECKING:
        ...
        from pysec2pri.parsers.example import ExampleParser   # ← add this

``__all__``
===========

``__all__`` declares the public API of the package. It is used by
``from pysec2pri.parsers import *``, ``dir(pysec2pri.parsers)``, and
documentation generators. A name absent from ``__all__`` is considered
private even if it can be imported.

.. code-block:: python

    __all__ = [
        ...
        "ExampleParser",   # ← add this
    ]


*******************************************************
4. ``downloads/example.py``: only if URLs are not fixed
*******************************************************

If Example publishes the same static URLs every release, skip this step.
The ``download_urls`` in the config are sufficient: each URL is fetched,
decompressed if needed, and handed to the parser.

Write a downloader only when:

* the URL contains a version string or species slug that must be substituted,
* you need to find or list available releases programmatically, or
* you need the real upstream release date for the SSSOM ``mapping_date`` field.

.. literalinclude:: example/example_downloads.py
   :language: python

``downloads/__init__.py`` registries
=====================================

After writing the downloader, register it in the three dicts in
``src/pysec2pri/downloads/__init__.py``.

**``CHECK_RELEASE``** needed by
:func:`mapkgsutils.download.get_latest_release_info`. Maps the datasource name
to the function that queries the upstream source and returns the latest
available release identifier.

.. code-block:: python

    CHECK_RELEASE = {
        ...
        "example": example.check_example_release,   # ← add this
    }

**``DOWNLOADERS``** needed by :func:`mapkgsutils.download.list_versions`.
Maps the datasource name to the ``BaseDownloader`` subclass.

.. code-block:: python

    DOWNLOADERS: dict[str, type[BaseDownloader]] = {
        ...
        "example": example.ExampleDownloader,   # ← add this (if versioned)
    }

**``URLS_AND_DATE``** needed by
:func:`mapkgsutils.download.resolve_datasource_urls`. Maps the datasource name
to the function that returns the concrete download URLs and the release date for
a given version.

.. code-block:: python

    URLS_AND_DATE = {
        ...
        "example": example.urls_and_date,   # add this
    }

Also add the import of the submodule and the downloader class at the top of
``downloads/__init__.py``, and include the class in ``__all__``.


**********
5. Tests
**********

Put a small sample of the real upstream files in ``tests/data/`` and write a
test that correctly parses them into mapping sets. The config tests in
``tests/test_config.py`` already verify automatically that ``parser_class``,
every ``method``, and every ``inputs`` argument named in the config are real attributes on the parser class, and that each ``inputs`` key corresponds to a downloadable file.
