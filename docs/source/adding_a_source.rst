################
Adding a source
################

Source is **Example**. Two files are needed: a config YAML and a parser.
:mod:`pysec2pri.api` and :mod:`pysec2pri.cli` read the config, so
they pick the newly added source up on their own.

``docs/source/example/`` shows a working example to copy.

Every source can produce:

.. list-table::
   :header-rows: 1

   * - kind
     - maps
     - example
   * - ``ids``
     - old ID -> current ID
     - ``HGNC:31139`` -> ``HGNC:11979``
   * - ``labels``
     - old symbol -> current symbol
     - ``TMEM8`` -> ``TMEM8A``


****************************
1. ``config/example.yaml``
****************************

Copy ``docs/source/example/example.yaml`` to
``src/pysec2pri/config/<source>.yaml`` and edit it.

.. literalinclude:: example/example.yaml
   :language: yaml

``inputs`` (see under ``mapping_sets``) specifies which files to download and which
argument each one becomes. Each key also becomes a command line option, so
``withdrawn`` gives you ``--withdrawn`` in the CLI.

Datasets within a release
=========================

Most sources publish one dataset per release. HGNC is one, so the source and
the version name it completely::

    https://sec2pri.github.io/hgnc/2024-01-01

Some publish several. ChEBI release 245 is both a ``3star`` subset and a
``complete`` one. Ensembl release 115 is a different set of mappings for every
species. Those are separate mapping sets, so they need separate IRIs, and the
version alone no longer tells them apart.

The **slug** is the extra bit of the IRI that does::

    https://sec2pri.github.io/chebi/245/3star
                                        ^^^^^ the slug

``products`` names the parser arguments the IRI is built from as slugs:

.. code-block:: yaml

    products: ["subset"]

Now each subset of release 245 gets its own IRIs:

.. list-table::
   :header-rows: 1

   * - ``--subset``
     - ``mapping_set_id``
     - ``record_id``
   * - ``3star``
     - ``https://sec2pri.github.io/chebi/245/3star``
     - ``sec2pri:chebi/245/3star/<hash>``
   * - ``complete``
     - ``https://sec2pri.github.io/chebi/245/complete``
     - ``sec2pri:chebi/245/complete/<hash>``

Name several and each adds a segment, in the order you list them:
``products: ["species", "subset"]`` gives ``.../115/9606/3star`` for ChEBI 115.

Declare an argument when two runs produce mappings that
should not be pooled, e.g. when generating species-specific mapping sets for Ensembl.

**************************
1. ``parsers/example.py``
**************************

Each source publishes their deprecation and synonym release in a different way. In this example, there are two files, and the parser uses them for different things:

* **the withdrawn file**: the retirements themselves. These become the
  mappings.
* **the complete set**: every current entry. Not mappings; the list of what
  is live right now.

The complete set of IDs and labels allows:

* **Completeness**: the withdrawn file only names a primary if something
  retired into it. The complete set names all of them, so ``to_pri_ids()``
  returns everything.
* **Ambiguity**: an ID that is retired *and* still current is ambiguous. You
  can only see that against the full current list.

Copy ``docs/source/example/example.py`` to
``src/pysec2pri/parsers/<source>.py``. Name the class to match
``parser_class``, and return ``self.create_mapping_set(...)``:

.. literalinclude:: example/example.py
   :language: python

Every method argument named in ``inputs`` gets its file.

Each name in ``products`` becomes a command line option too, so
``products: ["subset"]`` gives ``--subset``.

**********************************************************
1. ``downloads/example.py``: only if URLs are not fixed
**********************************************************

If Example publishes the same URLs every time, skip this. The ``download_urls``
in the config are the download logic: each one is fetched, unzipped if needed,
and handed to the parser.

Write a downloader only when:

* the URL needs a replaced string (version or a species, etc.) in it,
* you need to find the latest release, or list every release,
* you need the real release date for SSSOM ``mapping_date``.

.. literalinclude:: example/example_downloads.py
   :language: python

Then add it to the registries in ``downloads/__init__.py``.

***************************
1. ``parsers/__init__.py``
***************************

Export the parser to three files:

.. code-block:: python

    _LAZY_EXPORTS = {"ExampleParser": "pysec2pri.parsers.example"}
    # under TYPE_CHECKING:  from pysec2pri.parsers.example import ExampleParser
    __all__ = ["ExampleParser"]

**********
5. Tests
**********

Put a small sample of the real files in ``tests/data/`` and write a test to parse them. The
config tests already check ``parser_class``, every ``method``, and
every ``inputs`` argument are real, and that each ``inputs`` key is
downloadable.

************
Checklist
************

* ``config/example.yaml``: ``parser_class``, ``download_urls``, ``mapping_sets``
* ``parsers/example.py``: ``ExampleParser`` with the methods the config names
* ``downloads/example.py``: only if the URLs are not fixed
* ``parsers/__init__.py``: export the parser in all three places
* ``tests/data/`` / ``tests/test_parsers.py``: a sample file and a test
