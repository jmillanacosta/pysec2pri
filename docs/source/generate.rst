######################
 Generate mapping sets
######################

A mapping set says what a database retired. There are two kinds, and one
function each:

.. list-table::
   :header-rows: 1

   * - kind
     - maps
     - example
   * - ``ids``
     - old ID -> current ID
     - ``HGNC:31139`` -> ``HGNC:11979``
   * - ``labels``
     - old or alias label -> current label
     - ``TMEM8`` -> ``TMEM8A``

.. code-block:: python

    from pysec2pri import generate_ids, generate_labels, sources

    sources()          # every source
    sources("labels")  # sources with labels

    hgnc = generate_ids("hgnc")
    chebi = generate_labels("chebi", subset="3star")

Both return an SSSOM mapping set: an ``IdMappingSet`` or a ``LabelMappingSet``.
Subjects are secondary, objects are primary. See :doc:`exports` to write one
out, and :doc:`update_ids` to resolve your own data against it.

The same from the command line, one command per kind:

.. code-block:: bash

    pysec2pri hgnc ids
    pysec2pri hgnc labels

Source options
==============

Some sources publish several datasets per release, and an argument picks which:
``species`` for the multi-species ones, ``subset`` for ChEBI, ``entity_type``
for Wikidata. Each is a mapping set of its own, with its own ``mapping_set_id``.
Passing one to a source that has no such option does nothing, so the same call
works everywhere.

.. code-block:: python

    generate_ids("ensembl", version="115", species="9606")
    generate_ids("hgnc", species="9606")  # HGNC has no species; ignored

``pysec2pri <source> ids --help`` lists what a source takes.

Input files
===========

Files are downloaded unless you pass your own. Each input a source declares
gets its own option, named after the file:

.. code-block:: bash

    pysec2pri hgnc ids --withdrawn withdrawn.txt --complete hgnc_complete_set.txt

.. code-block:: python

    generate_ids("hgnc", inputs={"withdrawn": "withdrawn.txt", "complete": "hgnc_complete_set.txt"})

Looking further back
====================

A mapping set says what the source's current release says. Sources keep
different amounts of their own history, and some drop retired entries
altogether. ``--consolidate`` reads the source's past releases too: it finds
mappings the current release no longer mentions, and gives every mapping the
release it first appeared in.

.. code-block:: bash

    pysec2pri hgnc ids --consolidate -o hgnc.sssom.tsv

.. code-block:: python

    from pysec2pri import supports_consolidate

    supports_consolidate("hgnc", "ids")
    generate_ids("hgnc", consolidate=True)

This downloads every past release, so it is slow. It saves its progress, so
stopping it and running it again picks up where it left off. Use ``cache_dir``
to say where, ``force`` to start over, and ``from_version``/``to_version`` to
read part of the archive only.

Not every source can do this. :func:`~pysec2pri.api.supports_consolidate` says
which can, and ``--consolidate`` only appears in a source's ``--help`` when it
applies.

Full signatures: :func:`~pysec2pri.api.generate_ids`,
:func:`~pysec2pri.api.generate_labels`, :func:`~pysec2pri.api.sources`,
:func:`~pysec2pri.api.supports_consolidate`.
