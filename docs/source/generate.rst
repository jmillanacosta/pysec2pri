######################
 Generate mapping sets
######################

A mapping set contains identifier or label mappings between secondaries (withdrawn, split, merged, redirected, aliased) and primaries.

.. list-table::
   :header-rows: 1

   * - type
     - maps
     - example
   * - ``ids``
     - old ID -> current ID
     - ``HGNC:31139`` -> ``HGNC:11979``
   * - ``labels``
     - old or alias label -> current label
     - ``TMEM8`` -> ``TMEM8A``

Generating secondary-to-primary identifier and label mapping sets:

.. code-block:: python

    from pysec2pri import generate_ids, generate_labels, sources

    sources()          # every source
    sources("labels")  # sources that produce label mapping sets

    hgnc = generate_ids("hgnc")
    chebi = generate_labels("chebi", subset="3star")

The same from the command line, one command for IDs and another for labels:

.. code-block:: bash

    pysec2pri hgnc ids
    pysec2pri hgnc labels

Both return a SSSOM mapping set: an ``IdMappingSet`` or a ``LabelMappingSet``
where the secondary identifier or label is the subject of a statement about its relationship to the primary that replaces it, if any.

See :doc:`update` for instructions on how to check and update the
identifiers and labels in your own data against a generated mapping set.

Output format
=============

SSSOM is the default output, but other formats can be specified
(e.g. a plain secondary-to-primary TSV, or a bare list of current IDs).
``--format`` picks one; ``--format all`` writes every format the mapping set
supports to a directory:

.. code-block:: bash

    pysec2pri hgnc ids --format sec2pri -o hgnc_sec2pri.tsv
    pysec2pri hgnc ids --format all -o hgnc_2026-04-07/

From Python, the same mapping set can be saved with
:meth:`~pysec2pri.parsers.base.BaseMappingSet.save`, or read into memory
without writing anything via its ``to_<format>()`` methods:

.. code-block:: python

    hgnc.save("sec2pri", "hgnc_sec2pri.tsv")
    df = hgnc.to_sec2pri()  # same table, no file written

See :doc:`exports` for the full list of formats per mapping-set kind, and how
to discover which ones a given source supports.


Source options
==============

Some sources publish several datasets per release.

To subset the generation of mapping sets, arguments like  ```species`` for the multi-species ones, ``subset`` for ChEBI, ``entity_type``
for Wikidata result in a mapping set of its own, with its own ``mapping_set_id``.

.. code-block:: python

    generate_ids("ensembl", version="115", species="9606")
    generate_ids("hgnc", species="9606")  # HGNC has no species; ignored

``pysec2pri <source> ids --help`` lists what flags a source takes.

Input files
===========

Files are downloaded according to the specified ``version``:

.. code-block:: bash

    pysec2pri hgnc ids

.. code-block:: bash

    pysec2pri hgnc ids --withdrawn withdrawn.txt --complete hgnc_complete_set.txt

You can see the available versions for a database:

.. code-block:: bash

    pysec2pri list-versions hgnc

You can also provide the needed input files:

.. code-block:: python

    generate_ids("hgnc", inputs={"withdrawn": "withdrawn.txt", "complete": "hgnc_complete_set.txt"})

Consolidated mapping sets
=========================

A mapping set generated for a specific release of a database is limited to the deprecation information that can be extracted from that release.

``--consolidate`` reads the source's past releases too: it finds
mappings the current release might no longer mention, and gives every mapping the
release version and date it first appeared in.

.. code-block:: bash

    pysec2pri hgnc ids --consolidate -o hgnc.sssom.tsv

.. code-block:: python

    from pysec2pri import supports_consolidate

    supports_consolidate("hgnc", "ids")
    generate_ids("hgnc", consolidate=True)

