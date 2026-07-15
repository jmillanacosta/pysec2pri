#####################
Update IDs and labels
#####################

Take a table of your own and add a column with each value's current ID or
label.

.. code-block:: python

    import pandas as pd
    from pysec2pri import generate_ids, update_ids

    df = pd.read_csv("genes.tsv", sep="\t")
    out = update_ids(df, generate_ids("hgnc"), at="gene_id")

``gene_id_primary`` is the new column. :func:`~pysec2pri.update.update_labels`
is the same for symbols, and writes ``_current``.

.. code-block:: bash

    pysec2pri update-ids genes.tsv hgnc --at gene_id -o out.tsv
    pysec2pri update-labels genes.tsv hgnc --at symbol -o out.tsv


Every row keeps its original value, and the new column shows the solved value:

.. list-table::
   :header-rows: 1

   * - the value is
     - the new column holds
   * - retired
     - what it retired into
   * - current
     - itself
   * - not from this source
     - itself
   * - both retired and current
     - nothing: an empty cell

Ambiguity
=========

Ambiguity: ``HGNC:2`` can be a retired ID that
became ``HGNC:3``, *and* the current ID of a different gene. One row says
``HGNC:2`` and both readings are correct.

pysec2pri does not guess. The cell is left empty, and an empty cell is how you
find the rows worth looking at:

.. code-block:: python

    out[out["gene_id_primary"] == ""]

To resolve them, you can give a hint: another column of the same row that says which
entry it means.

Hints from a crosswalk
======================

``synonyms`` names a column of names for the row. ``xref`` names a column of
identifiers from another vocabulary, such as an Ensembl gene ID next to
the HGNC to be solved.

An xref hint needs a table saying which entry each of those identifiers belongs
to. Bring your own with :func:`~mapkgsutils.context.load_xref_mapping`. Any
SSSOM file works, and so does a plain TSV of ``subject_id`` and ``object_id``:

.. code-block:: text

    subject_id       object_id
    ENSG00000121410  HGNC:5
    ENSG00000175899  HGNC:7

.. code-block:: python

    from mapkgsutils.context import load_xref_mapping

    out = update_ids(
        df,
        generate_ids("hgnc"),
        at="gene_id",
        xref="ensembl",
        xref_mapping=load_xref_mapping("ensembl_to_hgnc.tsv"),
    )

Now an ambiguous ``gene_id`` is decided by that row's ``ensembl`` value: the
crosswalk says which gene it is, and that gene's ID is the answer. Hints are
only consulted for ambiguous rows, so a wrong one cannot spoil a row that was
already clear.

Using a hint also adds a ``gene_id_primary_id`` column, holding the ID the row
resolved to.

The same from the command line, with ``--xref-file`` for your own table:

.. code-block:: bash

    pysec2pri update-ids genes.tsv hgnc --at gene_id \
        --xref ensembl --xref-file ensembl_to_hgnc.tsv

``--xref-source`` downloads a crosswalk the source's config already lists,
instead of you supplying one. ``--xref-on`` says which vocabulary your column
holds:

.. code-block:: bash

    pysec2pri update-ids genes.tsv hgnc --at gene_id \
        --xref ensembl --xref-source hgnc_custom --xref-on ensembl

Report decisions
=========================

``report_path`` (``--report``) writes down every hint it considered, whether it
was accepted, and why:

.. code-block:: text

    stage         token   predicate_id  candidate  accepted  reason
    xref_filter   ENSG_B                HGNC:3     True      no predicate given, assumed equivalence

It is a TSV, so read it back with pandas:

.. code-block:: python

    out = update_ids(
        df,
        generate_ids("hgnc"),
        at="gene_id",
        xref="ensembl",
        xref_mapping=load_xref_mapping("ensembl_to_hgnc.tsv"),
        report_path="decisions.tsv",
    )

    report = pd.read_csv("decisions.tsv", sep="\t")
    report[~report["accepted"]]  # the hints that did not settle a row

One row per hint considered, not one per row of your data: a row with no
ambiguity never consults a hint and never appears here.

Any equivalence predicate is accepted by default, including records with none.
``xref_predicates`` (``--xref-predicate``) narrows that to the ones you name.

Reference
=========

.. automodule:: pysec2pri.update
    :undoc-members:
    :show-inheritance:
