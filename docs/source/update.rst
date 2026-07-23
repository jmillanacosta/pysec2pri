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

   * - ``column_name``
     - ``column_name_primary``
   * - secondary
     - its primary
   * - primary
     - itself
   * - not from this source
     - itself
   * - both retired and current (ambiguous)
     - nothing: an empty cell

(``update_labels`` is identical, but writes ``column_name_current`` by default.)

Ambiguity
=========

Ambiguity: ``ABC`` can be a retired symbol that
became ``DEF``, *and* the current symbol of another gene.

.. code-block:: python

    df[df["gene_sym_primary"] == ""]

Would return all ambiguous cases that could not be solved when updating `df["gene_sym"]`, so ``ABC`` becomes `None`.

To resolve these cases, you can give a hint: another column value of the same row that says which
entry it means.

Hints from a crosswalk
======================

``synonyms`` names a column of names for the row. ``xref`` names a column of
identifiers from another vocabulary, such as an Ensembl gene ID next to
the HGNC to be solved.

An xref hint needs a table stating which entry each of those identifiers belongs
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

Now an ambiguous ``gene_id`` is decided by that row's ``ensembl``. Hints are
only consulted for ambiguous rows, so a wrong one can't mess a row value that was
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
