Diff Module
===========

The diff module provides utilities for comparing mapping sets to identify
additions, removals, and changes between releases. It uses Polars for
efficient operations on large mapping sets.

Comparing Mapping Sets
----------------------

Compare two MappingSet objects:

.. code-block:: python

    from pysec2pri import parse_chebi, diff_mapping_sets, summarize_diff

    # Parse two versions
    old_mappings = parse_chebi("chebi_v220.sdf")
    new_mappings = parse_chebi("chebi_v221.sdf")

    # Compare
    diff = diff_mapping_sets(old_mappings, new_mappings)

    # Get summary
    summary = summarize_diff(diff)
    print(summary)

    # Access detailed dataframes
    print(f"Added mappings: {diff.added_count}")
    print(f"Removed mappings: {diff.removed_count}")

    # View the actual changes (Polars DataFrames)
    print(diff.added.head())
    print(diff.removed.head())

Comparing SSSOM Files
---------------------

Compare existing SSSOM files directly:

.. code-block:: python

    from pysec2pri import diff_sssom_files

    diff = diff_sssom_files(
        "mappings_v1.sssom.tsv",
        "mappings_v2.sssom.tsv"
    )

    # Save diff results
    diff.added.write_csv("added_mappings.tsv", separator="\t")
    diff.removed.write_csv("removed_mappings.tsv", separator="\t")

CLI Usage
---------

The diff functionality is also available via the command line:

.. code-block:: bash

    # Compare two SSSOM files
    pysec2pri diff compare old.sssom.tsv new.sssom.tsv

    # Generate detailed report
    pysec2pri diff compare old.sssom.tsv new.sssom.tsv --output-dir ./diff_results/

MappingDiff Structure
---------------------

The ``MappingDiff`` dataclass contains:

- ``added``: Polars DataFrame of mappings present in new but not in old
- ``removed``: Polars DataFrame of mappings present in old but not in new
- ``added_count``: Number of added mappings
- ``removed_count``: Number of removed mappings
- ``old_total``: Total mappings in old set
- ``new_total``: Total mappings in new set

API Reference
-------------

.. automodule:: pysec2pri.diff
   :members:
   :undoc-members:
   :show-inheritance:
