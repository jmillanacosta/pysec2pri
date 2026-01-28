Download Module
===============

The download module provides utilities for checking datasource releases
and downloading source files for parsing.

Release Checking
----------------

Check if a new release is available for a datasource:

.. code-block:: python

    from pysec2pri import check_release, get_latest_release_info

    # Check specific datasource
    info = check_release("chebi")
    print(f"Latest release: {info.version}")
    print(f"Last modified: {info.last_modified}")
    print(f"Download URL: {info.url}")

    # Get release info for any datasource
    info = get_latest_release_info("hmdb")

Downloading Files
-----------------

Download source files for parsing:

.. code-block:: python

    from pysec2pri import download_datasource, download_file

    # Download to specific directory
    filepath = download_datasource("chebi", output_dir="./data")

    # Download with custom filename
    filepath = download_file(
        "https://example.com/data.txt",
        output_path="./data/custom_name.txt"
    )

CLI Usage
---------

The download functionality is also available via the command line:

.. code-block:: bash

    # Check for new release
    pysec2pri download check chebi

    # Download source file
    pysec2pri download fetch chebi --output ./data/

    # Download all datasources
    pysec2pri download fetch all --output ./data/

API Reference
-------------

.. automodule:: pysec2pri.download
   :members:
   :undoc-members:
   :show-inheritance:
