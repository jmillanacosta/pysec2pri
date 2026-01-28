<!--
<p align="center">
  <img src="https://github.com/sec2pri/pysec2pri/raw/main/docs/source/logo.png" height="150">
</p>
-->

<h1 align="center">
  pySec2Pri
</h1>

<p align="center">
    <a href="https://github.com/sec2pri/pysec2pri/actions/workflows/tests.yml">
        <img alt="Tests" src="https://github.com/sec2pri/pysec2pri/actions/workflows/tests.yml/badge.svg" /></a>
    <a href="https://pypi.org/project/pysec2pri">
        <img alt="PyPI" src="https://img.shields.io/pypi/v/pysec2pri" /></a>
    <a href="https://pypi.org/project/pysec2pri">
        <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/pysec2pri" /></a>
    <a href="https://github.com/sec2pri/pysec2pri/blob/main/LICENSE">
        <img alt="PyPI - License" src="https://img.shields.io/pypi/l/pysec2pri" /></a>
    <a href='https://pysec2pri.readthedocs.io/en/latest/?badge=latest'>
        <img src='https://readthedocs.org/projects/pysec2pri/badge/?version=latest' alt='Documentation Status' /></a>
</p>

Create mapping files for secondary (retired/withdrawn) biological database
identifiers to primary (current) identifiers.

Outputs mappings in [SSSOM format](https://w3id.org/sssom).

## Supported Databases

| Database  | Description                              | Auto-download |
| --------- | ---------------------------------------- | ------------- |
| ChEBI     | Chemical Entities of Biological Interest | Supported            |
| HMDB      | Human Metabolome Database                | Supported            |
| HGNC      | HUGO Gene Nomenclature Committee         | Supported            |
| NCBI Gene | Entrez Gene                              | Supported            |
| UniProt   | Protein sequence database                | Supported            |
| Wikidata  | Redirect mappings via SPARQL             | Supported            |

## Installation

```console
pip install pysec2pri
```

Or install from source:

```console
pip install git+https://github.com/sec2pri/pysec2pri.git
```

## Quick Start

**Auto-download and process** (recommended):

```console
# Download latest ChEBI and generate mappings
pysec2pri chebi -o chebi_sec2pri.sssom.tsv

# Download latest HMDB and generate mappings
pysec2pri hmdb -o hmdb_sec2pri.sssom.tsv
```

**Use local files**:

```console
# Use a local file instead of downloading
pysec2pri chebi --input ChEBI_complete_3star.sdf -o chebi_sec2pri.sssom.tsv
```

## Command Line Interface

All commands auto-download source files by default. Use input options to
specify local files.

### ChEBI

```console
# Auto-download (recommended)
pysec2pri chebi -o chebi_sec2pri.sssom.tsv

# Use local file
pysec2pri chebi --input ChEBI_complete_3star.sdf -o chebi_sec2pri.sssom.tsv
```

### HMDB

```console
# Auto-download
pysec2pri hmdb -o hmdb_sec2pri.sssom.tsv

# Use local file
pysec2pri hmdb --input hmdb_metabolites.zip -o hmdb_sec2pri.sssom.tsv
```

### HGNC

```console
# Auto-download
pysec2pri hgnc -o hgnc_sec2pri.sssom.tsv

# Use local files
pysec2pri hgnc --withdrawn withdrawn.txt --complete-set hgnc_complete_set.txt
```

### NCBI Gene

```console
# Auto-download (human genes by default)
pysec2pri ncbi -o ncbi_sec2pri.sssom.tsv

# Different organism (mouse)
pysec2pri ncbi --tax-id 10090 -o ncbi_mouse_sec2pri.sssom.tsv

# Use local files
pysec2pri ncbi --history gene_history.gz --gene-info gene_info.gz
```

### UniProt

```console
# Auto-download
pysec2pri uniprot -o uniprot_sec2pri.sssom.tsv

# Use local files
pysec2pri uniprot --sec-ac sec_ac.txt --delac delac_sp.txt
```

### Wikidata (SPARQL query)

```console
# Query metabolite redirects
pysec2pri wikidata --type metabolites -o wikidata_metabolites.sssom.tsv

# Query gene redirects
pysec2pri wikidata --type genes -o wikidata_genes.sssom.tsv
```

### Utility Commands

```console
# Download files only (without processing)
pysec2pri download chebi -o ./data
pysec2pri download all -o ./data

# Check for new releases
pysec2pri check-release chebi
pysec2pri check-release all

# Compare two SSSOM files
pysec2pri diff old.sssom.tsv new.sssom.tsv
```

## Python API

```python
from pysec2pri import parse_chebi, write_sssom

# Parse source file
mapping_set = parse_chebi("ChEBI_complete_3star.sdf")

# Write SSSOM output
write_sssom(mapping_set, "chebi_sec2pri.sssom.tsv")
```

### Download and Parse

```python
from pysec2pri import parse_chebi, write_sssom
from pysec2pri.download import download_datasource

# Download latest files
files = download_datasource("chebi", output_dir="./data")

# Parse downloaded file
mapping_set = parse_chebi(files["sdf"])
write_sssom(mapping_set, "chebi_sec2pri.sssom.tsv")
```

### Compare Releases

```python
from pysec2pri import parse_chebi
from pysec2pri.diff import diff_mapping_sets, summarize_diff

old_set = parse_chebi("chebi_v220.sdf")
new_set = parse_chebi("chebi_v221.sdf")

diff = diff_mapping_sets(old_set, new_set)
print(summarize_diff(diff))
```

### Supported Functions

| Function               | Description                        |
| ---------------------- | ---------------------------------- |
| `parse_chebi()`        | Parse ChEBI SDF file               |
| `parse_hmdb()`         | Parse HMDB XML/ZIP file            |
| `parse_hgnc()`         | Parse HGNC withdrawn + complete    |
| `parse_ncbi()`         | Parse NCBI gene_history + gene_info|
| `parse_uniprot()`      | Parse UniProt sec_ac + delac files |
| `write_sssom()`        | Write MappingSet to SSSOM TSV      |
| `download_datasource()`| Download source files              |
| `check_release()`      | Check for new upstream releases    |
| `diff_mapping_sets()`  | Compare two MappingSets            |
| `diff_sssom_files()`   | Compare two SSSOM files            |

## Documentation

Full documentation: <https://pysec2pri.readthedocs.io/>

## License

MIT License. See [LICENSE](LICENSE) for details.
