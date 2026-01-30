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

| Datasource | license                                                                                                                            | citation                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| ---------- | ---------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| ChEBI      | [CC BY 4.0](docs/licenses/chebi/LICENSE).                                                                                          | Hastings J, Owen G, Dekker A, et al. ChEBI in 2016: Improved services and an expanding collection of metabolites. Nucleic Acids Research. 2016 Jan;44(D1):D1214-9. DOI: [10.1093/nar/gkv1031](https://doi.org/10.1093/nar/gkv1031). PMID: 26467479; PMCID: PMC4702775.                                                                                                                                                                                                                                                                                                                                        |
| HMDB       | [CC0](https://hmdb.ca/about#compliance:~:text=international%20scientific%20conferences.-,Citing%20the%20HMDB,-HMDB%20is%20offered) | Wishart DS, Guo A, Oler E, Wang F, Anjum A, Peters H, Dizon R, Sayeeda Z, Tian S, Lee BL, Berjanskii M, Mah R, Yamamoto M, Jovel J, Torres-Calzada C, Hiebert-Giesbrecht M, Lui VW, Varshavi D, Varshavi D, Allen D, Arndt D, Khetarpal N, Sivakumaran A, Harford K, Sanford S, Yee K, Cao X, Budinski Z, Liigand J, Zhang L, Zheng J, Mandal R, Karu N, Dambrova M, Schiöth HB, Greiner R, Gautam V. HMDB 5.0: the Human Metabolome Database for 2022. Nucleic Acids Res. 2022 Jan 7;50(D1):D622-D631. doi: [10.1093/nar/gkab1062](https://doi.org/10.1093/nar/gkab1062). PMID: 34986597; PMCID: PMC8728138. |
| HGNC       | [link](https://www.genenames.org/about/license/)                                                                                   | Seal RL, Braschi B, Gray K, Jones TEM, Tweedie S, Haim-Vilmovsky L, Bruford EA. Genenames.org: the HGNC resources in 2023. Nucleic Acids Res. 2023 Jan 6;51(D1):D1003-D1009. doi: [10.1093/nar/gkac888](https://doi.org/10.1093/nar/gkac888). PMID: 36243972; PMCID: PMC9825485.                                                                                                                                                                                                                                                                                                                              |
| NCBI       | [link](https://www.ncbi.nlm.nih.gov/home/about/policies/)                                                                          | Sayers EW, Bolton EE, Brister JR, Canese K, Chan J, Comeau DC, Connor R, Funk K, Kelly C, Kim S, Madej T, Marchler-Bauer A, Lanczycki C, Lathrop S, Lu Z, Thibaud-Nissen F, Murphy T, Phan L, Skripchenko Y, Tse T, Wang J, Williams R, Trawick BW, Pruitt KD, Sherry ST. Database resources of the national center for biotechnology information. Nucleic Acids Res. 2022 Jan 7;50(D1):D20-D26. doi: [10.1093/nar/gkab1112](https://doi.org/10.1093/nar/gkab1112). PMID: 34850941; PMCID: PMC8728269.                                                                                                        |
| UniProt    | [CC BY 4.0](https://ftp.uniprot.org/pub/docs/licenses/uniprot/current_release/knowledgebase/complete/LICENSE)                      | UniProt Consortium. UniProt: the universal protein knowledgebase in 2021. Nucleic Acids Res. 2021 Jan 8;49(D1):D480-D489. doi: [10.1093/nar/gkaa1100](https://doi.org/10.1093/nar/gkaa1100). PMID: 33237286; PMCID: PMC7778908.                                                                                                                                                                                                                                                                                                                                                                               |
| Wikidata   |                                                                                                                                    | Vrandecic, D., Krotzsch, M. Wikidata: a free collaborative knowledgebase. Communications of the ACM. 2014. doi: [10.1145/2629489](https://doi.org/10.1145/2629489).                                                                                                                                                                                                                                                                                                                                                                                                                                           |

## Installation

```console
uv pip install pysec2pri
```

Or install from source:

```console
uv pip install git+https://github.com/sec2pri/pysec2pri.git
```

## Quick Start

### Command Line Usage

To obtain the secondary to primary identifier SSSOM mapping set for ChEBI:

```bash
pysec2pri chebi
```

This will automatically download the latest ChEBI release and generate an SSSOM
mapping file in your current directory.

To process locally and specify the output:

```bash
pysec2pri chebi ChEBI_complete_3star.sdf --output my_mappings.sssom.tsv
```

For more options and help on any command:

```bash
pysec2pri --help
pysec2pri chebi --help
```

### SSSOM Output

The default output is in
[SSSOM](https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://mapping-commons.github.io/sssom/&ved=2ahUKEwiwwcOYjLKSAxXR2wIHHXXRB4gQFnoECA4QAQ&usg=AOvVaw3eQLVYhT0H9PXyznzUSExF)
(Simple Standard for Sharing Ontology Mappings) TSV format.

## Documentation

Full documentation: <https://pysec2pri.readthedocs.io/>

## License

MIT License. See [LICENSE](LICENSE) for details.
