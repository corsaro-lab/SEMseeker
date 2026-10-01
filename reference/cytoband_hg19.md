# cytoband_hg19

Cytogenetic band coordinates for the hg19 human genome assembly. Used by
[`anno_probe_annotation_build`](https://corsaro-lab.github.io/SEMseeker/reference/anno_probe_annotation_build.md)
to assign a `CHR_CYTOBAND` label (e.g. `"q12.2"`) to each CpG probe
based on its chromosomal position. Band boundaries are derived from
probe positions in the full Illumina EPIC annotation and cover all
autosomes and sex chromosomes.

## Usage

``` r
cytoband_hg19
```

## Format

A data frame with 829 rows and four columns:

- CHR:

  Chromosome identifier without `"chr"` prefix (e.g. `"1"`, `"X"`).

- START:

  Approximate start position of the cytogenetic band (bp).

- END:

  Approximate end position of the cytogenetic band (bp).

- CYTOBAND:

  ISCN band label (e.g. `"q12.2"`, `"p36.33"`).
