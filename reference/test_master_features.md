# test_master_features

Master probe-features fixture used by the test suite and by the
`getting-started` vignette as a runnable synthetic example. Contains
20,000 real Illumina EPIC (K850) probe IDs, sampled deterministically
around 58 curated human imprinting DMRs (KCNQ1OT1, H19/IGF2, MEG3/DLK1,
GNAS, PEG3, SNURF, PLAGL1, NNAT, MEST, ...).

## Usage

``` r
test_master_features
```

## Format

A data frame with 20,000 rows and 6 columns:

- PROBE:

  Illumina EPIC probe identifier (e.g. `"cg00000924"`).

- CHR:

  Chromosome label without `"chr"` prefix (e.g. `"11"`, `"X"`).

- START:

  1-based start position on hg19 (bp).

- END:

  End position; equal to `START` for array probes.

- ABSOLUTE:

  Concatenated chromosome/position key (`paste(CHR, START, sep = "_")`).

- DMR_LABEL:

  Imprinting DMR identifier (e.g. `"KCNQ1OT1:TSS-DMR"`) for the ~820
  *signal* probes; `NA` for flanking and background probes.

## Source

Built from
[`dmr_annotation`](https://corsaro-lab.github.io/SEMseeker/reference/dmr_annotation.md)
and the `IlluminaHumanMethylationEPICanno.ilm10b4.hg19` Bioconductor
annotation package with `set.seed(20210713)` (date of the v.0.1.9 Zenodo
software-archive release, DOI
[10.5281/zenodo.5095417](https://doi.org/10.5281/zenodo.5095417)).
Re-generate with `Rscript data-raw/build_test_master_features.R`.

## Details

Construction strategy (see `data-raw/build_test_master_features.R`):

1.  All unique probe IDs in
    [`dmr_annotation`](https://corsaro-lab.github.io/SEMseeker/reference/dmr_annotation.md)
    are kept as the biological *signal* layer (~820 probes labelled with
    their parent DMR in `DMR_LABEL`).

2.  Each imprinting DMR is expanded into a genomic window of \\\pm
    500\\\textrm{kb}\\ to capture the surrounding imprinted gene
    cluster; all EPIC probes within these windows are added.

3.  If the resulting pool is below 20,000, the remainder is filled with
    EPIC probes sampled deterministically from outside any imprinting
    window (background layer, `DMR_LABEL = NA`).

This makes the fixture biology-aware (BWS / SRS / PWS-AS / GNAS-related
regions are guaranteed to be represented) while preserving statistical
validity (20k probes is the minimum size used by the test suite for IQR
robustness). Because the fixture ships with real EPIC probe IDs, the
Bioconductor annotation package is not required at test-setup time –
eliminating the macOS `tcltk`/XQuartz segfault path triggered by
[`requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")`](https://bitbucket.com/kasperdanielhansen/Illumina_EPIC).
