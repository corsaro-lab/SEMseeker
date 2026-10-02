# test_samplesheet_gse133774

Sample sheet for
[`test_signal_gse133774`](https://corsaro-lab.github.io/SEMseeker/reference/test_signal_gse133774.md)
in the canonical SEMseeker three-class design (Reference / Control /
Case).

## Usage

``` r
test_samplesheet_gse133774
```

## Format

A data frame with 16 rows and 2 columns:

- Sample_ID:

  Sample identifier matching column names of
  [`test_signal_gse133774`](https://corsaro-lab.github.io/SEMseeker/reference/test_signal_gse133774.md)
  (e.g. `"CTRL01"`, `"L1"`).

- Sample_Group:

  One of `"Reference"`, `"Control"`, or `"Case"`.

## Source

Derived from GEO accession GSE133774 metadata.

## Details

Control samples (CTRL01-CTRL06) appear twice: once as `Reference`
(population baseline for IQR threshold estimation) and once as `Control`
(comparison group). Family samples (L1-L4) are `Case`; L1 is the BWS
proband. This is the *Reference-reuse pattern* documented in the
getting-started vignette.
