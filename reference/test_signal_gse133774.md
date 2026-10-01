# test_signal_gse133774

Real EPIC methylation beta-value matrix from GEO series GSE133774
(Infinium MethylationEPIC 850k, GPL21145), filtered to the probe IDs
present in
[`test_master_features`](https://corsaro-lab.github.io/SEMseeker/reference/test_master_features.md).
Contains 10 samples from a Beckwith-Wiedemann Syndrome (BWS) /
Multi-Locus Imprinting Disturbance (MLID) family study: 6 unrelated
controls and 4 family members (L1 = BWS proband, L2–L4 = siblings/parent
with NLRP5 compound heterozygous variants).

## Usage

``` r
test_signal_gse133774
```

## Format

A numeric matrix with 18,089 rows (EPIC probe IDs) and 10 columns
(samples: CTRL01–CTRL06, L1–L4). Values are beta coefficients in \\\[0,
1\]\\; a small number of probes may have `NA` (QC-filtered positions in
the original GEO submission).

## Source

GEO accession GSE133774, series matrix file parsed by
`data-raw/build_test_signal_fixture.R`. Original study: Docherty *et
al.* (2020), NLRP5 variants associated with MLID and BWS.

## Details

This fixture replaces all
[`rbeta()`](https://rdrr.io/r/stats/Beta.html)-based synthetic signal
generation in the test suite and vignette (AI-123). Because the data
contain real BWS epimutations at imprinting DMRs (KCNQ1OT1, H19/IGF2,
MEG3, ...), running the full SEMseeker pipeline on this matrix detects
biologically expected hypo-epimutation events in the Case samples
without any artificial injection.
