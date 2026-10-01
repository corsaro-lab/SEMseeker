# dmr_annotation

Differentially methylated region (DMR) annotations for CpG probes on
Illumina methylation arrays. Contains only probes that overlap at least
one known imprinted or disease-associated DMR. This lightweight table
(~1,600 rows) supplements the Bioconductor array annotation packages
used by
[`anno_probe_annotation_build`](https://corsaro-lab.github.io/SEMseeker/reference/anno_probe_annotation_build.md),
which do not carry DMR-level annotations.

## Usage

``` r
dmr_annotation
```

## Format

A data frame with three columns:

- PROBE:

  Illumina probe identifier (e.g. `"cg00000924"`).

- DMR_WHOLE:

  Genomic region label of the enclosing DMR (e.g. `"KCNQ1OT1:TSS-DMR"`).

- DMR_DMR:

  Fine-grained DMR label, identical to `DMR_WHOLE` for most probes.
