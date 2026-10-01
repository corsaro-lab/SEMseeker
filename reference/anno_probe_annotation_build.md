# Build the Illumina probe annotation table

Internal helper. Assembles the per-probe annotation (genomic position,
cytoband and DMR/area membership) for an Illumina methylation array
platform, joining the bundled
[`cytoband_hg19`](https://corsaro-lab.github.io/SEMseeker/reference/cytoband_hg19.md)
and
[`dmr_annotation`](https://corsaro-lab.github.io/SEMseeker/reference/dmr_annotation.md)
reference data.

## Usage

``` r
anno_probe_annotation_build(tech, force = FALSE)
```

## Arguments

- tech:

  Illumina platform identifier (e.g. "EPIC", "450k", "27k").

- force:

  Logical; rebuild even when a cached annotation is available.

## Value

A data frame of per-probe annotation columns.
