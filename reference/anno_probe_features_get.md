# Retrieve probe feature annotations for a given genomic area

Returns a `data.frame` of CpG probe coordinates and feature annotations
for the requested area/subarea combination.

## Usage

``` r
anno_probe_features_get(area_subarea)
```

## Arguments

- area_subarea:

  Character scalar: area and subarea joined by an underscore (e.g.
  `"GENE_BODY"`, `"CHR_WHOLE"`, `"ISLAND_N_SHORE"`). If no underscore is
  present, `"_WHOLE"` is appended automatically.

## Value

A `data.frame` with columns `PROBE`, `CHR`, `START`, `END`, and the
requested feature column.

## Details

For **Illumina** data (K850/K450/K27), annotations are built from
Bioconductor array annotation packages (see
[`anno_probe_annotation_build`](https://corsaro-lab.github.io/SEMseeker/reference/anno_probe_annotation_build.md))
and cached in the session environment.

For **WGBS** and **LONGREAD** data, coordinates are read directly from
the saved POSITION pivot parquet; semantic areas (GENE\_\*, ISLAND\_\*,
CHR_CYTOBAND, DMR\_\*) are resolved via
[`anno_area_granges_build`](https://corsaro-lab.github.io/SEMseeker/reference/anno_area_granges_build.md)
and
[`GenomicRanges::findOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html).

## PROBE_WHOLE vs POSITION_WHOLE - technology semantics

The area `PROBE_WHOLE` has different meanings depending on technology:

- Illumina:

  Each row identifies a specific *array probe* by its manufacturer ID
  (e.g. `cg00000029`). The statistical test is performed at the
  individual-probe level. Probe identity is meaningful here: two studies
  using the same array share the exact same set of probe IDs and can be
  directly compared.

- WGBS / LONGREAD:

  There are no probe IDs. `PROBE_WHOLE` is treated as
  **`POSITION_WHOLE`**: each row identifies a CpG by its genomic
  coordinate (`CHR\_START`, e.g. `"1\_10000"`). The statistical test is
  performed at the individual-position level. Two WGBS datasets can be
  compared only if they share the same reference genome
  (`ssEnv\$genome_build`) - mismatches are detected by the session
  provenance guard (C-06).

In both cases the downstream analysis pipeline is identical; the
distinction is purely in what the `PROBE` column *means* to the
researcher.

Probes on sex chromosomes are removed when `ssEnv$sex_chromosome_remove`
is `TRUE`.
