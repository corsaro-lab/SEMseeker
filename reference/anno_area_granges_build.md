# Build a GRanges object for a given genomic area/subarea

Constructs the same region boundaries used by Illumina array annotation
(TSS200, TSS1500, gene body, CpG islands, shores, shelves ...) from TxDb
packages and AnnotationHub, so that WGBS and long-read analyses share
identical region semantics with Illumina array analyses.

## Usage

``` r
anno_area_granges_build(area_subarea, genome_build = NULL)
```

## Arguments

- area_subarea:

  Character scalar: area and subarea joined by `"_"` (e.g.
  `"GENE_BODY"`, `"ISLAND_N_SHORE"`, `"CHR_CYTOBAND"`).

- genome_build:

  Character scalar: reference assembly. Defaults to `ssEnv$genome_build`
  (set by
  [`core_init_env`](https://corsaro-lab.github.io/SEMseeker/reference/core_init_env.md)),
  or `"hg19"` if the session is not initialised.

## Value

A `GRanges` object. `mcols(gr)$label` contains the subarea identifier
used downstream to group CpGs (gene symbol, island coordinate, cytoband
name, etc.).

## Details

Results are cached in memory for the duration of the R session to avoid
repeated package loads. CpG island tracks downloaded from AnnotationHub
are also cached on disk in `tools::R_user_dir("SEMseeker", "cache")`.

## PROBE_WHOLE semantics by technology

`PROBE_WHOLE` is **not** handled by this function. It is resolved inline
by
[`anno_probe_features_get`](https://corsaro-lab.github.io/SEMseeker/reference/anno_probe_features_get.md):

- **Illumina**: one row per array probe (manufacturer ID, e.g.
  `cg00000029`). Probe identity is meaningful and cross-study comparable
  for the same array platform.

- **WGBS / LONGREAD**: treated as `POSITION_WHOLE` - one row per genomic
  position encoded as `"CHR\_START"` (e.g. `"1\_10000"`). Cross-study
  comparisons require the same `genome_build`.

## Required packages

Install via
[`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html).

- GENE areas:

  `TxDb.Hsapiens.UCSC.hg19.knownGene` (or hg38/mm10), `GenomicFeatures`,
  `GenomicRanges`, `IRanges`. `org.Hs.eg.db` is required in practice:
  without it the label of a gene region stays the Entrez id that TxDb
  keys on, while the Illumina path labels with gene symbols, so the two
  backends name the same gene two different ways. The code falls back
  rather than failing, but the result is not comparable across backends.

- ISLAND areas:

  `AnnotationHub` (downloads track on first use, then caches locally).

- CHR / DMR areas:

  no extra packages needed.
