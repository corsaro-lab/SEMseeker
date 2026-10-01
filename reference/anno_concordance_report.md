# Cross-pipeline annotation concordance report (Illumina manifest vs WGBS)

Validates
[`anno_area_granges_build`](https://corsaro-lab.github.io/SEMseeker/reference/anno_area_granges_build.md)
(the WGBS / long-read annotation path) against the Illumina manifest
annotation used by
[`anno_probe_annotation_build`](https://corsaro-lab.github.io/SEMseeker/reference/anno_probe_annotation_build.md).
The two pipelines are expected to assign CpGs to the same semantic areas
(`GENE_*`, `ISLAND_*`, `CHR_CYTOBAND`, `DMR_*`). This function takes a
subset of Illumina probes of known annotation, runs the WGBS pipeline on
the same coordinates, and reports the concordance rate per area.

## Usage

``` r
anno_concordance_report(
  tech = "K850",
  n_probes = 1000L,
  genome_build = "hg19",
  areas = NULL,
  seed = 42L,
  csv_out = NULL
)
```

## Arguments

- tech:

  Character: one of `"K27"`, `"K450"`, `"K850"`. Default `"K850"`.

- n_probes:

  Integer (or `NULL` for all probes). Default 1000L.

- genome_build:

  Character: `"hg19"` (default), `"hg38"`, or `"mm10"`.

- areas:

  Character vector of `area_subarea` names to benchmark (e.g.
  `"GENE_BODY"`, `"ISLAND_N_SHORE"`). `NULL` (default) = all supported
  areas.

- seed:

  Integer, random seed for probe subsampling. Default 42L.

- csv_out:

  Optional path where the report is written via
  [`write.csv2()`](https://rdrr.io/r/utils/write.table.html). `NULL`
  (default) = no file.

## Value

A `data.frame` with one row per area and columns: `area`, `category`,
`n_probes`, `n_both_na`, `n_both_labeled`, `n_only_illumina`,
`n_only_wgbs`, `n_label_match_strict`, `n_label_match_intersection`,
`concordance_rate_strict`, `concordance_rate_intersection`.

## Details

Concordance categories (informational, returned in the `category`
column):

- `"bundled"` — both pipelines use the same bundled data
  (`cytoband_hg19.rda`, `dmr_annotation.rda`). Expected rate: 1.0.

- `"txdb"` — WGBS path uses `TxDb` UCSC `knownGene` while the Illumina
  manifest uses RefSeq. Some probes map to different gene symbols
  between the two sources.

- `"annotationhub"` — WGBS path uses CpG island BED from `AnnotationHub`
  (UCSC). Expected to be near-identical to the manifest but boundary
  probes may shift across shore/shelf strata.

Two concordance rates are returned:

- `concordance_rate_strict`: label sets (split on `;`) must match
  exactly (`setequal`).

- `concordance_rate_intersection`: non-empty set intersection (at least
  one shared gene / island between the two pipelines).

NA-only probes (missing in both pipelines) are counted in `n_both_na`
and excluded from the rate denominators to avoid inflating the score.

## See also

[`anno_probe_annotation_build`](https://corsaro-lab.github.io/SEMseeker/reference/anno_probe_annotation_build.md),
[`anno_area_granges_build`](https://corsaro-lab.github.io/SEMseeker/reference/anno_area_granges_build.md).
