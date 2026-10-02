# Compare enrichment results across studies

Compares enrichment (pathway/term) results across two or more
independent studies stored under \`result_folder\` and produces Venn
diagrams of the overlapping enriched terms. Part of the cross-study
replication workflow, alongside \[meta_association_overlaps_studies()\].

## Usage

``` r
meta_enrichment_compare_studies(result_folder, ...)
```

## Arguments

- result_folder:

  Path to the folder holding the per-study enrichment results to be
  compared.

- ...:

  Additional arguments forwarded to \[core_init_env()\] (e.g.
  \`maxResources\`, \`parallel_strategy\`).

## Value

Invisibly \`NULL\`; Venn diagrams and comparison tables are written
under \`result_folder\`.

## See also

\[meta_association_overlaps_studies()\],
\[meta_enrichment_overlaps_subsamples()\]

## Examples

``` r
# See vignette("pathway-analysis", package = "SEMseeker") for a runnable
# cross-study enrichment comparison workflow on real result folders.
invisible(NULL)
```
