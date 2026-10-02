# Enrichment stability across random subsamples of a study

Tests the stability of enrichment (pathway/term) results under random
subsampling of a single cohort, quantifying how robust the enriched
terms are to sample variability. The intra-study counterpart of
\[meta_enrichment_compare_studies()\].

## Usage

``` r
meta_enrichment_overlaps_subsamples(
  inference_details,
  pathways_sql_selection = "",
  old_label_samples_sql_condition = "",
  new_label_samples_sql_condition = "",
  old_label_association_results_sql_condition = "",
  new_label_association_results_sql_condition = "",
  run_prefix = "",
  pathway_package = "",
  association_pvalue_column = "",
  significance = TRUE,
  result_folder,
  ...
)
```

## Arguments

- inference_details:

  Data frame describing the enrichment analyses to run (one row per
  analysis).

- pathways_sql_selection:

  Optional SQL WHERE fragment to restrict the enriched terms considered.

- old_label_samples_sql_condition, new_label_samples_sql_condition:

  SQL conditions selecting the two sample sets to contrast.

- old_label_association_results_sql_condition,
  new_label_association_results_sql_condition:

  SQL conditions selecting the association results to contrast.

- run_prefix:

  Optional prefix for the run outputs.

- pathway_package:

  Enrichment backend to use (e.g. \`"WebGestalt"\`); \`""\` compares all
  available backends.

- association_pvalue_column:

  Name of the p-value column driving significance.

- significance:

  Logical; keep only significant terms when \`TRUE\`.

- result_folder:

  Path to the folder holding the study results.

- ...:

  Additional arguments forwarded to \[core_init_env()\].

## Value

Invisibly \`NULL\`; stability tables and plots are written under
\`result_folder\`.

## See also

\[meta_enrichment_compare_studies()\]

## Examples

``` r
# See vignette("pathway-analysis", package = "SEMseeker") for a runnable
# enrichment-stability (subsampling) workflow on real result folders.
invisible(NULL)
```
