# Compare association results across subsamples of one study

Same comparison as
[`meta_association_overlaps_studies`](https://corsaro-lab.github.io/SEMseeker/reference/meta_association_overlaps_studies.md)
but within a single study, across the subsamples produced by a
replication run. It joins the per-subsample inference results on the
area taxonomy and writes, per area, which subsamples called it
significant.

## Usage

``` r
meta_association_overlaps_subsamples(
  inference_details,
  alpha = 0.05,
  adjust_per_area = FALSE,
  adjust_globally = FALSE,
  pvalue_column = "PVALUE_ADJ_ALL_BH",
  statistic_parameter,
  adjustment_method = "BH",
  old_label = NULL,
  new_label = NULL,
  run_prefix = "",
  result_folder,
  ...
)
```

## Arguments

- inference_details:

  Inference specification rows identifying the requests whose results
  are compared.

- alpha:

  Significance threshold applied to `pvalue_column`.

- adjust_per_area:

  Adjust p-values within each area separately.

- adjust_globally:

  Adjust p-values across the whole result set.

- pvalue_column:

  Name of the p-value column to test against `alpha`. Defaults to the
  all-scope BH-adjusted column.

- statistic_parameter:

  Name of the effect-size column carried into the comparison table
  alongside the p-value.

- adjustment_method:

  Multiple-testing correction passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

- old_label, new_label:

  Optional relabelling of the subsample names in the output table.

- run_prefix:

  Prefix prepended to the output file names, to keep the results of
  several runs side by side.

- result_folder:

  Folder holding the subsamples' results and receiving the comparison
  output.

- ...:

  Passed through to the session setup.

## Value

Called for its side effect: CSV files written under the inference
folder. Returns `NULL` invisibly, and early if no results are found.
