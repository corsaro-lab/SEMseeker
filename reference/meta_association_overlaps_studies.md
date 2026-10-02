# Compare association results across studies

Reads the inference results of several studies from a shared result
folder, joins them on the area taxonomy (marker, figure, area, subarea)
and reports, per area, which studies called it significant. The
aggregated table and one table per marker are written as CSV under the
inference folder; the function is called for those files and returns
nothing.

## Usage

``` r
meta_association_overlaps_studies(
  inference_detail,
  studies,
  alpha = 0.05,
  adjust_per_area = FALSE,
  adjust_globally = FALSE,
  pvalue_column = "PVALUE_ADJ_ALL_BH",
  statistic_parameter,
  adjustment_method = "BH",
  result_folder,
  ...
)
```

## Arguments

- inference_detail:

  One row of the inference specification identifying the request whose
  results are compared. If more than one row is passed the first
  supplies the run-level parameters.

- studies:

  Character vector of study names to compare. Each must have its results
  already present under `result_folder`.

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

- result_folder:

  Folder holding the studies' results and receiving the comparison
  output.

- ...:

  Passed through to the session setup.

## Value

Called for its side effect: CSV files written under the inference
folder. Returns `NULL` invisibly, and early if no results are found.
