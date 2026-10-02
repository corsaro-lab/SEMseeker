# Title

Title

## Usage

``` r
assoc_apply_stat_model(
  tempDataFrame,
  g_start,
  family_test,
  covariates = NULL,
  key,
  transformation_y,
  session_folder,
  independent_variable,
  samples_sql_condition,
  inference_detail = NULL,
  ...
)
```

## Arguments

- tempDataFrame:

  data frame to apply association

- g_start:

  index of starting data

- family_test:

  family of test to run

- covariates:

  vector of covariates

- key:

  key to identify file to elaborate

- transformation_y:

  transformation_y to apply to covariates, burden and independent
  variable

- session_folder:

  where to save log file

- independent_variable:

  independent variable name

- samples_sql_condition:

  SQL condition string to filter samples

- inference_detail:

  one row of the inference specification, carrying the per-request
  fields (scope, aggregation, transformation_x) that the model call
  needs; NULL falls back to the defaults

- ...:

  extra parameters

## Value

A data.frame with one row per tested genomic area, including columns for
p-value, adjusted p-value, test statistic, AIC, residuals, and model
metadata; returns `NULL` if no results could be computed.
