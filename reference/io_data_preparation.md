# Title

Title

## Usage

``` r
io_data_preparation(
  family_test,
  transformation_y,
  tempDataFrame,
  independent_variable,
  g_start,
  g_end,
  covariates,
  key,
  transformation_x = "none"
)
```

## Arguments

- family_test:

  test or regression to apply

- transformation_y:

  transformation_y to apply to data

- tempDataFrame:

  data frame to use for test/regression

- independent_variable:

  regressor

- g_start:

  index of the first burden column in tempDataFrame

- g_end:

  index of the last burden column in tempDataFrame

- covariates:

  vector of covariates to be found in the sample sheet

- key:

  named list with AREA, SUBAREA, MARKER and FIGURE identifiers, used to
  name the artefact in the log when degenerate columns are dropped

- transformation_x:

  transformation to apply to the independent variable before the fit;
  "none" leaves it untouched

## Value

A named list with two elements: `tempDataFrame` (the prepared and
optionally transformed data.frame) and `independent_variableLevels` (the
factor levels of the independent variable, or `NULL` for continuous
outcomes).
