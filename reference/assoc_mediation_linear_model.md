# Mediation analysis using linear models

Fits a causal mediation model via `mediate`, testing whether the effect
of `treatment` on `outcome` methods::is (partially) mediated by a
`mediator` variable. The formula must follow the convention
`mediator ~ outcome + treatment [+ covariates]`.

## Usage

``` r
assoc_mediation_linear_model(
  family_test,
  tempDataFrame,
  sig.formula,
  transformation_y,
  plot,
  samples_sql_condition = samples_sql_condition,
  key
)
```

## Arguments

- family_test:

  Character string encoding the model type and permutation counts,
  formatted as `"mediation_<permutations_test>_<permutations>"`.

- tempDataFrame:

  `data.frame` containing all variables referenced in `sig.formula`.

- sig.formula:

  Formula of the form `mediator ~ outcome + treatment [+ covariates]`.

- transformation_y:

  Character scalar: transformation applied to the dependent variable
  (e.g. `"none"`, `"log"`).

- plot:

  Logical; if `TRUE`, generate diagnostic plots.

- samples_sql_condition:

  Character scalar: SQL `WHERE` clause used to subset samples
  (propagated to file-naming helpers).

- key:

  Named list with elements `AREA`, `SUBAREA`, `MARKER`, and `FIGURE`
  identifying this test instance.

## Value

A `data.frame` with one row containing mediation results: ACME (average
causal mediation effect), ADE (average direct effect), total effect,
proportion mediated, their confidence intervals and p-values, plus the
number of permutations used.

## Details

The two-round permutation scheme mirrors the approach used in
`assoc_mean_permutation` and `assoc_spearman_permutation`: a fast first
round (`permutations_test`) establishes significance; if significant, a
slower second round (`permutations`) refines the estimate.
