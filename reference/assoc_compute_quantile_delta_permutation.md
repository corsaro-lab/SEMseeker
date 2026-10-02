# Compute quantile delta permutation (CPU)

Compute quantile delta permutation (CPU)

## Usage

``` r
assoc_compute_quantile_delta_permutation(
  sig.formula,
  df,
  shuffle = FALSE,
  quantile = 0.5
)
```

## Arguments

- sig.formula:

  formula to apply

- df:

  dataframe to use

- shuffle:

  logical; if TRUE, permute the independent variable before fitting

- quantile:

  quantile level (0-1) at which to compute the group difference; default
  0.5 (median)

## Value

A numeric scalar: the observed quantile difference between groups (used
as the test statistic in the permutation test).
