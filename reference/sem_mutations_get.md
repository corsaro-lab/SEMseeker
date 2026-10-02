# sem_mutations_get

sem_mutations_get

## Usage

``` r
sem_mutations_get(values, figure, thresholds, sampleName)
```

## Arguments

- values:

  values of methylation - data.frame with columns CHR, START, END and a
  fourth numeric VALUE column.

- figure:

  figure to get Mutations of HYPO or HYPER methylation

- thresholds:

  threshold to use for comparison - data.frame with columns CHR, START,
  END, signal_inferior_thresholds, signal_superior_thresholds.

- sampleName:

  name of the sample

## Value

mutations data.frame with columns CHR, START, END, MUTATIONS (0/1). Only
positions present in BOTH values and thresholds are returned (inner join
on CHR, START, END via util_join_values_to_thresholds).
