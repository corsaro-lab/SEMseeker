# Calculate stochastic epi mutations from a methylation dataset as outcome report of pivot

Calculate stochastic epi mutations from a methylation dataset as outcome
report of pivot

## Usage

``` r
sem_analyze_population(
  signal_data,
  sample_sheet,
  signal_thresholds,
  probe_features
)
```

## Arguments

- signal_data:

  whole matrix of data to analyze.

- sample_sheet:

  name of samplesheet's column to use as control population selector
  followed by selection value,

- signal_thresholds:

  thresholds defined to calculate epimutations

- probe_features:

  probe_features detail from 27 to EPIC illumina dataset

## Value

files into the result folder with pivot table and bedgraph. A BANNER is
logged once per batch before the per-sample loop showing:
input_positions, beta_range_positions, covered_by_inner_join - allows
immediate audit of cross-run coverage (e.g. Nanopore sample vs Illumina
reference).
