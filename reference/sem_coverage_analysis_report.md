# Probe coverage analysis report

Generate a coverage analysis report for Illumina methylation array data,
summarising probe representation across genomic regions. WGBS data are
not supported and will cause the function to stop with an informative
error.

## Usage

``` r
sem_coverage_analysis_report(
  signal_data,
  result_folder,
  maxResources = 90,
  parallel_strategy = "multicore",
  ...
)
```

## Arguments

- signal_data:

  character. Path to the signal parquet file (e.g.
  `Data/SIGNAL_MEAN_PROBE_WHOLE.parquet`) or a data.frame already loaded
  into memory.

- result_folder:

  character. Path to the SEMseeker result folder.

- maxResources:

  numeric. Maximum percentage of CPU cores to use (default 90).

- parallel_strategy:

  character. Parallelisation backend passed to `future` (default
  `"multicore"`).

- ...:

  Additional named arguments passed to
  [`core_init_env()`](https://corsaro-lab.github.io/SEMseeker/reference/core_init_env.md).

## Value

Invisibly `NULL`. Coverage tables and charts are written to the result
folder.

## Examples

``` r
result_dir <- tempdir()
if (FALSE) { # \dontrun{
sem_coverage_analysis_report(
  signal_data   = "~/semseeker_results/Data/SIGNAL_MEAN_PROBE_WHOLE.parquet",
  result_folder = "~/semseeker_results/"
)
} # }
```
