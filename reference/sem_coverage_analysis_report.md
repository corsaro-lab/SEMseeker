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
  parallel_strategy = "multisession",
  ...
)
```

## Arguments

- signal_data:

  character. Path to the signal parquet file under
  `Data/Pivots/SIGNAL/`, or a data.frame already loaded into memory. The
  file is named `SIGNAL_<FIGURE>_PROBE_WHOLE_<GENOME_BUILD>.parquet`,
  where `FIGURE` is the scale the run was written on - `BETA` for a
  proportion bounded in \[0,1\], `MVALUE` for the logit-transformed one.

- result_folder:

  character. Path to the SEMseeker result folder.

- maxResources:

  numeric. Maximum percentage of CPU cores to use (default 90).

- parallel_strategy:

  character. Parallelisation backend passed to `future` (default
  `"multisession"`). Asking for `"multicore"` is accepted and converted
  to `"multisession"`: it means fork(), which is unsafe with this
  package's native thread pool on every platform that offers it, and
  absent on Windows. The conversion is logged.

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
# BETA is the scale this particular run was written on; a run on the
# logit scale produces SIGNAL_MVALUE_PROBE_WHOLE_HG19.parquet instead.
sem_coverage_analysis_report(
  signal_data   = "~/semseeker_results/Data/Pivots/SIGNAL/SIGNAL_BETA_PROBE_WHOLE_HG19.parquet",
  result_folder = "~/semseeker_results/"
)
} # }
```
