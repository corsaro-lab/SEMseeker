# Bayesian posterior probability analysis of SEMseeker mutations and lesions

Computes P(case \| epimutated) via empirical Bayes for each marker and
genomic area, providing a probabilistic complement to the frequentist
association analysis.

## Usage

``` r
assoc_bayes_analysis(
  result_folder,
  independent_variable = "Sample_Group",
  maxResources = 90,
  parallel_strategy = "multisession",
  bayes_case_threshold = 0.9,
  bayes_control_threshold = 0.1,
  ...
)
```

## Arguments

- result_folder:

  character. Path to the SEMseeker result folder.

- independent_variable:

  character. Sample sheet column defining the case/control grouping
  variable (default `"Sample_Group"`).

- maxResources:

  numeric. Maximum percentage of CPU cores to use (default 90).

- parallel_strategy:

  character. Parallelisation backend; one of `"multisession"`,
  `"sequential"`, `"cluster"` (default `"multisession"`). Asking for
  `"multicore"` is accepted and converted to `"multisession"`: it means
  fork(), which is unsafe with this package's native thread pool on
  every platform that offers it, and absent on Windows. The conversion
  is logged.

- bayes_case_threshold:

  numeric. Minimum P(case \| epimutated) required to report a hit
  (default 0.9).

- bayes_control_threshold:

  numeric. Maximum P(control \| epimutated) allowed to report a hit
  (default 0.1).

- ...:

  Additional arguments passed to
  [`core_init_env()`](https://corsaro-lab.github.io/SEMseeker/reference/core_init_env.md).

## Value

Invisibly `NULL`. Bayesian posterior probability tables
(`bayes_analysis_*.csv`) are written to the `Euristic/` sub-folder of
`result_folder`.

## Examples

``` r
result_dir <- tempdir()
if (FALSE) { # \dontrun{
assoc_bayes_analysis(
  result_folder        = "~/semseeker_results/",
  independent_variable = "Sample_Group"
)
} # }
```
