# Diagnostic performance (sensitivity & specificity) of SEMseeker mutations and lesions

Diagnostic performance (sensitivity & specificity) of SEMseeker
mutations and lesions

## Usage

``` r
diagnostic_performance(
  samples_sql_selection = "",
  combinations,
  result_folder,
  independent_variable = "Sample_Group",
  maxResources = 90,
  parallel_strategy = "multisession",
  ...
)
```

## Arguments

- samples_sql_selection:

  character. SQL-style filter expression to restrict which samples are
  included (default `""`, all samples).

- combinations:

  list of character vectors. Each element is a vector of sample group
  labels to compare (e.g. `list(c("CASE_A", "CASE_B"))`).

- result_folder:

  character. Path to the SEMseeker result folder.

- independent_variable:

  character. Sample sheet column defining the grouping variable (default
  `"Sample_Group"`).

- maxResources:

  numeric. Maximum percentage of CPU cores to use (default 90).

- parallel_strategy:

  character. Parallelisation backend; one of `"multisession"`,
  `"sequential"`, `"cluster"` (default `"multisession"`). Asking for
  `"multicore"` is accepted and converted to `"multisession"`: it means
  fork(), which is unsafe with this package's native thread pool on
  every platform that offers it, and absent on Windows. The conversion
  is logged.

- ...:

  Additional arguments passed to
  [`core_init_env()`](https://corsaro-lab.github.io/SEMseeker/reference/core_init_env.md).

## Value

Invisibly `NULL`. Sensitivity / specificity tables are written to the
`Euristic/` sub-folder of `result_folder`.

## Examples

``` r
result_dir <- tempdir()
if (FALSE) { # \dontrun{
diagnostic_performance(
  combinations         = list(c("CASE_A", "CASE_B")),
  result_folder        = "~/semseeker_results/",
  independent_variable = "Sample_Group"
)
} # }
```
