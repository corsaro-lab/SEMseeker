# Intra-study association (focused replication across models)

Re-run association analysis focusing on genomic regions that are
statistically significant in a reference study, then merge results
across two datasets for comparative inference.

## Usage

``` r
meta_association_replication(
  inference_details_origin,
  inference_details,
  result_folder,
  maxResources = 90,
  parallel_strategy = "multisession",
  start_fresh = FALSE,
  ...
)
```

## Arguments

- inference_details_origin:

  data.frame. Inference parameters for the reference (origin) study
  whose significant regions define the search space.

- inference_details:

  data.frame. Inference parameters for the target study on which the
  focused association analysis is performed.

- result_folder:

  character. Path to the SEMseeker result folder.

- maxResources:

  numeric. Maximum percentage of CPU cores to use (default 90).

- parallel_strategy:

  character. Parallelisation backend passed to `future`; one of
  `"multisession"`, `"sequential"`, `"cluster"` (default
  `"multisession"`). Asking for `"multicore"` is accepted and converted
  to `"multisession"`: it means fork(), which is unsafe with this
  package's native thread pool on every platform that offers it, and
  absent on Windows. The conversion is logged.

- start_fresh:

  logical. If `TRUE`, delete previous results before running (default
  `FALSE`).

- ...:

  Additional named arguments passed to
  [`core_init_env()`](https://corsaro-lab.github.io/SEMseeker/reference/core_init_env.md).

## Value

Invisibly `NULL`. Results are written to the inference sub-folder of
`result_folder`.

## Examples

``` r
result_dir <- tempdir()
if (FALSE) { # \dontrun{
meta_association_replication(
  inference_details_origin = inference_study1,
  inference_details        = inference_study2,
  result_folder            = "~/semseeker_comparison/"
)
} # }
```
