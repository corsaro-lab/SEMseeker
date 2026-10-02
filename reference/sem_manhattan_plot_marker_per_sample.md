# Manhattan plot of SEM markers for a single sample

Manhattan plot of SEM markers for a single sample

## Usage

``` r
sem_manhattan_plot_marker_per_sample(
  sample_name = "NAME",
  probes_range = 1000:2000,
  hyper_color = "blue",
  hypo_color = "orange",
  non_outlier_color = "grey",
  limit_label_color = ssEnv$color_palette[1],
  result_folder = result_folder,
  maxResources = maxResources,
  parallel_strategy = parallel_strategy,
  ...
)
```

## Arguments

- sample_name:

  character. Sample identifier as it appears in the sample sheet
  (default `"NAME"`).

- probes_range:

  integer vector. Probe index range to display (default `1000:2000`).

- hyper_color:

  character. Colour for hypermethylated probes (default `"blue"`).

- hypo_color:

  character. Colour for hypomethylated probes (default `"orange"`).

- non_outlier_color:

  character. Colour for non-outlier probes (default `"grey"`).

- limit_label_color:

  character. Colour for threshold labels.

- result_folder:

  character. Path to the SEMseeker result folder.

- maxResources:

  numeric. Maximum percentage of CPU cores to use (default 90).

- parallel_strategy:

  character. Parallelisation backend (default `"multisession"`). Asking
  for `"multicore"` is accepted and converted to `"multisession"`: it
  means fork(), which is unsafe with this package's native thread pool
  on every platform that offers it, and absent on Windows. The
  conversion is logged.

- ...:

  Additional arguments passed to
  [`core_init_env()`](https://corsaro-lab.github.io/SEMseeker/reference/core_init_env.md).

## Value

Invisibly `NULL`. A PNG plot is saved under `Charts/MARKER_PER_SAMPLE/`
in `result_folder`.

## Examples

``` r
result_dir <- tempdir()
if (FALSE) { # \dontrun{
sem_manhattan_plot_marker_per_sample(
  sample_name  = "CASE_001",
  probes_range = 1:5000,
  result_folder = result_dir
)
} # }
```
