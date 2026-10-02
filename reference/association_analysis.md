# Association analysis of SEMseeker results

Run statistical association models between SEM metrics and a phenotype
variable. Supports group tests (Wilcoxon, t-test), GLM families
(gaussian, poisson, binomial), quantile regression, correlations
(Pearson, Kendall, Spearman), and multi-covariate formulas (e.g.
`MUTATIONS_* ~ covariate1 + covariate2`).

## Usage

``` r
association_analysis(
  inference_details,
  result_folder,
  maxResources = 90,
  parallel_strategy = "multisession",
  start_fresh = FALSE,
  ...
)
```

## Arguments

- inference_details:

  data.frame. Each row defines one analysis run. Required columns:

  independent_variable

  :   Sample sheet column used as grouping / covariate variable.

  family_test

  :   Statistical model: `"wilcoxon"`, `"stats::t.test"`, `"gaussian"`,
      `"poisson"`, `"binomial"`, `"pearson"`, `"kendall"`, `"spearman"`,
      or quantile regression as `"quantreg_<tau>_<runs>"` (e.g.
      `"quantreg_0.25_2000"`).

  transformation_y

  :   Transformation applied to the dependent variable: `"none"`,
      `"scale"`, `"log"`, `"log2"`, `"log10"`, `"exp"`, or
      `"quantile_<n>"` (e.g. `"quantile_3"`).

  scope

  :   Required. Which of the two aggregations to run, over the region
      classes of the call (`areas`, `subareas`). `"SAMPLE"` reduces the
      positions of each class to one number per sample: the burden, or
      its density, or a descriptor of the signal. `"INSTANCE"` reduces
      them to one number per instance of the class: one row per gene,
      per island, per cytoband, per probe.

      The two are **mutually exclusive**. They used to be produced
      together and written into the same file, which made every result
      the union of two different questions; a request that wants both
      now writes two rows, and each carries its own `SCOPE` in the
      output.

      It has no default on purpose. Any default would answer half of
      what was asked and say nothing about the other half.

  aggregation

  :   Required. How the positions are reduced to the one number the
      model is fitted on: `"SUM"`, `"MEAN"`, `"MEDIAN"`, `"VARIANCE"`,
      `"IQR"`, `"MODELOW"` or `"MODEHIGH"`. While every marker admitted
      exactly one operator this could stay implicit; a scope now carries
      several, so the request has to name the one it wants. Which are
      admissible depends on the marker: a count carries `SUM` (the
      burden) and `MEAN` (the density, the form comparable across
      regions of different size), while its median and IQR are
      degenerate; the two modes exist only for the signal on the beta
      scale. A request no marker of the run admits is dropped with a
      warning naming it, not answered with an empty result.

  The markers are **not** named in `inference_details` either. They are
  chosen by the `markers` argument of this call and every row of the
  request is tested on all of them, one result file per marker. A
  `marker` column was documented here for a long time and never existed
  in the vocabulary, so a request written from that description was
  rejected as carrying an unknown column.

  The region classes are **not** named in `inference_details`: they are
  the `(AREA, SUBAREA)` pairs of the run, declared with the `areas` and
  `subareas` arguments of this call and built at runtime. Both scopes
  range over the same pairs, so `scope = "SAMPLE"` with
  `areas = c("GENE", "PROBE")` gives the burden over the gene probes and
  the burden over the whole sample, and `scope = "INSTANCE"` gives one
  row per gene and one row per probe. The artefact is built on the way
  in if it does not exist yet, so a class no previous run foresaw costs
  one scan of the position pivot rather than a rerun.

  Result rows carry the coordinates as columns: `MARKER`, `FIGURE`,
  `SCOPE`, `AREA`, `SUBAREA`, `AGGREGATION`, never a class name squashed
  into one of them.

- result_folder:

  character. Path to the SEMseeker result folder.

- maxResources:

  numeric. Maximum percentage of CPU cores to use (default 90).

- parallel_strategy:

  character. Parallelisation backend; one of `"multisession"`,
  `"sequential"`, `"cluster"` (default `"multisession"`). Asking for
  `"multicore"` is accepted and converted to `"multisession"`: it means
  fork(), which is unsafe with this package's native thread pool on
  every platform that offers it, and absent on Windows. The conversion
  is logged.

- start_fresh:

  logical. If `TRUE`, delete previous inference results before running
  (default `FALSE`).

- ...:

  Additional arguments passed to
  [`core_init_env()`](https://corsaro-lab.github.io/SEMseeker/reference/core_init_env.md).

## Value

Invisibly `NULL`. Inference result CSV files are written to the
`Inference/` sub-folder of `result_folder`, one file per
marker/area/family combination defined in `inference_details`.

## Examples

``` r
result_dir <- tempdir()
if (FALSE) { # \dontrun{
association_analysis(
  inference_details = data.frame(
    independent_variable = "Sample_Group",
    family_test          = "wilcoxon",
    transformation_y     = "none",
    aggregation          = "MEAN",
    scope                = "INSTANCE"
  ),
  result_folder     = "~/semseeker_results/",
  markers           = "DELTARP",
  areas             = "GENE",
  subareas          = "WHOLE",
  multiple_test_adj = "BH"
)
} # }
```
