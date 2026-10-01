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
  parallel_strategy = "multicore",
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

  marker

  :   SEM metric column prefix (e.g. `"DELTARP"`, `"MUTATIONS"`).

  depth_analysis

  :   Integer depth: `1` = sample level, `2` = type level (gene, DMR,
      CpG island), `3` = genomic area (TSS1550, WHOLE, TSS200, …).

- result_folder:

  character. Path to the SEMseeker result folder.

- maxResources:

  numeric. Maximum percentage of CPU cores to use (default 90).

- parallel_strategy:

  character. Parallelisation backend; possible values: `"none"`,
  `"multisession"`, `"sequential"`, `"multicore"`, `"cluster"` (default
  `"multicore"`).

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
    marker               = "DELTARP",
    areas                = "GENE"
  ),
  result_folder     = "~/semseeker_results/",
  multiple_test_adj = "BH"
)
} # }
```
