# Enrichment analysis (pathway + phenotype)

Dispatcher per l'analisi di arricchimento di pathway e phenotype sui
geni con epimutazioni significative identificate da
\`association_analysis\`. Supporta multipli backend opt-in: WebGestalt,
STRINGdb, pathfindR, CTDbase (ctdR), e Phenolyzer.

**Phenolyzer requirement:** the `phenolyzer` and `Phenolyzer_STRINGdb` /
`Phenolyzer_WebGestalt` backends require Phenolyzer to be installed and
configured at the operating-system level (see
<https://phenolyzer.wglab.org>). Supply the path to the Phenolyzer
binary directory via `phenolyzer_folder_bin`.

## Usage

``` r
enrichment_analysis(
  inference_details,
  adjust_per_area_s,
  adjust_globally_s,
  pvalue_columns,
  adjustment_methods,
  alphas,
  study,
  significance,
  statistic_parameter,
  path_dbs,
  phenolyzer_folder_bin,
  disease,
  phenolyzer = FALSE,
  WebGestalt = FALSE,
  pathfindr = FALSE,
  STRINGdb = FALSE,
  Phenolyzer_STRINGdb = FALSE,
  Phenolyzer_WebGestalt = FALSE,
  ctdR = FALSE,
  result_folder,
  maxResources = 90,
  parallel_strategy = "multicore",
  ...
)
```

## Arguments

- inference_details:

  data.frame. Inference parameter table (must contain a `depth_analysis`
  column; rows with `depth_analysis == 3` are processed).

- adjust_per_area_s:

  logical vector. Whether to adjust p-values per area for each
  `pvalue_columns` entry.

- adjust_globally_s:

  logical vector. Whether to adjust p-values globally for each
  `pvalue_columns` entry.

- pvalue_columns:

  character vector. Column name(s) in the inference results to use as
  p-value filter.

- adjustment_methods:

  character vector. Multiple-testing correction method(s) applied (e.g.
  `"BH"`).

- alphas:

  numeric vector. Significance thresholds to iterate over.

- study:

  character. Study identifier used for labelling outputs.

- significance:

  logical. Whether to filter by significance.

- statistic_parameter:

  character. Statistic column used for ranking.

- path_dbs:

  character. Path database(s) for pathfindR.

- phenolyzer_folder_bin:

  character. Path to the Phenolyzer binary directory (required when
  `phenolyzer = TRUE` or related flags).

- disease:

  character. Disease keyword passed to Phenolyzer.

- phenolyzer:

  logical. Enable Phenolyzer backend (default `FALSE`). Requires
  Phenolyzer installed system-wide.

- WebGestalt:

  logical. Enable WebGestalt backend (default `FALSE`).

- pathfindr:

  logical. Enable pathfindR backend (default `FALSE`).

- STRINGdb:

  logical. Enable STRINGdb backend (default `FALSE`).

- Phenolyzer_STRINGdb:

  logical. Enable combined Phenolyzer+STRINGdb backend (default
  `FALSE`).

- Phenolyzer_WebGestalt:

  logical. Enable combined Phenolyzer+WebGestalt backend (default
  `FALSE`).

- ctdR:

  logical. Enable CTDbase backend via the ctdR package (default
  `FALSE`).

- result_folder:

  character. Path to the SEMseeker result folder.

- maxResources:

  numeric. Maximum percentage of CPU cores to use (default 90).

- parallel_strategy:

  character. Parallelisation backend (default `"multicore"`).

- ...:

  Additional named arguments passed to
  [`core_init_env()`](https://corsaro-lab.github.io/SEMseeker/reference/core_init_env.md).

## Value

Invisibly `NULL`. Pathway enrichment results are written to the pathway
sub-folder of `result_folder`.

## Examples

``` r
result_dir <- tempdir()
if (FALSE) { # \dontrun{
enrichment_analysis(
  inference_details  = inference_df,
  adjust_per_area_s  = TRUE,
  adjust_globally_s  = FALSE,
  pvalue_columns     = "PVALUE_ADJ_ALL_BH",
  adjustment_methods = "BH",
  alphas             = 0.05,
  result_folder      = "~/semseeker_results/",
  WebGestalt         = TRUE
)
} # }
```
