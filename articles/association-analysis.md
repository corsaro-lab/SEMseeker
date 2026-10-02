# Association analysis — statistical models, correlation and group tests

## Overview

[`association_analysis()`](https://corsaro-lab.github.io/SEMseeker/reference/association_analysis.md)
is the downstream statistical engine of semseeker. After
[`semseeker()`](https://corsaro-lab.github.io/SEMseeker/reference/semseeker.md)
has computed epimutation counts, lesion counts and delta metrics
(`DELTAS`, `DELTAR`, and their ranked variants `DELTAP`/`DELTARP` and
`DELTAQ`/`DELTARQ`), this function tests whether any of those values are
associated with phenotypic or clinical variables.

The function supports **three broad families** of models, all run per
genomic region in parallel using `future` (Bengtsson 2021):

| Family | Methods |
|----|----|
| **Correlation** | Pearson, Spearman, Kendall, permuted Spearman, Jensen–Shannon divergence |
| **Group tests** | Wilcoxon (Mann-Whitney), paired Wilcoxon, Student’s t-test, paired t-test, Kruskal-Wallis, chi-squared, Fisher’s exact, Bartlett, mean permutation |
| **Regression** | Gaussian/Poisson/Binomial GLM, quantile regression (standard and permuted), polynomial, multinomial, nonlinear (log₁₀, pow₁₀, log, exp), mediation (quantreg, ridge, linear) |

Multiple comparisons are adjusted with the Benjamini-Hochberg FDR
procedure (Benjamini and Hochberg 1995). All analyses are run across
every combination of genomic area, marker and sample group defined
during the
[`semseeker()`](https://corsaro-lab.github.io/SEMseeker/reference/semseeker.md)
call.

------------------------------------------------------------------------

## Prerequisites

You need a completed semseeker run (see the [Getting started
vignette](https://corsaro-lab.github.io/SEMseeker/articles/getting-started.md)).
The sample sheet used in the original
[`semseeker()`](https://corsaro-lab.github.io/SEMseeker/reference/semseeker.md)
call must already contain the clinical variable(s) you want to test:

``` r

# In practice: sample_sheet <- read.csv2("sample_sheet.csv").
# Here we use the synthetic sample sheet defined in the hidden setup chunk.

# Required: Sample_ID, Sample_Group
# Plus any clinical covariates you plan to test:
# age, bmi, treatment, exposure_level, ...
head(sample_sheet)
#>   Sample_ID Sample_Group age  bmi sex
#> 1      S001    Reference  34 22.1   F
#> 2      S002    Reference  41 25.4   M
#> 3      S003    Reference  29 19.8   F
#> 4      S004         Case  52 28.6   M
#> 5      S005         Case  47 26.3   F
#> 6      S006         Case  61 31.1   M
```

------------------------------------------------------------------------

## The `inference_details` data frame

Every analysis is specified through an `inference_details` data frame.
Each row is one analysis to run. The columns are:

| Column | Required | Description |
|----|----|----|
| `formula` | ✅ | R-style formula string: `"DELTARP ~ age"` or `"MUTATIONS_HYPO ~ age + sex"` |
| `family` | ✅ | Model/test type (see full table below) |
| `transformation_y` | ✅ | Pre-transformation on the dependent variable |
| `depth_analysis` | ✅ | Granularity of results (1–3, see below) |
| `filter_p_value` | ✅ | `TRUE` = keep only FDR-significant results |
| `samples_sql_condition` | ⬜ | SQL WHERE fragment to subset samples (e.g. `"Sample_Group = 'Case'"`) |
| `covariates` | ⬜ | Comma-separated string of covariate column names |
| `covariates_dummy` | ⬜ | Covariates to one-hot encode before modelling |
| `collinearity_check` | ⬜ | `TRUE` = remove collinear covariates automatically |
| `covariates_pca` | ⬜ | `TRUE` = replace covariates with their first PCs |

------------------------------------------------------------------------

## Dependent variable: what can be tested

Any semseeker output metric can be the dependent variable in `formula`:

| Metric | Description |
|----|----|
| `MUTATIONS_HYPO` | Count of hypomethylated mutations per region |
| `MUTATIONS_HYPER` | Count of hypermethylated mutations per region |
| `LESIONS_HYPO` | Count of hypomethylated lesions per region |
| `LESIONS_HYPER` | Count of hypermethylated lesions per region |
| `DELTAS` | Delta Signal — signed deviation of the value beyond the reference threshold |
| `DELTAR` | Delta Ratio — `DELTAS` divided by the reference interval width (superior − inferior) |
| `DELTAP` | `DELTAS` discretised into equal-width bins, each assigned an integer rank weight |
| `DELTARP` | `DELTAR` discretised into equal-width bins, each assigned an integer rank weight (default 4 bins) |
| `DELTAQ` | `DELTAS` discretised into quantile bins, each assigned an integer rank weight |
| `DELTARQ` | `DELTAR` discretised into quantile bins, each assigned an integer rank weight (default 4 bins) |
| `SIGNAL` | Per-probe methylation signal; its figure is the scale of the value — `BETA` for a proportion bounded in \[0,1\], `MVALUE` for the logit-transformed one |

------------------------------------------------------------------------

## All supported model types (`family`)

### Correlation

``` r

# Spearman rank correlation (recommended for non-normal SEM metrics)
inf <- data.frame(formula="DELTARP ~ age", family="spearman",
                  transformation_y="none", depth_analysis=3,
                  filter_p_value=TRUE)

# Pearson linear correlation
inf <- data.frame(formula="DELTARP ~ bmi", family="pearson",
                  transformation_y="log", depth_analysis=3,
                  filter_p_value=FALSE)

# Kendall's tau (more robust with tied ranks)
inf <- data.frame(formula="MUTATIONS_HYPO ~ age", family="kendall",
                  transformation_y="none", depth_analysis=2,
                  filter_p_value=TRUE)

# Spearman with permutation-based p-values (more accurate for small N)
# Syntax: "spearman-permutation_<n_permutations>"
inf <- data.frame(formula="DELTARP ~ age", family="spearman-permutation_1000",
                  transformation_y="none", depth_analysis=3,
                  filter_p_value=TRUE)

# Jensen-Shannon divergence (compares distributions)
inf <- data.frame(formula="DELTARP ~ Sample_Group", family="jsd",
                  transformation_y="none", depth_analysis=2,
                  filter_p_value=TRUE)

# Run the analysis (uncomment when you have a real semseeker results folder):
# association_analysis(inf, result_folder = "~/semseeker_results/")
```

Spearman (Spearman 1904) is the default recommendation for SEM metrics
because they are often skewed and non-normally distributed. Pearson
assumes bivariate normality. Kendall’s tau (Kendall 1938) is preferable
when there are many tied values.

------------------------------------------------------------------------

### Group comparison tests

``` r

# Wilcoxon rank-sum / Mann-Whitney (two independent groups)
inf <- data.frame(formula="MUTATIONS_HYPO ~ Sample_Group", family="wilcoxon",
                  transformation_y="none", depth_analysis=2,
                  filter_p_value=TRUE)

# Paired Wilcoxon (repeated measures / matched samples)
inf <- data.frame(formula="DELTARP ~ timepoint", family="wilcoxon.paired_timepoint",
                  transformation_y="none", depth_analysis=2,
                  filter_p_value=TRUE)

# Student's t-test
inf <- data.frame(formula="DELTARP ~ Sample_Group", family="t.test",
                  transformation_y="scale", depth_analysis=3,
                  filter_p_value=TRUE)

# Paired t-test
inf <- data.frame(formula="DELTARP ~ timepoint", family="t.test.paired_timepoint",
                  transformation_y="none", depth_analysis=3,
                  filter_p_value=TRUE)

# Kruskal-Wallis (> 2 groups, non-parametric ANOVA)
inf <- data.frame(formula="MUTATIONS_HYPO ~ disease_stage", family="kruskal.test",
                  transformation_y="none", depth_analysis=2,
                  filter_p_value=TRUE)

# Chi-squared test (categorical dependent variable)
inf <- data.frame(formula="MUTATIONS_HYPO ~ exposure", family="chisq.test",
                  transformation_y="none", depth_analysis=1,
                  filter_p_value=TRUE)

# Fisher's exact test (small samples, 2×2 tables)
inf <- data.frame(formula="MUTATIONS_HYPO ~ binary_exposure", family="fisher.test",
                  transformation_y="none", depth_analysis=1,
                  filter_p_value=TRUE)

# Bartlett's test for equality of variances
inf <- data.frame(formula="DELTARP ~ Sample_Group", family="bartlett.test",
                  transformation_y="none", depth_analysis=2,
                  filter_p_value=TRUE)

# Mean permutation test (non-parametric, permutation-based)
# Syntax: "mean-permutation_<n_permutations>"
inf <- data.frame(formula="DELTARP ~ Sample_Group", family="mean-permutation_1000",
                  transformation_y="none", depth_analysis=3,
                  filter_p_value=TRUE)

# Uncomment when you have a real semseeker results folder:
# association_analysis(inf, result_folder = "~/semseeker_results/")
```

The Wilcoxon rank-sum test (Wilcoxon 1945) is the recommended
non-parametric alternative to the t-test when the normality assumption
cannot be verified, which is common for mutation count data.

------------------------------------------------------------------------

### Regression models (GLM)

``` r

# Gaussian linear regression (continuous outcome)
# Formula can include multiple predictors
inf <- data.frame(formula="DELTARP ~ age + sex + bmi",
                  family="gaussian",
                  transformation_y="log",
                  depth_analysis=3,
                  filter_p_value=TRUE,
                  covariates="sex,bmi")   # declare covariates explicitly

# Poisson regression (count data — mutations, lesions)
inf <- data.frame(formula="MUTATIONS_HYPO ~ exposure_level",
                  family="poisson",
                  transformation_y="none",
                  depth_analysis=2,
                  filter_p_value=TRUE)

# Binomial / logistic regression (binary outcome)
inf <- data.frame(formula="case_control ~ DELTARP",
                  family="binomial",
                  transformation_y="none",
                  depth_analysis=3,
                  filter_p_value=TRUE)

# Multinomial regression (> 2 categories)
inf <- data.frame(formula="disease_stage ~ DELTARP",
                  family="multinomial",
                  transformation_y="none",
                  depth_analysis=3,
                  filter_p_value=TRUE)

# Uncomment when you have a real semseeker results folder:
# association_analysis(inf, result_folder = "~/semseeker_results/")
```

GLMs (McCullagh and Nelder 1989) are appropriate when the relationship
between the SEM metric and the clinical variable is expected to be
linear (Gaussian), when the outcome is a count (Poisson), or when
predicting a binary endpoint from SEM scores (Binomial).

------------------------------------------------------------------------

### Quantile regression

``` r

# Standard quantile regression at quantile τ
# Syntax: "quantreg_<tau>"
inf <- data.frame(formula="DELTARP ~ age",
                  family="quantreg_0.5",       # median regression
                  transformation_y="none",
                  depth_analysis=3,
                  filter_p_value=TRUE)

# Quantile regression with permutation-based inference
# Syntax: "quantreg-permutation_<tau>_<n_permutations>"
inf <- data.frame(formula="DELTARP ~ age",
                  family="quantreg-permutation_0.5_1000",
                  transformation_y="none",
                  depth_analysis=3,
                  filter_p_value=TRUE)

# Quantile permutation model
# Syntax: "quantile-permutation_<n_permutations>"
inf <- data.frame(formula="DELTARP ~ Sample_Group",
                  family="quantile-permutation_1000",
                  transformation_y="none",
                  depth_analysis=2,
                  filter_p_value=TRUE)

# Uncomment when you have a real semseeker results folder:
# association_analysis(inf, result_folder = "~/semseeker_results/")
```

Quantile regression (Koenker 2005; Geraci and Bottai 2014) is
particularly useful for SEM metrics: unlike mean regression it is robust
to outliers and captures the full conditional distribution, which
matters when only high-quantile epigenetic variability is biologically
relevant.

------------------------------------------------------------------------

### Nonlinear regression (NLS)

``` r

# Log₁₀ dose-response model
# Syntax: "log10_<scale_factor>"
inf <- data.frame(formula="DELTARP ~ dose",
                  family="log10_1",
                  transformation_y="none",
                  depth_analysis=3,
                  filter_p_value=TRUE)

# Power of 10 model
inf <- data.frame(formula="DELTARP ~ dose",
                  family="pow10_1",
                  transformation_y="none",
                  depth_analysis=3,
                  filter_p_value=TRUE)

# Natural log model
inf <- data.frame(formula="DELTARP ~ concentration",
                  family="log_1",
                  transformation_y="none",
                  depth_analysis=3,
                  filter_p_value=TRUE)

# Uncomment when you have a real semseeker results folder:
# association_analysis(inf, result_folder = "~/semseeker_results/")
```

------------------------------------------------------------------------

### Polynomial regression

``` r

# Syntax: "polynomial_<degree>"
inf <- data.frame(formula="DELTARP ~ age",
                  family="polynomial_2",       # quadratic
                  transformation_y="none",
                  depth_analysis=3,
                  filter_p_value=TRUE)

# Uncomment when you have a real semseeker results folder:
# association_analysis(inf, result_folder = "~/semseeker_results/")
```

------------------------------------------------------------------------

### Mediation analysis

Mediation models test whether a third variable (mediator) explains the
relationship between the SEM metric and an outcome:

``` r

# Mediation via quantile regression
inf <- data.frame(formula="DELTARP ~ exposure + mediator",
                  family="mediation-quantreg_0.5",
                  transformation_y="none",
                  depth_analysis=3,
                  filter_p_value=TRUE)

# Mediation via ridge regression (regularised, handles multicollinearity)
inf <- data.frame(formula="DELTARP ~ exposure + mediator",
                  family="mediation-ridge",
                  transformation_y="none",
                  depth_analysis=3,
                  filter_p_value=TRUE)

# Mediation via linear model
inf <- data.frame(formula="DELTARP ~ exposure + mediator",
                  family="mediation-linear",
                  transformation_y="none",
                  depth_analysis=3,
                  filter_p_value=TRUE)

# Uncomment when you have a real semseeker results folder:
# association_analysis(inf, result_folder = "~/semseeker_results/")
```

------------------------------------------------------------------------

## Covariates and confounders

Include covariates directly in the formula and declare them in the
`covariates` column so semseeker correctly separates the independent
variable of interest from confounders:

``` r

inf <- data.frame(
  formula          = "DELTARP ~ tcdd_exposure + age + sex + bmi",
  family           = "gaussian",
  transformation_y = "log",
  depth_analysis   = 3,
  filter_p_value   = TRUE,
  covariates       = "age,sex,bmi",        # confounders to adjust for
  collinearity_check = TRUE,               # auto-remove collinear covariates
  covariates_pca   = FALSE                 # TRUE = replace with PCs
)

# Uncomment when you have a real semseeker results folder:
# association_analysis(inf, result_folder = "~/semseeker_results/")
```

When `collinearity_check = TRUE`, semseeker automatically detects and
removes collinear covariates using VIF before fitting the model. When
`covariates_pca = TRUE`, covariates are replaced by their principal
components before fitting.

------------------------------------------------------------------------

## Filtering by sample subset (`samples_sql_condition`)

Run the analysis on a subset of samples using a SQL WHERE clause:

``` r

# Only Case samples
inf <- data.frame(
  formula               = "DELTARP ~ age",
  family                = "spearman",
  transformation_y      = "none",
  depth_analysis        = 3,
  filter_p_value        = TRUE,
  samples_sql_condition = "Sample_Group = 'Case'"
)

# Only samples with age > 40
inf <- data.frame(
  formula               = "DELTARP ~ bmi",
  family                = "pearson",
  transformation_y      = "none",
  depth_analysis        = 3,
  filter_p_value        = TRUE,
  samples_sql_condition = "age > 40"
)

# Uncomment when you have a real semseeker results folder:
# association_analysis(inf, result_folder = "~/semseeker_results/")
```

------------------------------------------------------------------------

## Running multiple analyses at once

Combine multiple rows in `inference_details` to run all analyses in a
single call:

``` r

inference_details <- rbind(
  data.frame(formula="DELTARP ~ age",             family="spearman",
             transformation_y="none",             depth_analysis=3,
             filter_p_value=TRUE,                 covariates=NA_character_,
             stringsAsFactors=FALSE),
  data.frame(formula="MUTATIONS_HYPO ~ Sample_Group", family="wilcoxon",
             transformation_y="none",             depth_analysis=2,
             filter_p_value=TRUE,                 covariates=NA_character_,
             stringsAsFactors=FALSE),
  data.frame(formula="DELTARP ~ bmi",             family="quantreg_0.5",
             transformation_y="log",              depth_analysis=3,
             filter_p_value=TRUE,                 covariates=NA_character_,
             stringsAsFactors=FALSE),
  data.frame(formula="MUTATIONS_HYPO ~ exposure + age + sex", family="poisson",
             transformation_y="none",             depth_analysis=3,
             filter_p_value=TRUE,                 covariates="age,sex",
             stringsAsFactors=FALSE)
)

# Uncomment when you have a real semseeker results folder:
# association_analysis(
#   inference_details = inference_details,
#   result_folder     = "~/semseeker_results/",
#   parallel_strategy = "multicore",
#   maxResources      = 90
# )
```

------------------------------------------------------------------------

## Transformation options (`transformation_y`)

Pre-transformations are applied to the SEM metric before testing:

| Value | Transformation | When to use |
|----|----|----|
| `"none"` | No transformation | Metric already approximately normal |
| `"log"` | Natural log: ln(y + ε) | Right-skewed counts |
| `"log2"` | Log₂ | Fold-change interpretation |
| `"log10"` | Log₁₀ | Wide dynamic range |
| `"scale"` | Z-score: (y − μ) / σ | Comparing coefficients across regions |
| `"exp"` | Exponential: e^y | Log-scale input |
| `"quantile_3"` | Quantile normalisation, 3 bins | Non-parametric normalisation |

------------------------------------------------------------------------

## Depth levels (`depth_analysis`)

| Value | Breakdown | Use case |
|----|----|----|
| `1` | Sample-level aggregate only | Quick per-sample overview |
| `2` | Per marker type (MUTATIONS, LESIONS, DELTARP, …) | Identify which metric drives the association |
| `3` | Per genomic area (GENE body, TSS200, TSS1500, island, …) | Fine-grained genomic localisation |

Higher depth levels include all lower levels.

------------------------------------------------------------------------

## Output files

Results are written to `result_folder/Inference/`. File names encode the
full analysis specification:

    Inference/
    ├── DELTARP_GENE_spearman_age_depth3.csv
    ├── MUTATIONS_HYPO_GENE_wilcoxon_Sample_Group_depth2.csv
    ├── DELTARP_GENE_gaussian_tcdd_exposure_depth3.csv
    └── ...

Each CSV contains:

| Column | Description |
|----|----|
| `CHR`, `START`, `END` | Genomic coordinates |
| `AREA`, `SUBAREA` | Aggregation level and sub-region |
| `MARKER` | SEM metric tested |
| `SCORE` | Test statistic (ρ, W, t, β, …) |
| `P_Value` | Nominal p-value |
| `PVALUE_ADJ_KEY_<m>` | Adjusted within the identity key — the instances tested for one `(MARKER, FIGURE, SCOPE, AREA, SUBAREA, AGGREGATION)` |
| `PVALUE_ADJ_SCOPE_<m>` | Adjusted across every row of the same `SCOPE` — all its region classes, figures and aggregations |
| `PVALUE_ADJ_ALL_<m>` | Adjusted across the whole file — everything tested on this marker under this model |
| `SIGNIFICATIVE_ADJ` | `TRUE` if the `ALL` level is below α (default 0.05) |
| `FAMILY_TEST` | Model used |
| `SAMPLE_GROUP` | Sample subset tested |

The three adjusted columns differ **only** by the family the correction
was controlled over, and `<m>` is the method the run declared
(`multiple_test_adj`). Which one to report is an analytical choice that
belongs in your Methods:

- at `SCOPE = INSTANCE` the `KEY` level is the usual one — its members
  are the genes (or islands, or cytobands) tested for one measurement;
- at `SCOPE = SAMPLE` the key holds a **single row**, so
  `PVALUE_ADJ_KEY_<m>` equals the raw p-value by arithmetic. The `SCOPE`
  level is the one with members to count there;
- the three are **not** a severity ladder. Benjamini-Hochberg is
  adaptive, so a wider family is usually — but not necessarily — more
  conservative.

------------------------------------------------------------------------

## Reading and exploring results

``` r

# Read a Spearman result file
results <- read.csv(
  list.files("~/semseeker_results/Inference",
             pattern = "spearman.*age.*depth3\\.csv$",
             full.names = TRUE)[1]
)

# Top significant genes by adjusted p-value
sig <- subset(results, SIGNIFICATIVE_ADJ == TRUE)
head(sig[order(sig$Q), c("CHR","START","END","AREA","SCORE","P_Value","Q")], 20)

# Volcano-style: effect size vs -log10(Q)
plot(-log10(results$Q), results$SCORE,
     xlab = "-log10(Q)", ylab = "Spearman rho",
     col  = ifelse(results$SIGNIFICATIVE_ADJ, "red", "grey"),
     pch  = 16, main = "Spearman rho vs FDR")
abline(h = 0, lty = 2)
```

------------------------------------------------------------------------

## Cross-study and cross-subsample analyses

After running
[`association_analysis()`](https://corsaro-lab.github.io/SEMseeker/reference/association_analysis.md)
on two or more independent cohorts stored in separate `result_folder`s,
semseeker provides functions to compare and meta-analyse results across
studies:

``` r

# Overlap significant hits across two studies
meta_association_overlaps_studies(
  studies       = list("~/study_A/", "~/study_B/"),
  result_folder = "~/meta_results/",
  pvalue_column = "PVALUE_ADJ_ALL_BH"
)

# Meta-analysis (combines p-values / effect sizes)
meta_association_across_studies(
  studies       = list("~/study_A/", "~/study_B/"),
  result_folder = "~/meta_results/"
)

# Subsample stability analysis (bootstrap cohort sub-sampling)
meta_association_overlaps_subsamples(
  result_folder = "~/semseeker_results/",
  n_subsamples  = 100,
  fraction      = 0.8
)
```

------------------------------------------------------------------------

## Complementary: Bayesian posterior probability

While
[`association_analysis()`](https://corsaro-lab.github.io/SEMseeker/reference/association_analysis.md)
reports frequentist p-values / q-values for each region,
[`assoc_bayes_analysis()`](https://corsaro-lab.github.io/SEMseeker/reference/assoc_bayes_analysis.md)
provides the *probabilistic complement*: per-region posterior
probabilities computed via empirical Bayes.

For each `(MARKER, AREA, SUBAREA, AREA_OF_TEST)` it estimates:

- `P(case | epimutated)` — how likely a sample is a case given the SEM
  hit
- `P(control | not-epimutated)` — how likely a sample is a control given
  absence of the SEM hit

A region is reported as a candidate biomarker when both posteriors
exceed the configured thresholds (`bayes_case_threshold` and
`bayes_control_threshold`). Output is written as CSV files under
`<result_folder>/Euristic/`, sorted by posterior probability.

``` r

assoc_bayes_analysis(
  result_folder           = "~/semseeker_results/",
  independent_variable    = "Sample_Group",
  bayes_case_threshold    = 0.9,   # min P(case | epimutated)
  bayes_control_threshold = 0.1,   # max P(control | epimutated)
  parallel_strategy       = "multicore",
  maxResources            = 90
)
```

This is most useful when the frequentist q-value ranking is similar
across many regions and you want a probability-based criterion to triage
candidates. The two analyses are independent — running both gives you a
frequentist and a Bayesian view of the same data.

------------------------------------------------------------------------

## Diagnostic performance: sensitivity and specificity

For binary outcomes (case vs control),
[`diagnostic_performance()`](https://corsaro-lab.github.io/SEMseeker/reference/diagnostic_performance.md)
computes per-region sensitivity, specificity, Jensen-Shannon distance,
and a combined diagnostic score — treating each SEM mutation as a
candidate diagnostic marker.

``` r

diagnostic_performance(
  combinations         = list(c("CASE_A", "CASE_B")),
  result_folder        = "~/semseeker_results/",
  independent_variable = "Sample_Group",
  parallel_strategy    = "multicore"
)
```

Output (per marker) is sorted descending by `SCORE`; the top rows are
the regions with the most promising diagnostic profile in the studied
cohort. CSV files are written to `<result_folder>/Euristic/`.

------------------------------------------------------------------------

## Cross-study replication

[`meta_association_replication()`](https://corsaro-lab.github.io/SEMseeker/reference/meta_association_replication.md)
reruns the association on a target study restricted to the regions that
were already significant in a reference study, then merges the two
result sets. It’s the standard “replication / validation” workflow when
you have a discovery cohort and a validation cohort.

``` r

meta_association_replication(
  inference_details_origin = inference_discovery,
  inference_details        = inference_validation,
  result_folder            = "~/semseeker_replication/",
  parallel_strategy        = "multicore"
)
```

------------------------------------------------------------------------

## Downstream: pathway and enrichment analysis

Significant genomic regions from
[`association_analysis()`](https://corsaro-lab.github.io/SEMseeker/reference/association_analysis.md)
can be fed directly into pathway enrichment tools:

``` r

# WebGestaltR — ORA or GSEA
enrich_WebGestalt(
  study            = "my_study",
  types            = c("BP", "MF"),
  enrich_methods   = c("ORA"),
  inference_detail = inf,
  significance     = TRUE
)

# STRINGdb — protein network analysis
enrich_STRINGdb(
  study            = "my_study",
  inference_detail = inf,
  significance     = TRUE
)

# pathfindR — active subnetwork enrichment
enrich_pathfindR(
  study            = "my_study",
  inference_detail = inf
)
```

See the [Pathway and enrichment analysis
vignette](https://corsaro-lab.github.io/SEMseeker/articles/pathway-analysis.md)
for full details.

------------------------------------------------------------------------

## R packages used

Statistical models are implemented using R (R Core Team 2024) base
`stats` functions, with parallelisation via `future` (Bengtsson 2021)
and `doFuture`. Quantile regression uses `lqmm` (Geraci and Bottai 2014)
and the `quantreg` package (Koenker 2005).

## Session info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: aarch64-apple-darwin23
#> Running under: macOS Tahoe 26.6.2
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: UTC
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.61         cachem_1.1.0      knitr_1.52        htmltools_0.5.9  
#>  [9] rmarkdown_2.32    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.1     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.6.1    tools_4.6.1       ragg_1.5.2        bslib_0.12.0     
#> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
#> [25] rlang_1.3.0       fs_2.1.0          htmlwidgets_1.6.4
```

------------------------------------------------------------------------

## References

Bengtsson, Henrik. 2021. “A Unifying Framework for Parallel and
Distributed Processing in R Using Futures.” *The R Journal* 13 (2):
208–27. <https://doi.org/10.32614/RJ-2021-048>.

Benjamini, Yoav, and Yosef Hochberg. 1995. “Controlling the False
Discovery Rate: A Practical and Powerful Approach to Multiple Testing.”
*Journal of the Royal Statistical Society: Series B* 57 (1): 289–300.
<https://doi.org/10.1111/j.2517-6161.1995.tb02031.x>.

Geraci, Marco, and Matteo Bottai. 2014. “Linear Quantile Mixed Models.”
*Statistics and Computing* 24 (3): 461–79.
<https://doi.org/10.1007/s11222-013-9381-9>.

Kendall, Maurice G. 1938. “A New Measure of Rank Correlation.”
*Biometrika* 30 (1/2): 81–93. <https://doi.org/10.2307/2332226>.

Koenker, Roger. 2005. *Quantile Regression*. Cambridge University Press.
<https://doi.org/10.1017/CBO9780511754098>.

McCullagh, Peter, and John A. Nelder. 1989. *Generalized Linear Models*.
2nd ed. Chapman; Hall. <https://doi.org/10.1007/978-1-4899-3242-6>.

R Core Team. 2024. *R: A Language and Environment for Statistical
Computing*. R Foundation for Statistical Computing.
<https://www.R-project.org/>.

Spearman, Charles. 1904. “The Proof and Measurement of Association
Between Two Things.” *The American Journal of Psychology* 15 (1):
72–101. <https://doi.org/10.2307/1412159>.

Wilcoxon, Frank. 1945. “Individual Comparisons by Ranking Methods.”
*Biometrics Bulletin* 1 (6): 80–83. <https://doi.org/10.2307/3001968>.
