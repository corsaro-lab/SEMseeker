# Tests for pure (side-effect-free) helper functions
#
# Covered:
#  - core_name_cleaning()                      string sanitiser
#  - assoc_is_family_dicotomic()                family-test classifier
#  - assoc_sig_formula_vars()                   formula parser
#  - assoc_apply_stat_model_sig_formula()       formula builder
#  - assoc_exact_pvalue()                       permutation p-value calculator
#  - assoc_model_performance()                  regression metric calculator
#  - anno_sort_by_chr_and_start()              genomic-coordinate sorter
#  - assoc_compute_mean_delta_permutation_cpu() mean-difference permutation statistic
#  - assoc_compute_spearman_permutation()       Spearman permutation statistic
#  - assoc_compute_quantile_delta_permutation() quantile-difference permutation statistic
#  - util_describe_dataframe()                 per-column summary statistics
#  - assoc_polynomial_formula_build()           polynomial regression formula builder

# ---------------------------------------------------------------------------
# 1. core_name_cleaning
# ---------------------------------------------------------------------------

test_that("core_name_cleaning uppercases and replaces spaces/dots/dashes with underscores", {
  expect_equal(SEMseeker:::core_name_cleaning("hello world"), "HELLO_WORLD")
  expect_equal(SEMseeker:::core_name_cleaning("a.b.c"),       "A_B_C")
  expect_equal(SEMseeker:::core_name_cleaning("a-b-c"),       "A_B_C")
  expect_equal(SEMseeker:::core_name_cleaning("abc"),         "ABC")
})

test_that("core_name_cleaning replaces comparison operators", {
  expect_equal(SEMseeker:::core_name_cleaning("x>=1"),  "X_GTE_1")
  expect_equal(SEMseeker:::core_name_cleaning("x<=1"),  "X_LTE_1")
  expect_equal(SEMseeker:::core_name_cleaning("x==1"),  "X_EQ_1")
  expect_equal(SEMseeker:::core_name_cleaning("x!=1"),  "X_DISEQ_1")
  expect_equal(SEMseeker:::core_name_cleaning("x>1"),   "X_GT_1")
  expect_equal(SEMseeker:::core_name_cleaning("x<1"),   "X_LT_1")
})

test_that("core_name_cleaning strips leading and trailing underscores", {
  expect_false(startsWith(SEMseeker:::core_name_cleaning("_abc"), "_"))
  expect_false(endsWith(SEMseeker:::core_name_cleaning("abc_"),   "_"))
})

test_that("core_name_cleaning collapses double underscores", {
  result <- SEMseeker:::core_name_cleaning("a__b")
  expect_false(grepl("__", result))
})

test_that("core_name_cleaning handles empty strings with placeholder", {
  expect_equal(SEMseeker:::core_name_cleaning("", empty_place_holder = "EMPTY"), "EMPTY")
})

test_that("core_name_cleaning handles vectors", {
  result <- SEMseeker:::core_name_cleaning(c("hello world", "foo.bar"))
  expect_equal(result, c("HELLO_WORLD", "FOO_BAR"))
})

# ---------------------------------------------------------------------------
# 2. assoc_is_family_dicotomic
# ---------------------------------------------------------------------------

test_that("assoc_is_family_dicotomic returns FALSE for continuous families", {
  expect_false(SEMseeker:::assoc_is_family_dicotomic("gaussian"))
  expect_false(SEMseeker:::assoc_is_family_dicotomic("spearman"))
  expect_false(SEMseeker:::assoc_is_family_dicotomic("pearson"))
  expect_false(SEMseeker:::assoc_is_family_dicotomic("kendall"))
  expect_false(SEMseeker:::assoc_is_family_dicotomic("poisson"))
  expect_false(SEMseeker:::assoc_is_family_dicotomic("quantreg_0.5"))
  expect_false(SEMseeker:::assoc_is_family_dicotomic("quantreg-permutation_0.5_5_10_0.9"))
  expect_false(SEMseeker:::assoc_is_family_dicotomic("mean-permutation_100_1000_0.95"))
  expect_false(SEMseeker:::assoc_is_family_dicotomic("spearman-permutation"))
  expect_false(SEMseeker:::assoc_is_family_dicotomic("quantile-permutation"))
  expect_false(SEMseeker:::assoc_is_family_dicotomic("polynomial_4_1"))
  expect_false(SEMseeker:::assoc_is_family_dicotomic("log_1"))
})

test_that("assoc_is_family_dicotomic returns TRUE for dichotomous families", {
  expect_true(SEMseeker:::assoc_is_family_dicotomic("wilcoxon"))
  expect_true(SEMseeker:::assoc_is_family_dicotomic("t.test"))
  expect_true(SEMseeker:::assoc_is_family_dicotomic("kruskal.test"))
  expect_true(SEMseeker:::assoc_is_family_dicotomic("chisq.test"))
  expect_true(SEMseeker:::assoc_is_family_dicotomic("fisher.test"))
  expect_true(SEMseeker:::assoc_is_family_dicotomic("binomial"))
  expect_true(SEMseeker:::assoc_is_family_dicotomic("multinomial"))
  expect_true(SEMseeker:::assoc_is_family_dicotomic("jsd"))
})

# ---------------------------------------------------------------------------
# 3. assoc_sig_formula_vars
# ---------------------------------------------------------------------------

test_that("assoc_sig_formula_vars extracts dependent and independent variables", {
  f <- stats::as.formula("burden ~ group")
  res <- SEMseeker:::assoc_sig_formula_vars(f)
  expect_equal(res$dependent_variable,   "burden")
  expect_equal(res$independent_variable, "group")
  expect_length(res$covariates, 0)
})

test_that("assoc_sig_formula_vars extracts covariates correctly", {
  f <- stats::as.formula("burden ~ group + age + sex")
  res <- SEMseeker:::assoc_sig_formula_vars(f)
  expect_equal(res$dependent_variable,   "burden")
  expect_equal(res$independent_variable, "group")
  expect_equal(sort(res$covariates), c("age", "sex"))
})

# ---------------------------------------------------------------------------
# 4. assoc_apply_stat_model_sig_formula
# ---------------------------------------------------------------------------

test_that("assoc_apply_stat_model_sig_formula builds correct formula for wilcoxon", {
  f <- SEMseeker:::assoc_apply_stat_model_sig_formula("wilcoxon", "BURDEN", "GROUP", NULL)
  expect_equal(deparse(f), "BURDEN ~ GROUP")
})

test_that("assoc_apply_stat_model_sig_formula builds correct formula for gaussian with covariates", {
  f <- SEMseeker:::assoc_apply_stat_model_sig_formula("gaussian", "BURDEN", "GROUP", c("AGE", "SEX"))
  expect_true(grepl("GROUP", deparse(f)))
  expect_true(grepl("AGE",   deparse(f)))
  expect_true(grepl("SEX",   deparse(f)))
})

test_that("assoc_apply_stat_model_sig_formula inverts roles for binomial", {
  f <- SEMseeker:::assoc_apply_stat_model_sig_formula("binomial", "BURDEN", "GROUP", NULL)
  # for binomial: GROUP ~ BURDEN
  expect_equal(deparse(f), "GROUP ~ BURDEN")
})

test_that("assoc_apply_stat_model_sig_formula builds correct formula for pearson", {
  f <- SEMseeker:::assoc_apply_stat_model_sig_formula("pearson", "BURDEN", "PHENO", NULL)
  expect_equal(deparse(f), "BURDEN ~ PHENO")
})

test_that("assoc_apply_stat_model_sig_formula builds correct formula for mean-permutation", {
  f <- SEMseeker:::assoc_apply_stat_model_sig_formula("mean-permutation_100_1000_0.95", "CG_AREA", "SAMPLE_GROUP", NULL)
  expect_equal(deparse(f), "CG_AREA ~ SAMPLE_GROUP")
})

# ---------------------------------------------------------------------------
# 5. assoc_exact_pvalue
# ---------------------------------------------------------------------------

test_that("assoc_exact_pvalue returns vector of length 4", {
  set.seed(42)
  perm <- rnorm(1000, mean = 0, sd = 1)
  result <- SEMseeker:::assoc_exact_pvalue(perm, statistic_parameter = 3, conf.level = 0.95)
  expect_length(result, 4)
})

test_that("assoc_exact_pvalue: large statistic gives small p-value", {
  set.seed(42)
  perm <- rnorm(10000, mean = 0, sd = 1)
  result <- SEMseeker:::assoc_exact_pvalue(perm, statistic_parameter = 5, conf.level = 0.95)
  p_value <- result[3]
  expect_lt(p_value, 0.01)
})

test_that("assoc_exact_pvalue: statistic = median gives p-value ~ 0.5", {
  set.seed(42)
  perm <- rnorm(10000, mean = 0, sd = 1)
  med <- median(perm)
  result <- SEMseeker:::assoc_exact_pvalue(perm, statistic_parameter = med, conf.level = 0.95)
  p_value <- result[3]
  expect_lt(p_value, 0.6)
  expect_gt(p_value, 0.4)
})

test_that("assoc_exact_pvalue: CI lower < CI upper", {
  set.seed(42)
  perm <- rnorm(1000)
  result <- SEMseeker:::assoc_exact_pvalue(perm, statistic_parameter = 0, conf.level = 0.95)
  expect_lt(result[1], result[2])
})

test_that("assoc_exact_pvalue: pvalue_limit = 1 - conf.level", {
  set.seed(42)
  perm <- rnorm(1000)
  conf_level <- 0.90
  result <- SEMseeker:::assoc_exact_pvalue(perm, statistic_parameter = 0, conf.level = conf_level)
  expect_equal(unname(result[4]), 1 - conf_level)
})

# ---------------------------------------------------------------------------
# 6. assoc_model_performance
# ---------------------------------------------------------------------------

test_that("assoc_model_performance returns a data.frame with expected columns", {
  fitted   <- c(1.1, 2.0, 3.1, 4.0)
  expected <- c(1.0, 2.0, 3.0, 4.0)
  res <- SEMseeker:::assoc_model_performance(fitted, expected, c(), c())
  expect_s3_class(res, "data.frame")
  expect_true(all(c("mse", "rmse", "mae", "r_squared") %in% colnames(res)))
})

test_that("assoc_model_performance: perfect fit gives r_squared = 1", {
  x <- c(1, 2, 3, 4, 5)
  res <- SEMseeker:::assoc_model_performance(x, x, c(), c())
  expect_equal(res$mse, 0)
  expect_equal(res$mae, 0)
  expect_equal(res$rmse, 0)
  expect_equal(res$r_squared, 1)
})

test_that("assoc_model_performance: mse is mean of squared residuals", {
  fitted   <- c(2, 4)
  expected <- c(1, 3)
  res <- SEMseeker:::assoc_model_performance(fitted, expected, c(), c())
  expect_equal(res$mse, 1)   # both residuals = -1, mean((1)^2) = 1
})

test_that("assoc_model_performance: rmse = sqrt(mse)", {
  fitted   <- c(2, 2, 6)
  expected <- c(1, 3, 5)
  res <- SEMseeker:::assoc_model_performance(fitted, expected, c(), c())
  expect_equal(res$rmse, sqrt(res$mse))
})

# ---------------------------------------------------------------------------
# 7. anno_sort_by_chr_and_start
# ---------------------------------------------------------------------------

test_that("anno_sort_by_chr_and_start orders chromosomes correctly", {
  df <- data.frame(
    CHR   = c("chr2", "chr1", "chr10", "chrX"),
    START = c(100, 200, 50, 300),
    VAL   = 1:4
  )
  sorted <- SEMseeker:::anno_sort_by_chr_and_start(df)
  expect_equal(as.character(sorted$CHR[1]), "chr1")
  expect_equal(as.character(sorted$CHR[2]), "chr2")
  expect_equal(as.character(sorted$CHR[3]), "chr10")
})

test_that("anno_sort_by_chr_and_start sorts by START within same chromosome", {
  df <- data.frame(
    CHR   = c("chr1", "chr1", "chr1"),
    START = c(300, 100, 200),
    VAL   = 1:3
  )
  sorted <- SEMseeker:::anno_sort_by_chr_and_start(df)
  expect_equal(sorted$START, c(100, 200, 300))
})

test_that("anno_sort_by_chr_and_start handles empty data.frame", {
  df <- data.frame(CHR = character(0), START = integer(0))
  result <- SEMseeker:::anno_sort_by_chr_and_start(df)
  expect_equal(nrow(result), 0)
})

# ---------------------------------------------------------------------------
# 8. assoc_compute_mean_delta_permutation_cpu
# ---------------------------------------------------------------------------

test_that("assoc_compute_mean_delta_permutation_cpu returns a numeric scalar", {
  df <- data.frame(
    burden = c(1.0, 1.2, 1.1, 2.0, 2.2, 2.1),
    group  = c(0, 0, 0, 1, 1, 1)
  )
  f   <- stats::as.formula("burden ~ group")
  res <- SEMseeker:::assoc_compute_mean_delta_permutation_cpu(f, df, shuffle = FALSE)
  expect_length(res, 1)
  expect_true(is.numeric(res))
})

test_that("assoc_compute_mean_delta_permutation_cpu: group 1 > group 0 gives negative delta", {
  # group 1 has higher values → mean(g1) - mean(g0) > 0 → function returns -1 * diff
  df <- data.frame(
    burden = c(1, 1, 1, 3, 3, 3),
    group  = c(0, 0, 0, 1, 1, 1)
  )
  f   <- stats::as.formula("burden ~ group")
  res <- SEMseeker:::assoc_compute_mean_delta_permutation_cpu(f, df, shuffle = FALSE)
  expect_lt(res, 0)    # -1 * (mean(g1) - mean(g0)) = -1 * 2 = -2
})

test_that("assoc_compute_mean_delta_permutation_cpu: shuffle changes the result", {
  set.seed(123)
  df <- data.frame(
    burden = c(1, 1, 1, 10, 10, 10),
    group  = c(0, 0, 0,  1,  1,  1)
  )
  f      <- stats::as.formula("burden ~ group")
  obs    <- SEMseeker:::assoc_compute_mean_delta_permutation_cpu(f, df, shuffle = FALSE)
  perms  <- replicate(50, SEMseeker:::assoc_compute_mean_delta_permutation_cpu(f, df, shuffle = TRUE))
  # Shuffled results should not all equal the observed statistic
  expect_false(all(perms == obs))
})

# ---------------------------------------------------------------------------
# 9. assoc_compute_spearman_permutation
# ---------------------------------------------------------------------------

test_that("assoc_compute_spearman_permutation returns a numeric scalar in [-1, 1]", {
  df <- data.frame(
    burden = 1:10,
    group  = 1:10
  )
  f   <- stats::as.formula("burden ~ group")
  res <- SEMseeker:::assoc_compute_spearman_permutation(f, df, shuffle = FALSE)
  expect_length(res, 1)
  expect_true(is.numeric(res))
  expect_gte(res, -1)
  expect_lte(res,  1)
})

test_that("assoc_compute_spearman_permutation: perfectly correlated data gives rho = 1", {
  df <- data.frame(
    burden = 1:20,
    group  = 1:20
  )
  f   <- stats::as.formula("burden ~ group")
  res <- SEMseeker:::assoc_compute_spearman_permutation(f, df, shuffle = FALSE)
  expect_equal(unname(res), 1, tolerance = 1e-10)
})

test_that("assoc_compute_spearman_permutation: shuffle changes the result", {
  set.seed(42)
  df <- data.frame(
    burden = 1:30,
    group  = 1:30
  )
  f     <- stats::as.formula("burden ~ group")
  obs   <- SEMseeker:::assoc_compute_spearman_permutation(f, df, shuffle = FALSE)
  perms <- replicate(30, SEMseeker:::assoc_compute_spearman_permutation(f, df, shuffle = TRUE))
  expect_false(all(perms == obs))
})

# ---------------------------------------------------------------------------
# 10. assoc_compute_quantile_delta_permutation
# ---------------------------------------------------------------------------

test_that("assoc_compute_quantile_delta_permutation returns a numeric scalar", {
  df <- data.frame(
    burden = c(1, 2, 3, 10, 11, 12),
    group  = c(0, 0, 0,  1,  1,  1)
  )
  f   <- stats::as.formula("burden ~ group")
  res <- SEMseeker:::assoc_compute_quantile_delta_permutation(f, df, shuffle = FALSE, quantile = 0.5)
  expect_length(res, 1)
  expect_true(is.numeric(res))
})

test_that("assoc_compute_quantile_delta_permutation: Q50 of group 1 > Q50 of group 0", {
  df <- data.frame(
    burden = c(1, 2, 3, 10, 11, 12),
    group  = c(0, 0, 0,  1,  1,  1)
  )
  f   <- stats::as.formula("burden ~ group")
  res <- SEMseeker:::assoc_compute_quantile_delta_permutation(f, df, shuffle = FALSE, quantile = 0.5)
  expect_gt(res, 0)   # median(g1) - median(g0) = 11 - 2 = 9 > 0
})

test_that("assoc_compute_quantile_delta_permutation: Q10 and Q90 give different results", {
  set.seed(1)
  df <- data.frame(
    burden = c(stats::rnorm(20, 0), stats::rnorm(20, 5)),
    group  = c(rep(0, 20), rep(1, 20))
  )
  f   <- stats::as.formula("burden ~ group")
  r10 <- SEMseeker:::assoc_compute_quantile_delta_permutation(f, df, shuffle = FALSE, quantile = 0.10)
  r90 <- SEMseeker:::assoc_compute_quantile_delta_permutation(f, df, shuffle = FALSE, quantile = 0.90)
  expect_false(isTRUE(all.equal(r10, r90)))
})

test_that("assoc_compute_quantile_delta_permutation: shuffle changes the result", {
  set.seed(99)
  df <- data.frame(
    burden = c(1, 2, 3, 100, 101, 102),
    group  = c(0, 0, 0,   1,   1,   1)
  )
  f     <- stats::as.formula("burden ~ group")
  obs   <- SEMseeker:::assoc_compute_quantile_delta_permutation(f, df, shuffle = FALSE, quantile = 0.5)
  perms <- replicate(30, SEMseeker:::assoc_compute_quantile_delta_permutation(f, df, shuffle = TRUE, quantile = 0.5))
  expect_false(all(perms == obs))
})

# ---------------------------------------------------------------------------
# 11. util_describe_dataframe
# ---------------------------------------------------------------------------

test_that("util_describe_dataframe: returns one row per column", {
  df <- data.frame(x = 1:5, y = letters[1:5], z = c(1.1, NA, 2.2, NA, 3.3))
  res <- SEMseeker:::util_describe_dataframe(df)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)
})

test_that("util_describe_dataframe: Variable column matches input column names", {
  df <- data.frame(alpha = 1:3, beta = c("a", "b", "c"))
  res <- SEMseeker:::util_describe_dataframe(df)
  expect_equal(as.character(res$Variable), c("alpha", "beta"))
})

test_that("util_describe_dataframe: Missing_Values counts NAs correctly", {
  df <- data.frame(x = c(1, NA, 3, NA, 5))
  res <- SEMseeker:::util_describe_dataframe(df)
  expect_equal(res$Missing_Values, 2)
})

test_that("util_describe_dataframe: Missing_Values_Percent is 100*NA/n", {
  df <- data.frame(x = c(NA, NA, 1, 1))  # 50% missing
  res <- SEMseeker:::util_describe_dataframe(df)
  expect_equal(res$Missing_Values_Percent, 50)
})

test_that("util_describe_dataframe: Mean/Median/Min/Max are NA for non-numeric columns", {
  df <- data.frame(label = c("a", "b", "c"))
  res <- SEMseeker:::util_describe_dataframe(df)
  expect_true(is.na(res$Mean))
  expect_true(is.na(res$Median))
  expect_true(is.na(res$Min))
  expect_true(is.na(res$Max))
})

test_that("util_describe_dataframe: numeric statistics are correct", {
  df <- data.frame(v = c(1, 2, 3, 4, 5))
  res <- SEMseeker:::util_describe_dataframe(df)
  expect_equal(res$Mean,   3)
  expect_equal(res$Median, 3)
  expect_equal(res$Min,    1)
  expect_equal(res$Max,    5)
})

# ---------------------------------------------------------------------------
# 12. assoc_polynomial_formula_build
# ---------------------------------------------------------------------------

test_that("assoc_polynomial_formula_build: degree-1 without covariates gives y ~ I(x^1)", {
  f <- SEMseeker:::assoc_polynomial_formula_build("y", "x", 1, character(0))
  expect_equal(deparse(f), "y ~ I(x^1)")
})

test_that("assoc_polynomial_formula_build: degree-3 without covariates has 3 terms", {
  f    <- SEMseeker:::assoc_polynomial_formula_build("y", "x", 3, character(0))
  fstr <- paste(deparse(f), collapse = " ")
  expect_true(grepl("I(x^1)", fstr, fixed = TRUE))
  expect_true(grepl("I(x^2)", fstr, fixed = TRUE))
  expect_true(grepl("I(x^3)", fstr, fixed = TRUE))
})

# The test that used to stand here asserted the defect. It read
#   expect_true(grepl("I(x^1):age", ...))
# and it passed, because it asked the code to confirm what the code did. An
# expected value taken from the implementation agrees with the implementation
# even when the implementation is wrong, which is why nothing flagged a model
# that adjusted for nothing. These assert the model instead.

test_that("assoc_polynomial_formula_build: covariates are main effects, not interactions", {
  f    <- SEMseeker:::assoc_polynomial_formula_build("y", "x", 2, c("age", "sex"))
  fstr <- paste(deparse(f), collapse = " ")  # deparse may split long formulas

  expect_equal(fstr, "y ~ I(x^1) + I(x^2) + age + sex")

  # Said as the property rather than as the string, so a reformatting of the
  # formula cannot make this pass while the model changes: no term multiplies a
  # covariate by a power of the predictor.
  expect_false(grepl(":", fstr, fixed = TRUE))
  expect_true(all(c("age", "sex") %in% labels(stats::terms(f))))
})

test_that("assoc_polynomial_formula_build: the covariate can shift the level at x = 0", {
  # The property the missing main effect took away, stated where it can be
  # checked: with a covariate term the fit at x = 0 depends on the covariate,
  # and under the interaction-only model it could not.
  f <- SEMseeker:::assoc_polynomial_formula_build("y", "x", 2, "cov")
  mm <- stats::model.matrix(f, data.frame(y = 0, x = 0, cov = c(0, 1)))

  expect_true("cov" %in% colnames(mm))
  # at x = 0 the two rows differ only through the covariate column
  expect_equal(unname(mm[, "cov"]), c(0, 1))
  expect_equal(unname(mm[, "I(x^1)"]), c(0, 0))
})

test_that("assoc_polynomial_formula_build: an empty covariate name is dropped, not pasted", {
  # assoc_sig_formula_vars() can return "". It used to reach the formula as a
  # bare "I(x^1):" and stop the run on a parse error, so this is a path that
  # failed rather than a cosmetic case.
  f <- SEMseeker:::assoc_polynomial_formula_build("y", "x", 2, "")
  expect_equal(paste(deparse(f), collapse = " "), "y ~ I(x^1) + I(x^2)")

  f2 <- SEMseeker:::assoc_polynomial_formula_build("y", "x", 2, c("age", ""))
  expect_equal(paste(deparse(f2), collapse = " "), "y ~ I(x^1) + I(x^2) + age")
})

test_that("assoc_polynomial_formula_build: returns a formula object", {
  f <- SEMseeker:::assoc_polynomial_formula_build("outcome", "predictor", 2, character(0))
  expect_s3_class(f, "formula")
})

test_that("assoc_polynomial_formula_build: LHS is the dependent variable", {
  f <- SEMseeker:::assoc_polynomial_formula_build("response", "x", 2, character(0))
  expect_equal(as.character(f[[2]]), "response")
})
