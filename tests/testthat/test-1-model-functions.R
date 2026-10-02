## test-1-model-functions.R
## Session-based tests for model functions that require an active SEMseeker session.
##
## Covered:
##   assoc_quantreg_model          — quantile regression (lqmm), tau in result
##   assoc_mean_permutation        — CPU permutation test, p-value in result
##   assoc_test_model_paired       — wilcoxon.paired branch
##   assoc_covariates_model        — no-op pass-through (no scaling, no PCA, no dummies)
##   assoc_model_polynomial — polynomial lm, degree in result  [requires caret]
##
## Each test takes its own session folder from sem_test_folder(), to avoid
## collisions with other test files.

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

.make_key2 <- function() {
  list(AREA = "GENE", SUBAREA = "TSS200", MARKER = "MUTATIONS", FIGURE = "K850")
}

.two_group_cont <- function(n = 15L, seed = 2L) {
  set.seed(seed)
  data.frame(
    BURDEN = c(stats::rnorm(n, mean = 1), stats::rnorm(n, mean = 3)),
    GROUP  = factor(c(rep("ctrl", n), rep("case", n)))
  )
}

# ---------------------------------------------------------------------------
# assoc_quantreg_model
# ---------------------------------------------------------------------------

test_that("assoc_quantreg_model returns a data.frame with tau column", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  set.seed(10)
  df  <- data.frame(x = stats::rnorm(40), y = stats::rnorm(40))
  f   <- stats::as.formula("y ~ x")
  key <- .make_key2()

  res <- SEMseeker:::assoc_quantreg_model(
    family_test           = "quantreg_0.5",
    sig.formula           = f,
    tempDataFrame         = df,
    independent_variable  = "x",
    transformation_y      = "",
    plot                  = FALSE,
    samples_sql_condition = "",
    key                   = key
  )

  expect_s3_class(res, "data.frame")
  expect_true("tau" %in% colnames(res))
  expect_equal(res$tau, 0.5)
})

test_that("assoc_quantreg_model: tau is preserved correctly for 0.25 quantile", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  set.seed(11)
  df  <- data.frame(x = stats::rnorm(40), y = stats::rnorm(40))
  f   <- stats::as.formula("y ~ x")
  key <- .make_key2()

  res <- SEMseeker:::assoc_quantreg_model(
    family_test           = "quantreg_0.25",
    sig.formula           = f,
    tempDataFrame         = df,
    independent_variable  = "x",
    transformation_y      = "",
    plot                  = FALSE,
    samples_sql_condition = "",
    key                   = key
  )

  expect_s3_class(res, "data.frame")
  expect_equal(res$tau, 0.25)
})

# ---------------------------------------------------------------------------
# assoc_mean_permutation
# ---------------------------------------------------------------------------

test_that("assoc_mean_permutation returns a data.frame with pvalue", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  set.seed(5)
  df  <- .two_group_cont(n = 15L)
  f   <- stats::as.formula("BURDEN ~ GROUP")
  key <- .make_key2()

  res <- SEMseeker:::assoc_mean_permutation(
    family_test           = "mean-permutation_20_20_0.95",
    sig.formula           = f,
    tempDataFrame         = df,
    independent_variable  = "GROUP",
    plot                  = FALSE,
    samples_sql_condition = "",
    key                   = key
  )

  expect_s3_class(res, "data.frame")
  expect_true("pvalue" %in% colnames(res))
  expect_true(is.numeric(res$pvalue))
})

test_that("assoc_mean_permutation: well-separated groups give small p-value", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # Group 0 = 0, group 1 = 10: delta should be extreme and unique across permutations
  df <- data.frame(
    BURDEN = c(rep(0, 20), rep(10, 20)),
    GROUP  = factor(c(rep("ctrl", 20), rep("case", 20)))
  )
  f   <- stats::as.formula("BURDEN ~ GROUP")
  key <- .make_key2()

  res <- SEMseeker:::assoc_mean_permutation(
    family_test           = "mean-permutation_50_50_0.95",
    sig.formula           = f,
    tempDataFrame         = df,
    independent_variable  = "GROUP",
    plot                  = FALSE,
    samples_sql_condition = "",
    key                   = key
  )

  expect_lt(res$pvalue, 0.05)
})

# ---------------------------------------------------------------------------
# assoc_test_model_paired  (wilcoxon.paired branch)
# ---------------------------------------------------------------------------

test_that("assoc_test_model_paired wilcoxon.paired returns data.frame with pvalue", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # 5 patients, pre/post measurements
  df <- data.frame(
    BURDEN     = c(1, 2, 3, 4, 5,   3, 4, 5, 6, 7),
    GROUP      = factor(c(rep("pre", 5), rep("post", 5))),
    PATIENT_ID = c(1, 2, 3, 4, 5,   1, 2, 3, 4, 5)
  )
  key <- .make_key2()
  f   <- stats::as.formula("BURDEN ~ GROUP")

  res <- SEMseeker:::assoc_test_model_paired(
    family_test           = "wilcoxon.paired@PATIENT_ID",
    tempDataFrame         = df,
    sig.formula           = f,
    burdenValue           = "BURDEN",
    independent_variable  = "GROUP",
    transformation_y      = "",
    plot                  = FALSE,
    samples_sql_condition = "",
    key                   = key
  )

  expect_s3_class(res, "data.frame")
  expect_true("pvalue" %in% colnames(res))
  expect_true(is.numeric(res$pvalue))
})

test_that("assoc_test_model_paired wilcoxon.paired: pre/post shift gives small p-value", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # 20 patients, post = pre + 5 (strong systematic shift)
  set.seed(42)
  n <- 20
  pre  <- stats::rnorm(n, mean = 0)
  post <- pre + 5
  df <- data.frame(
    BURDEN     = c(pre, post),
    GROUP      = factor(c(rep("pre", n), rep("post", n))),
    PATIENT_ID = c(seq_len(n), seq_len(n))
  )
  key <- .make_key2()
  f   <- stats::as.formula("BURDEN ~ GROUP")

  res <- SEMseeker:::assoc_test_model_paired(
    family_test           = "wilcoxon.paired@PATIENT_ID",
    tempDataFrame         = df,
    sig.formula           = f,
    burdenValue           = "BURDEN",
    independent_variable  = "GROUP",
    transformation_y      = "",
    plot                  = FALSE,
    samples_sql_condition = "",
    key                   = key
  )

  expect_lt(res$pvalue, 0.01)
})

test_that("assoc_test_model_paired: >2 group levels returns NA pvalue early", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df <- data.frame(
    BURDEN     = 1:15,
    GROUP      = factor(c(rep("a", 5), rep("b", 5), rep("c", 5))),
    PATIENT_ID = 1:15
  )
  key <- .make_key2()
  f   <- stats::as.formula("BURDEN ~ GROUP")

  res <- SEMseeker:::assoc_test_model_paired(
    family_test           = "wilcoxon.paired@PATIENT_ID",
    tempDataFrame         = df,
    sig.formula           = f,
    burdenValue           = "BURDEN",
    independent_variable  = "GROUP",
    transformation_y      = "",
    plot                  = FALSE,
    samples_sql_condition = "",
    key                   = key
  )

  expect_true(is.na(res$pvalue))
})

# ---------------------------------------------------------------------------
# assoc_covariates_model  (no-op: no scaling, no PCA, no collinearity, no dummies)
# ---------------------------------------------------------------------------

test_that("assoc_covariates_model: no-op returns list with covariates and study_summary", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  set.seed(7)
  study_summary <- data.frame(
    Sample_ID = paste0("S", seq_len(20)),
    GROUP     = c(rep("ctrl", 10), rep("case", 10)),
    BURDEN    = c(stats::rnorm(10, 1), stats::rnorm(10, 3)),
    Cov1      = stats::rnorm(20),
    stringsAsFactors = FALSE
  )

  inference_detail <- list(
    collinearity_check    = FALSE,
    covariates_dummy      = "",
    covariates_pca        = FALSE,
    covariates            = "Cov1",
    independent_variable  = "GROUP",
    transformation_x      = "none",
    family_test           = "wilcoxon",
    transformation_y      = "none",
    samples_sql_condition = NULL
  )

  result <- SEMseeker:::assoc_covariates_model(inference_detail, study_summary)

  expect_type(result, "list")
  expect_true("covariates" %in% names(result))
  expect_true("study_summary" %in% names(result))
  expect_equal(result$covariates, "Cov1")
})

# ---------------------------------------------------------------------------
# assoc_model_polynomial  [requires caret]
# ---------------------------------------------------------------------------

test_that("assoc_model_polynomial: degree-2 no-covariate returns PL_DEGREE", {
  skip_if_not_installed("caret")

  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  set.seed(3)
  df  <- data.frame(x = seq_len(40), y = seq_len(40) + stats::rnorm(40))
  f   <- stats::as.formula("y ~ x")
  key <- .make_key2()

  res <- SEMseeker:::assoc_model_polynomial(
    family_test           = "polynomial_2_0.8",
    tempDataFrame         = df,
    sig.formula           = f,
    transformation_y      = "",
    plot                  = FALSE,
    samples_sql_condition = "",
    key                   = key
  )

  expect_s3_class(res, "data.frame")
  expect_equal(res$PL_DEGREE, 2)
  expect_equal(res$PL_PERC,   0.8)
})

test_that("assoc_model_polynomial: with covariate exercises assoc_polynomial_formula_build", {
  skip_if_not_installed("caret")

  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  set.seed(9)
  n  <- 40
  df <- data.frame(
    y   = seq_len(n) + stats::rnorm(n),
    x   = seq_len(n),
    cov = factor(c(rep("A", n / 2), rep("B", n / 2)))
  )
  f   <- stats::as.formula("y ~ x + cov")
  key <- .make_key2()

  res <- SEMseeker:::assoc_model_polynomial(
    family_test           = "polynomial_2_0.8",
    tempDataFrame         = df,
    sig.formula           = f,
    transformation_y      = "",
    plot                  = FALSE,
    samples_sql_condition = "",
    key                   = key
  )

  expect_s3_class(res, "data.frame")
  expect_equal(res$PL_DEGREE, 2)
})
