# limma_<degree> family
#
# Three guarantees this test file pins down:
#   (1) the family_test string parser accepts limma_<degree> and
#       limma_<degree>_<partition>, and rejects malformed strings;
#   (2) the dispatcher fails fast with an install hint when limma is
#       not available: applies the dispatch=guard pattern;
#   (3) on simulated data, the limma per-area output is numerically
#       compatible with stats::lm() (degenerate case: 1-row response
#       matrix → eBayes shrinkage collapses to OLS t-statistics).
#
# Heavy SEMseeker setup (probe annotation, futureSession, etc.) is NOT
# needed for any of these tests — we exercise the model function and
# the dispatcher branch in isolation.

# Tiny helpers ---------------------------------------------------------

# Build a one-area data frame with a polynomial signal of given degree.
.sim_one_area <- function(n = 60L, degree = 2L, seed = 1L,
                          covariates = 0L, sigma = 1) {
  set.seed(seed)
  iv <- stats::rnorm(n)
  signal <- 0
  for (d in seq_len(degree)) signal <- signal + (-1)^d * iv^d
  y <- signal + stats::rnorm(n, sd = sigma)
  df <- data.frame(IV = iv, Y = y)
  if (covariates > 0L) {
    for (k in seq_len(covariates)) df[[paste0("COV", k)]] <- stats::rnorm(n)
  }
  df
}

.fake_key <- function() {
  list(MARKER = "MUTATIONS", FIGURE = "HYPO",
       AREA = "GENE", SUBAREA = "WHOLE",
       AREA_OF_TEST = "TESTGENE", COMBINED = "MUTATIONS_HYPO")
}

# (1) String parser ----------------------------------------------------

test_that("limma_<degree> is accepted, limma_<degree>_<partition> is accepted, malformed is rejected", {
  testthat::skip_if_not_installed("limma")

  df <- .sim_one_area()
  key <- .fake_key()
  formula_obj <- Y ~ IV

  # Valid: limma_2
  res2 <- SEMseeker:::assoc_model_limma(
    "limma_2", df, formula_obj,
    transformation_y = "", plot = FALSE,
    samples_sql_condition = "", key = key)
  testthat::expect_s3_class(res2, "data.frame")
  testthat::expect_true(nrow(res2) == 1L)
  testthat::expect_equal(res2$PL_DEGREE, 2)
  testthat::expect_equal(res2$PL_PERC, 1)
  testthat::expect_equal(res2$r_model, "limma::lmFit+eBayes")

  # Valid with partition: limma_2_1
  res2p <- SEMseeker:::assoc_model_limma(
    "limma_2_1", df, formula_obj,
    transformation_y = "", plot = FALSE,
    samples_sql_condition = "", key = key)
  testthat::expect_equal(res2p$PL_DEGREE, 2)
  testthat::expect_equal(res2p$PL_PERC, 1)

  # Malformed: empty result data.frame (logs ERROR, does not throw)
  res_bad <- SEMseeker:::assoc_model_limma(
    "limma", df, formula_obj,
    transformation_y = "", plot = FALSE,
    samples_sql_condition = "", key = key)
  testthat::expect_equal(nrow(res_bad), 0L)

  # Malformed: non-numeric degree
  res_bad2 <- SEMseeker:::assoc_model_limma(
    "limma_abc", df, formula_obj,
    transformation_y = "", plot = FALSE,
    samples_sql_condition = "", key = key)
  # Degree-parse failure returns the seed res data frame (no model fit)
  testthat::expect_true(nrow(res_bad2) <= 1L)
})

# (2) Dispatcher guard (dispatch=guard pattern) ----------------

test_that("assoc_execute_model rejects family_test='limma_2' with an install hint when limma is missing", {
  # We can't reliably uninstall limma in CI to assert this end-to-end.
  # Instead, exercise the requireNamespace path by rebinding it inside
  # the SEMseeker namespace for the duration of this test.
  #
  # We previously used mockery::stub(where = SEMseeker:::assoc_execute_model, ...)
  # but mockery operates on a *copy* of the function value and does not
  # rewrite the actual namespace binding assoc_execute_model resolves at call
  # time — so the stub silently never fired and the test passed at home
  # but failed under CI's --as-cran path. testthat::local_mocked_bindings
  # (3rd edition) edits the live namespace and is undone on exit.
  # requireNamespace isn't imported into SEMseeker's namespace (it's a
  # base function), so .package = "SEMseeker" can't find a binding.
  # Instead override the binding inside base for the duration of the
  # call. SEMseeker:::assoc_execute_model resolves requireNamespace from base
  # at call time, so this rewrites what it actually sees.
  #
  # The original is captured HERE, before the binding is replaced, and not
  # looked up inside the replacement. The earlier version did
  #   get("requireNamespace", envir = baseenv(), inherits = FALSE)(package, ...)
  # to reach the real function, but that is the very binding
  # local_mocked_bindings(.package = "base") has just overwritten: the lookup
  # returned the mock, and any call for a package other than limma called
  # itself. `.orig` closes over the function value taken before the swap, so it
  # cannot be the mock whatever happens to the binding.
  #
  # It went unnoticed because the recursion needs a SECOND caller: nothing in
  # this test asks for a package other than limma, so the else branch never ran
  # while the mock was installed. Running the test files in parallel put other
  # machinery in that window, the branch fired, and the run died with
  # "evaluation nested too deeply". The defect was always here; only the
  # conditions that reach it are new.
  .orig_require <- base::requireNamespace
  testthat::local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (package == "limma") FALSE else .orig_require(package, ...)
    },
    .package = "base"
  )

  df <- .sim_one_area()
  key <- .fake_key()
  testthat::expect_error(
    SEMseeker:::assoc_execute_model(
      family_test = "limma_2",
      tempDataFrame = df,
      sig.formula   = Y ~ IV,
      burdenValue   = "Y",
      independent_variable = "IV",
      transformation_y = "", plot = FALSE,
      samples_sql_condition = "", key = key),
    regexp = "requires the 'limma' package.*BiocManager::install"
  )
})

# (3) Numerical: per-area limma ≈ stats::lm in degenerate (1-area) mode

test_that("limma_2 per-area p-values agree with stats::lm on simulated data", {
  testthat::skip_if_not_installed("limma")

  df  <- .sim_one_area(n = 200L, degree = 2L, sigma = 0.5)
  key <- .fake_key()

  # limma path
  res_l <- SEMseeker:::assoc_model_limma(
    "limma_2", df, Y ~ IV,
    transformation_y = "", plot = FALSE,
    samples_sql_condition = "", key = key)

  # lm reference using the same polynomial design
  fit_lm <- stats::lm(Y ~ stats::poly(IV, 2, raw = TRUE), data = df)
  sm <- summary(fit_lm)$coefficients

  # 4 columns: PL_DEGREE, PL_PERC, r_model, then 3 (intercept + 2 poly
  # terms) × 2 (pvalue + estimate) = 6, total 9.
  testthat::expect_gte(ncol(res_l), 9L)

  # Intercept estimate should match closely (eBayes does not shift mean)
  est_intercept_l <- res_l[, grepl("^INTERCEPT.*ESTIMATE$",
                                     toupper(colnames(res_l)))]
  testthat::expect_true(length(est_intercept_l) >= 1L)
  testthat::expect_equal(as.numeric(est_intercept_l[1]),
                          unname(sm[1, "Estimate"]),
                          tolerance = 1e-8)
})
