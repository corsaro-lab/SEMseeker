## Session-based coverage for assoc_glm_model() (previously untested).
##
## assoc_glm_model() fits stats::glm() for a given family and returns a
## one-row data.frame of AIC/BIC, per-coefficient p-values and estimates,
## model-performance metrics and (gaussian only) residual diagnostics.
## Uses synthetic data — this exercises the statistical plumbing, not a
## biological claim, so a controlled random design is appropriate.
##

.glm_key <- function() {
  list(AREA = "GENE", SUBAREA = "TSS200", MARKER = "MUTATIONS", FIGURE = "K850")
}

test_that("assoc_glm_model (gaussian) returns AIC/BIC, p-values and residual diagnostics", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  set.seed(1)
  df  <- data.frame(x = stats::rnorm(40), y = stats::rnorm(40))

  res <- SEMseeker:::assoc_glm_model(
    family_test           = "gaussian",
    tempDataFrame         = df,
    sig.formula           = stats::as.formula("y ~ x"),
    transformation_y      = "",
    plot                  = FALSE,
    samples_sql_condition = "",
    key                   = .glm_key()
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_true(all(c("AIC_VALUE", "BIC_VALUE") %in% colnames(res)))
  expect_type(res$AIC_VALUE, "double")
  expect_true(any(grepl("PVALUE",  colnames(res))))    # per-coefficient p-values
  expect_true(any(grepl("ESTIMATE", colnames(res))))
  expect_true("residuals_normality" %in% colnames(res))
  expect_equal(res$r_model, "stats::glm")
})

test_that("assoc_glm_model (binomial) fits a logistic model and reports AIC", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  set.seed(2)
  n  <- 40L
  x  <- stats::rnorm(n)
  df <- data.frame(x = x, y = stats::rbinom(n, 1, stats::plogis(x)))

  res <- SEMseeker:::assoc_glm_model(
    family_test           = "binomial",
    tempDataFrame         = df,
    sig.formula           = stats::as.formula("y ~ x"),
    transformation_y      = "",
    plot                  = FALSE,
    samples_sql_condition = "",
    key                   = .glm_key()
  )

  expect_s3_class(res, "data.frame")
  expect_true("AIC_VALUE" %in% colnames(res))
  expect_true(any(grepl("PVALUE", colnames(res))))
  expect_equal(res$r_model, "stats::glm")
})
