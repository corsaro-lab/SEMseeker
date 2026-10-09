# =============================================================================
# transformation_x reaches the model, and the burden columns are found by name.
#
# transformation_x used to be applied in io_data_preparation() and lost there:
# the transformed values were written into the table, and the table was then
# rebuilt from a copy taken before the transformation. Only "scale" worked, by
# another road (a <IV>_SCALED column made upstream). It is now applied upstream
# for every transformation, into a column of its own, as each covariate is.
#
# Expected values are computed here with base R on the input, not by running
# the package's transformation code again.
# =============================================================================

.x_detail <- function(transformation_x, family_test = "polynomial_2_1") {
  list(collinearity_check = FALSE, covariates_dummy = "", covariates_pca = FALSE,
       covariates = "", independent_variable = "STAGE",
       transformation_x = transformation_x, family_test = family_test,
       transformation_y = "none", samples_sql_condition = NULL)
}

.x_summary <- function() {
  data.frame(Sample_ID = paste0("S", seq_len(12)),
             STAGE = rep(c(1, 2, 4, 8), 3),
             stringsAsFactors = FALSE)
}

.x_session <- function(env = parent.frame()) {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE,
                            parallel_strategy = "sequential")
  withr::defer({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) },
               envir = env)
  tf
}

test_that("transformation_x makes a column of its own and the model is pointed at it", {
  .x_session()
  ss <- .x_summary()

  res <- SEMseeker:::assoc_predictors_prepare(.x_detail("log10"), ss)
  expect_equal(res$inference_detail$independent_variable, "STAGE_LOG10")
  expect_equal(res$study_summary$STAGE_LOG10, log10(ss$STAGE))

  res <- SEMseeker:::assoc_predictors_prepare(.x_detail("pow_2"), ss)
  expect_equal(res$inference_detail$independent_variable, "STAGE_POW2")
  expect_equal(res$study_summary$STAGE_POW2, ss$STAGE^2)

  # scale keeps the name it always had
  res <- SEMseeker:::assoc_predictors_prepare(.x_detail("scale"), ss)
  expect_equal(res$inference_detail$independent_variable, "STAGE_SCALED")
  expect_equal(res$study_summary$STAGE_SCALED,
               (ss$STAGE - mean(ss$STAGE)) / stats::sd(ss$STAGE))

  # none leaves the request alone
  res <- SEMseeker:::assoc_predictors_prepare(.x_detail("none"), ss)
  expect_equal(res$inference_detail$independent_variable, "STAGE")
})

test_that("transformation_x = factor is refused for a family that fits a number", {
  base <- data.frame(independent_variable = "STAGE", transformation_y = "none",
                     transformation_x = "factor", stringsAsFactors = FALSE)

  regression <- base; regression$family_test <- "gaussian"
  expect_error(SEMseeker:::assoc_validate_transformation(regression),
               "fits the independent variable as a number")

  polynomial <- base; polynomial$family_test <- "polynomial_2_1"
  expect_error(SEMseeker:::assoc_validate_transformation(polynomial),
               "fits the independent variable as a number")

  # a k-group family takes the variable as categorical already: accepted
  kgroup <- base; kgroup$family_test <- "kruskal.test"
  expect_silent(SEMseeker:::assoc_validate_transformation(kgroup))
})

test_that("a single burden column keeps its name", {
  .x_session()
  df <- data.frame(AGE = c(30, 40, 50, 60, 70), MEDIAN = c(1, 3, 2, 5, 4))
  key <- data.frame(MARKER = "SIGNAL", FIGURE = "BETA", AREA = "PROBE", SUBAREA = "WHOLE")

  out <- SEMseeker:::io_data_preparation("spearman", "none", df, "AGE", 2L, 2L,
                                         NULL, key)
  expect_equal(colnames(out$tempDataFrame), c("AGE", "MEDIAN"))
  expect_equal(out$burden_columns, "MEDIAN")
})

test_that("a layout that does not match g_start:g_end is refused", {
  .x_session()
  # a covariate after the burdens: by position it would be tested as an area
  df <- data.frame(AGE = 1:6, GENE1 = c(1, 2, 2, 3, 5, 8), BMI = 6:1)
  key <- data.frame(MARKER = "MUTATIONS", FIGURE = "HYPO", AREA = "GENE", SUBAREA = "WHOLE")

  expect_error(
    SEMseeker:::io_data_preparation("gaussian", "none", df, "AGE", 2L, 3L, "BMI", key),
    "g_start:g_end points at")
})

test_that("NaN in a burden is found also when the independent variable is a factor", {
  tf <- .x_session()
  ssEnv <- SEMseeker:::core_get_session_info()
  # log of a negative burden is NaN; the frame also carries a factor IV, which
  # used to turn the whole-table NaN check into a character matrix that found
  # nothing.
  df <- data.frame(GROUP = c("a", "a", "a", "b", "b", "b"),
                   GENE1 = c(-5, 1, 2, 3, 4, 6), GENE2 = c(1, 2, 3, 4, 5, 6),
                   stringsAsFactors = FALSE)
  key <- data.frame(MARKER = "MUTATIONS", FIGURE = "HYPO", AREA = "GENE", SUBAREA = "WHOLE")

  out <- suppressWarnings(
    SEMseeker:::io_data_preparation("kruskal.test", "log", df, "GROUP", 2L, 3L, NULL, key))
  expect_true(is.factor(out$tempDataFrame$GROUP))
  expect_false(any(is.nan(out$tempDataFrame$GENE1)))
  expect_length(list.files(ssEnv$session_folder, pattern = "^lost_data_"), 1L)
})
