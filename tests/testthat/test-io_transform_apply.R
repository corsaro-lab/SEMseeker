# io_transform_apply() / io_transform_suffix(): the vocabulary of transformations,
# and one transformation per covariate.
#
# The vocabulary was written twice, for the dependent and for the independent
# variable, with quantile_<n> handled beside the second copy rather than inside
# it. These pin it as one list, in one place, before a third caller uses it.

.ta <- function(...) SEMseeker:::io_transform_apply(...)
.ts <- function(...) SEMseeker:::io_transform_suffix(...)

# ---- the vocabulary ---------------------------------------------------------

test_that("io_transform_apply: each name does what it says, on a vector", {
  v <- c(1, 2, 4, 8)
  expect_equal(.ta(v, "none"),  v)
  expect_equal(.ta(v, "log"),   log(v))
  expect_equal(.ta(v, "log2"),  log2(v))
  expect_equal(.ta(v, "log10"), log10(v))
  expect_equal(.ta(v, "exp"),   exp(v))
  expect_equal(as.numeric(.ta(v, "scale")), as.numeric(scale(v)))
  expect_s3_class(.ta(c("a", "b"), "factor"), "factor")
})

test_that("io_transform_apply: absent, NA and empty all mean none", {
  v <- c(1, 2, 3)
  for (absent in list(NULL, NA, NA_character_, "", "  "))
    expect_equal(.ta(v, absent), v, info = paste(format(absent), collapse = ""))
})

test_that("io_transform_apply: pow_<n> raises to the power, fractional included", {
  v <- c(1, 2, 3, 4)
  expect_equal(.ta(v, "pow_2"),   v^2)
  expect_equal(.ta(v, "pow_3"),   v^3)
  expect_equal(.ta(v, "pow_0.5"), sqrt(v))
  # pow_1 is the identity and is allowed: a request that says it is explicit
  # about not transforming, which is not the same statement as "none".
  expect_equal(.ta(v, "pow_1"), v)
})

test_that("io_transform_apply: quantile_<n> keeps the behaviour it had before", {
  v <- c(1, 2, 3, 4, 5, 6)
  expect_equal(.ta(v, "quantile_3"), as.numeric(dplyr::ntile(v, 3)))
  # Fewer distinct values than tiles gave zeros where this lived before, beside
  # the switch rather than in it. Kept, so a request that used it reads the same.
  expect_equal(.ta(c(5, 5, 5), "quantile_3"), rep(0, 3))
})

test_that("io_transform_apply: a data.frame is transformed column-wise", {
  d <- data.frame(a = c(1, 2, 4), b = c(10, 100, 1000))
  got <- .ta(d, "log10")
  expect_equal(got$a, log10(d$a))
  expect_equal(got$b, c(1, 2, 3))
  # the dependent variable arrives as one column per sample, and quantile has to
  # tile each column on its own rather than the pooled values
  qd <- .ta(data.frame(a = c(1, 2, 3), b = c(30, 20, 10)), "quantile_3")
  expect_equal(qd$a, c(1, 2, 3))
  expect_equal(qd$b, c(3, 2, 1))
})

# ---- an unknown name is refused, where it used to be ignored ----------------

test_that("io_transform_apply: an unknown name is refused and the known ones listed", {
  # It used to fall through switch()'s default and return the values untouched.
  expect_error(.ta(c(1, 2), "lgo10"), "is not one this package knows")
  expect_error(.ta(c(1, 2), "lgo10"), "pow_<n>", fixed = TRUE)
})

test_that("io_transform_apply: a malformed parameter is refused, not guessed", {
  expect_error(.ta(c(1, 2), "pow_x"),      "needs a number")
  expect_error(.ta(c(1, 2), "quantile_x"), "needs an")
  expect_error(.ta(c(1, 2), "quantile_1"), "2 or more")
})

# ---- the suffix -------------------------------------------------------------

test_that("io_transform_suffix: the name of the column says the transformation", {
  expect_equal(.ts("exp"),        "EXP")
  expect_equal(.ts("log10"),      "LOG10")
  expect_equal(.ts("pow_2"),      "POW2")
  expect_equal(.ts("quantile_3"), "QUANTILE3")
  # scale keeps the suffix it had before this function existed, because
  # io_inference_file_name() and enrich_phenotype_analysis_name() both strip and
  # re-add it by that exact name.
  expect_equal(.ts("scale"), "SCALED")
  # none makes no column, so it has no suffix
  expect_true(is.na(.ts("none")))
  expect_true(is.na(.ts("")))
})

# ---- one transformation per covariate ---------------------------------------

.cov_detail <- function(covariates, transformation = NULL, transformation_x = "none") {
  d <- list(collinearity_check = FALSE, covariates_dummy = "", covariates_pca = FALSE,
            covariates = covariates, independent_variable = "STAGE",
            transformation_x = transformation_x, family_test = "polynomial_2_1",
            transformation_y = "none", samples_sql_condition = NULL)
  if (!is.null(transformation)) d$covariates_transformation <- transformation
  d
}

.cov_summary <- function() {
  set.seed(3)
  data.frame(Sample_ID = paste0("S", seq_len(20)),
             STAGE  = rep(1:4, 5),
             BURDEN = stats::rnorm(20, 10),
             AGE    = stats::runif(20, 40, 80),
             BMI    = stats::runif(20, 20, 32),
             SMOKING = rep(c(0, 1), 10),
             stringsAsFactors = FALSE)
}

test_that("assoc_predictors_prepare: each covariate takes its own transformation", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  ss <- .cov_summary()
  res <- SEMseeker:::assoc_predictors_prepare(
    .cov_detail("AGE + BMI + SMOKING", "exp + pow_2 + none"), ss)

  # the name moves with the values, so the model's coefficient is reported under
  # the name of the quantity it belongs to
  expect_equal(res$covariates, c("AGE_EXP", "BMI_POW2", "SMOKING"))
  expect_true(all(c("AGE_EXP", "BMI_POW2") %in% colnames(res$study_summary)))
  expect_equal(res$study_summary$AGE_EXP,  exp(ss$AGE))
  expect_equal(res$study_summary$BMI_POW2, ss$BMI^2)
  # the untransformed one is untouched and keeps its name
  expect_equal(res$study_summary$SMOKING, ss$SMOKING)
})

test_that("assoc_predictors_prepare: two lengths that disagree are refused, both named", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # R recycles silently, which would apply exp to SMOKING and name the column
  # for the covariate it came from.
  expect_error(
    SEMseeker:::assoc_predictors_prepare(
      .cov_detail("AGE + BMI + SMOKING", "exp + exp"), .cov_summary()),
    "has 2 entries and covariates has 3")
})

test_that("assoc_predictors_prepare: transformation_x no longer reaches the covariates", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  ss  <- .cov_summary()
  res <- SEMseeker:::assoc_predictors_prepare(
    .cov_detail("AGE + BMI", transformation_x = "scale"), ss)

  # the independent variable is scaled and renamed, as before
  expect_true("STAGE_SCALED" %in% colnames(res$study_summary))
  # the covariates are not, which they were: a transformation asked for one
  # variable was being applied to several
  expect_equal(res$covariates, c("AGE", "BMI"))
  expect_false(any(c("AGE_SCALED", "BMI_SCALED") %in% colnames(res$study_summary)))
})

test_that("assoc_predictors_prepare: running twice is the same as running once", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  d <- .cov_detail("AGE + BMI", "exp + pow_2")
  once  <- SEMseeker:::assoc_predictors_prepare(d, .cov_summary())
  twice <- SEMseeker:::assoc_predictors_prepare(d, once$study_summary)

  expect_equal(twice$covariates, once$covariates)
  expect_equal(sort(colnames(twice$study_summary)), sort(colnames(once$study_summary)))
  expect_equal(twice$study_summary$AGE_EXP, once$study_summary$AGE_EXP)
})

# ---- the door ---------------------------------------------------------------
# Refusing inside io_transform_apply() is not enough for the dependent and
# independent variables: those calls sit inside try(), which swallows the error
# and leaves the old behaviour. The refusal has to happen before that.

.req <- function(...) {
  d <- list(independent_variable = "STAGE", family_test = "kendall",
            aggregation = "SUM", scope = "SAMPLE",
            transformation_y = "none", transformation_x = "none")
  as.data.frame(utils::modifyList(d, list(...)), stringsAsFactors = FALSE)
}
.vt <- function(...) SEMseeker:::assoc_validate_transformation(...)

test_that("assoc_validate_transformation: the known names pass, including parameterised", {
  for (t in c("none", "scale", "log", "log2", "log10", "exp", "factor",
              "pow_2", "pow_0.5", "quantile_4"))
    expect_s3_class(.vt(.req(transformation_y = t)), "data.frame")
})

test_that("assoc_validate_transformation: an unknown name is refused at the door", {
  expect_error(.vt(.req(transformation_y = "lgo10")), "not a transformation this package knows")
  expect_error(.vt(.req(transformation_x = "sclae")), "transformation_x")
  # the message names the row, so a multi-row request says which one
  expect_error(.vt(.req(transformation_y = "nope")), "row 1")
})

test_that("assoc_validate_transformation: a malformed parameter is refused", {
  expect_error(.vt(.req(transformation_y = "pow_x")),      "not a transformation")
  expect_error(.vt(.req(transformation_y = "quantile_1")), "not a transformation")
})

test_that("assoc_validate_transformation: covariates_transformation is checked entry by entry", {
  expect_s3_class(.vt(.req(covariates = "AGE + BMI",
                           covariates_transformation = "exp + pow_2")), "data.frame")
  # the second entry is the bad one, and duplicates are not collapsed on the way
  expect_error(.vt(.req(covariates = "AGE + BMI",
                        covariates_transformation = "exp + lgo10")), "lgo10")
  expect_s3_class(.vt(.req(covariates = "AGE + BMI",
                           covariates_transformation = "exp + exp")), "data.frame")
})

test_that("assoc_validate_transformation: factor on a covariate is refused by name", {
  # Admissible for y and x, and not here: covariates_dummy encodes a categorical
  # covariate, where relabelling it as a factor does not.
  expect_s3_class(.vt(.req(transformation_y = "factor")), "data.frame")
  expect_error(.vt(.req(covariates = "SMOKING", covariates_transformation = "factor")),
               "covariates_dummy")
})

test_that("assoc_validate_transformation: a request that names no transformation passes", {
  d <- as.data.frame(list(independent_variable = "STAGE", family_test = "kendall",
                          aggregation = "SUM", scope = "SAMPLE"),
                     stringsAsFactors = FALSE)
  expect_s3_class(.vt(d), "data.frame")
  expect_null(.vt(NULL))
})
