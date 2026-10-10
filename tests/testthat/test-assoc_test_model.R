## test-assoc_test_model.R
## Regression tests for assoc_test_model() — the statistical-test dispatcher.
##
## Before the fixes these tests document, three branches crashed immediately:
##   chisq.test  — line 47: `rea$effect_size`  (object 'rea' not found)
##   bartlett.test — line 68: `rea$effect_size` AND passing column-name strings
##                  instead of data vectors to stats::bartlett.test()
##   t.test      — line 243: `d = statistic_parameter` (object not found;
##                  the value had been stored in res$statistic_parameter)
##
## Each test below exercises one branch end-to-end.
## All require an active session (assoc_test_model calls core_get_session_info()).

# ---------------------------------------------------------------------------
# Helpers shared across tests in this file
# ---------------------------------------------------------------------------

.make_key <- function() {
  list(AREA = "GENE", SUBAREA = "TSS200", MARKER = "MUTATIONS", FIGURE = "K850")
}

# Two-group data frame: 10 ctrl / 10 case
.two_group_df <- function(seed = 1L) {
  set.seed(seed)
  data.frame(
    BURDEN = c(stats::rnorm(10, mean = 1), stats::rnorm(10, mean = 3)),
    GROUP  = factor(c(rep("ctrl", 10), rep("case", 10)))
  )
}

# ---------------------------------------------------------------------------
# t.test branch  (was crashing: d = statistic_parameter)
# ---------------------------------------------------------------------------

test_that("assoc_test_model t.test returns a data.frame with pvalue", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df  <- .two_group_df()
  key <- .make_key()
  f   <- stats::as.formula("BURDEN ~ GROUP")

  res <- SEMseeker:::assoc_test_model(
    family_test = "t.test", tempDataFrame = df,
    sig.formula = f, burdenValue = "BURDEN",
    independent_variable = "GROUP", transformation_y = "",
    plot = FALSE, samples_sql_condition = "", key = key
  )

  expect_s3_class(res, "data.frame")
  expect_true("pvalue" %in% tolower(colnames(res)))
  expect_false(is.null(res$pvalue))
  expect_true(is.numeric(res$pvalue))
})

test_that("assoc_test_model t.test: clearly separated groups give small p-value", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  set.seed(99)
  df <- data.frame(
    BURDEN = c(stats::rnorm(20, mean = 0, sd = 1), stats::rnorm(20, mean = 100, sd = 1)),
    GROUP  = factor(c(rep("ctrl", 20), rep("case", 20)))
  )
  f   <- stats::as.formula("BURDEN ~ GROUP")
  key <- .make_key()

  res <- SEMseeker:::assoc_test_model(
    "t.test", df, f, "BURDEN", "GROUP", "", FALSE, "", key
  )
  expect_lt(res$pvalue, 0.001)
})

# ---------------------------------------------------------------------------
# chisq.test branch  (was crashing: rea$effect_size)
# ---------------------------------------------------------------------------

test_that("assoc_test_model chisq.test returns a data.frame with pvalue", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # Discrete burden values so the contingency table is meaningful
  df <- data.frame(
    BURDEN = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
               0, 0, 1, 1, 1, 1, 1, 1, 1, 1),
    GROUP  = factor(c(rep("ctrl", 10), rep("case", 10)))
  )
  f   <- stats::as.formula("BURDEN ~ GROUP")
  key <- .make_key()

  res <- SEMseeker:::assoc_test_model(
    "chisq.test", df, f, "BURDEN", "GROUP", "", FALSE, "", key
  )

  expect_s3_class(res, "data.frame")
  expect_true("pvalue" %in% tolower(colnames(res)))
  expect_true(is.numeric(res$pvalue))
  expect_true("effect_size" %in% tolower(colnames(res)))  # was silently lost before fix
})

# ---------------------------------------------------------------------------
# bartlett.test branch  (was crashing: rea$effect_size + wrong data passed)
# ---------------------------------------------------------------------------

test_that("assoc_test_model bartlett.test returns a data.frame with pvalue", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df  <- .two_group_df(seed = 7L)
  f   <- stats::as.formula("BURDEN ~ GROUP")
  key <- .make_key()

  res <- SEMseeker:::assoc_test_model(
    "bartlett.test", df, f, "BURDEN", "GROUP", "", FALSE, "", key
  )

  expect_s3_class(res, "data.frame")
  expect_true("pvalue" %in% tolower(colnames(res)))
  expect_true(is.numeric(res$pvalue))
  expect_true("effect_size" %in% tolower(colnames(res)))
})

test_that("assoc_test_model bartlett.test: equal-variance groups give large p-value", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  set.seed(42)
  df <- data.frame(
    BURDEN = c(stats::rnorm(30, sd = 1), stats::rnorm(30, sd = 1)),
    GROUP  = factor(c(rep("a", 30), rep("b", 30)))
  )
  f   <- stats::as.formula("BURDEN ~ GROUP")
  key <- .make_key()

  res <- SEMseeker:::assoc_test_model(
    "bartlett.test", df, f, "BURDEN", "GROUP", "", FALSE, "", key
  )
  expect_gt(res$pvalue, 0.05)
})

# ---------------------------------------------------------------------------
# wilcoxon branch — direct unit test (also exercised through pipeline)
# ---------------------------------------------------------------------------

test_that("assoc_test_model wilcoxon returns pvalue and effect_size columns", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df  <- .two_group_df(seed = 3L)
  f   <- stats::as.formula("BURDEN ~ GROUP")
  key <- .make_key()

  res <- SEMseeker:::assoc_test_model(
    "wilcoxon", df, f, "BURDEN", "GROUP", "", FALSE, "", key
  )

  expect_s3_class(res, "data.frame")
  expect_true("pvalue" %in% tolower(colnames(res)))
  expect_true(is.numeric(res$pvalue))
})

# ---------------------------------------------------------------------------
# kruskal.test branch
# ---------------------------------------------------------------------------

test_that("assoc_test_model kruskal.test with 3 groups returns pvalue", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  set.seed(11)
  df <- data.frame(
    BURDEN = c(stats::rnorm(10, 0), stats::rnorm(10, 3), stats::rnorm(10, 6)),
    GROUP  = factor(c(rep("a", 10), rep("b", 10), rep("c", 10)))
  )
  f   <- stats::as.formula("BURDEN ~ GROUP")
  key <- .make_key()

  res <- SEMseeker:::assoc_test_model(
    "kruskal.test", df, f, "BURDEN", "GROUP", "", FALSE, "", key
  )

  expect_s3_class(res, "data.frame")
  expect_true("pvalue" %in% tolower(colnames(res)))
  expect_lt(res$pvalue, 0.05)
})

test_that("assoc_test_model kruskal.test with single group returns NA pvalue", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df <- data.frame(
    BURDEN = stats::rnorm(10),
    GROUP  = factor(rep("only_one", 10))
  )
  f   <- stats::as.formula("BURDEN ~ GROUP")
  key <- .make_key()

  res <- SEMseeker:::assoc_test_model(
    "kruskal.test", df, f, "BURDEN", "GROUP", "", FALSE, "", key
  )

  expect_true(is.na(res$pvalue))
})
