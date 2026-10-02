# Tests for additional pure / deterministic helper functions
#
# Covered:
#  - util_join_values_to_thresholds()  Polars positional inner join (join_values_to_thresholds.R)
#  - sem_metrics_ranking()            ranking helper using metrics_properties (metrics_ranking.R)
#  - assoc_model_performance()          train+test path and overfitting detection (model_performance.R)
#  - assoc_compute_quantreg_permutation() quantile regression single-permutation draw

# ---------------------------------------------------------------------------
# 1. util_join_values_to_thresholds
# ---------------------------------------------------------------------------

# Helpers ----------------------------------------------------------------
.jvt_values <- function(chr, starts, vals) {
  data.frame(
    CHR   = as.character(chr),
    START = as.integer(starts),
    END   = as.integer(starts),
    VALUE = as.numeric(vals),
    stringsAsFactors = FALSE
  )
}

.jvt_thresholds <- function(chr, starts, inf, sup, extra = NULL) {
  d <- data.frame(
    CHR                        = as.character(chr),
    START                      = as.integer(starts),
    END                        = as.integer(starts),
    signal_inferior_thresholds = as.numeric(inf),
    signal_superior_thresholds = as.numeric(sup),
    stringsAsFactors           = FALSE
  )
  if (!is.null(extra)) d <- cbind(d, extra)
  d
}
# ------------------------------------------------------------------------

test_that("util_join_values_to_thresholds: full overlap returns all rows with correct columns", {
  v <- .jvt_values("1", 1:5 * 1000L, seq(0.1, 0.5, by = 0.1))
  t <- .jvt_thresholds("1", 1:5 * 1000L, rep(0.2, 5), rep(0.8, 5))
  res <- SEMseeker:::util_join_values_to_thresholds(v, t)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 5L)
  expect_true(all(c("CHR", "START", "END", "VALUE",
                    "signal_inferior_thresholds",
                    "signal_superior_thresholds") %in% colnames(res)))
})

test_that("util_join_values_to_thresholds: partial overlap — only shared positions returned", {
  v <- .jvt_values("1", c(1000L, 2000L, 3000L), c(0.5, 0.5, 0.5))
  t <- .jvt_thresholds("1", c(2000L, 3000L, 4000L), rep(0.2, 3), rep(0.8, 3))
  res <- SEMseeker:::util_join_values_to_thresholds(v, t)
  expect_equal(nrow(res), 2L)           # 2000 and 3000 shared
  expect_equal(sort(res$START), c(2000L, 3000L))
})

test_that("util_join_values_to_thresholds: zero overlap returns 0-row data.frame without crash", {
  v <- .jvt_values("chr1", c(1000L, 2000L), c(0.5, 0.6))
  t <- .jvt_thresholds("chr2", c(1000L, 2000L), rep(0.2, 2), rep(0.8, 2))
  res <- SEMseeker:::util_join_values_to_thresholds(v, t)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0L)
})

test_that("util_join_values_to_thresholds: VALUE column is preserved correctly", {
  v <- .jvt_values("1", c(1000L, 2000L), c(0.3, 0.9))
  t <- .jvt_thresholds("1", c(1000L, 2000L), c(0.2, 0.2), c(0.8, 0.8))
  res <- SEMseeker:::util_join_values_to_thresholds(v, t)
  res_sorted <- res[order(res$START), ]
  expect_equal(res_sorted$VALUE, c(0.3, 0.9), tolerance = 1e-10)
})

test_that("util_join_values_to_thresholds: threshold columns are preserved", {
  v <- .jvt_values("1", 1000L, 0.5)
  t <- .jvt_thresholds("1", 1000L, 0.1, 0.9)
  res <- SEMseeker:::util_join_values_to_thresholds(v, t)
  expect_equal(res$signal_inferior_thresholds, 0.1, tolerance = 1e-10)
  expect_equal(res$signal_superior_thresholds, 0.9, tolerance = 1e-10)
})

test_that("util_join_values_to_thresholds: CHR prefix is normalised — '1' matches 'chr1'", {
  # See E-13 in R/join_values_to_thresholds.R: strip_chr() normalises the CHR
  # column so bed-file values (chr1) match threshold values (1). The join is
  # semantically chr-prefix-insensitive, not exact-string.
  v <- .jvt_values("1",    1000L, 0.5)
  t <- .jvt_thresholds("chr1", 1000L, 0.1, 0.9)
  res <- SEMseeker:::util_join_values_to_thresholds(v, t)
  expect_equal(nrow(res), 1L)
})

test_that("util_join_values_to_thresholds: genuinely different CHR gives zero rows", {
  # The join still rejects CHR pairs that differ on more than the 'chr' prefix.
  v <- .jvt_values("1", 1000L, 0.5)
  t <- .jvt_thresholds("2", 1000L, 0.1, 0.9)
  res <- SEMseeker:::util_join_values_to_thresholds(v, t)
  expect_equal(nrow(res), 0L)
})

test_that("util_join_values_to_thresholds: large (genome-scale) position values work", {
  big_starts <- c(100000000L, 150000000L, 200000000L)
  v <- .jvt_values("1", big_starts, c(0.2, 0.5, 0.8))
  t <- .jvt_thresholds("1", big_starts, rep(0.1, 3), rep(0.9, 3))
  res <- SEMseeker:::util_join_values_to_thresholds(v, t)
  expect_equal(nrow(res), 3L)
})

test_that("util_join_values_to_thresholds: optional columns (iqr, q1, q3) are passed through", {
  v <- .jvt_values("1", c(1000L, 2000L), c(0.4, 0.6))
  t <- .jvt_thresholds("1", c(1000L, 2000L), c(0.1, 0.1), c(0.9, 0.9),
                        extra = data.frame(iqr = c(0.2, 0.2),
                                           q1  = c(0.3, 0.3),
                                           q3  = c(0.7, 0.7)))
  res <- SEMseeker:::util_join_values_to_thresholds(v, t)
  expect_equal(nrow(res), 2L)
  expect_true(all(c("iqr", "q1", "q3") %in% colnames(res)))
})

test_that("util_join_values_to_thresholds: extra columns not in keep-list are dropped", {
  v <- .jvt_values("1", 1000L, 0.5)
  t <- .jvt_thresholds("1", 1000L, 0.1, 0.9,
                        extra = data.frame(GENE_TSS200 = "BRCA1",  # annotation column
                                           stringsAsFactors = FALSE))
  res <- SEMseeker:::util_join_values_to_thresholds(v, t)
  # Annotation columns not in the keep_cols list must be absent
  expect_false("GENE_TSS200" %in% colnames(res))
})

test_that("util_join_values_to_thresholds: multiple chromosomes handled correctly", {
  starts <- rep(1000L, 6)
  v <- data.frame(
    CHR   = c("1","1","2","2","X","X"),
    START = starts,
    END   = starts,
    VALUE = seq(0.1, 0.6, by = 0.1),
    stringsAsFactors = FALSE
  )
  t <- data.frame(
    CHR                        = c("1","2","X"),
    START                      = 1000L,
    END                        = 1000L,
    signal_inferior_thresholds = 0.2,
    signal_superior_thresholds = 0.8,
    stringsAsFactors           = FALSE
  )
  res <- SEMseeker:::util_join_values_to_thresholds(v, t)
  # Two value rows per chromosome but only one threshold row → 3 match rows total
  # (polars inner join: each chr1 value with each chr1 threshold → cross-join within chr)
  # Actually: positions are unique per chr since all starts are 1000 → 2 values match 1 threshold each
  expect_equal(nrow(res), 6L)
})

# ---------------------------------------------------------------------------
# 2. sem_metrics_ranking
# ---------------------------------------------------------------------------

test_that("sem_metrics_ranking: returns data.frame with SCORE and METRIC columns", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df <- data.frame(REBASED = c(0.01, 0.05, 0.001, 0.20), stringsAsFactors = FALSE)
  result <- SEMseeker:::sem_metrics_ranking("PVALUE", df, column_to_rank = "REBASED")
  expect_s3_class(result, "data.frame")
  expect_true("SCORE"  %in% colnames(result))
  expect_true("METRIC" %in% colnames(result))
})

test_that("sem_metrics_ranking: PVALUE is lower-is-better — lower p-value gets higher score", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df     <- data.frame(REBASED = c(0.001, 0.05, 0.5, 0.9), stringsAsFactors = FALSE)
  result <- SEMseeker:::sem_metrics_ranking("PVALUE", df, column_to_rank = "REBASED")
  best_idx  <- which.min(df$REBASED)   # p = 0.001
  worst_idx <- which.max(df$REBASED)   # p = 0.9
  expect_gt(result$SCORE[best_idx], result$SCORE[worst_idx])
})

test_that("sem_metrics_ranking: R_SQUARED is higher-is-better — higher value gets higher score", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  skip_if(!("R_SQUARED" %in% toupper(SEMseeker::metrics_properties$Metric)),
          "R_SQUARED not in metrics_properties")

  df     <- data.frame(REBASED = c(0.1, 0.5, 0.9, 0.95), stringsAsFactors = FALSE)
  result <- SEMseeker:::sem_metrics_ranking("R_SQUARED", df, column_to_rank = "REBASED")
  best_idx  <- which.max(df$REBASED)   # 0.95
  worst_idx <- which.min(df$REBASED)   # 0.1
  expect_gt(result$SCORE[best_idx], result$SCORE[worst_idx])
})

test_that("sem_metrics_ranking: uniform input → all scores equal 1", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df     <- data.frame(REBASED = c(0.5, 0.5, 0.5), stringsAsFactors = FALSE)
  result <- SEMseeker:::sem_metrics_ranking("PVALUE", df, column_to_rank = "REBASED")
  expect_true(all(result$SCORE == 1))
})

test_that("sem_metrics_ranking: METRIC column is set to the requested metric name", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df     <- data.frame(REBASED = c(0.01, 0.05, 0.10), stringsAsFactors = FALSE)
  result <- SEMseeker:::sem_metrics_ranking("PVALUE", df, column_to_rank = "REBASED")
  expect_true(all(result$METRIC == "PVALUE"))
})

test_that("sem_metrics_ranking: PVALUE_ADJ is recognised as lower-is-better via grep", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df <- data.frame(REBASED = c(0.001, 0.01, 0.05, 0.20), stringsAsFactors = FALSE)
  # Any metric containing "PVALUE" (by grep) is treated as lower-is-better
  result_pv    <- SEMseeker:::sem_metrics_ranking("PVALUE",     df, column_to_rank = "REBASED")
  result_pvadj <- SEMseeker:::sem_metrics_ranking("PVALUE_ADJ", df, column_to_rank = "REBASED")
  # Both should rank the same way
  expect_equal(result_pv$SCORE, result_pvadj$SCORE)
})

test_that("sem_metrics_ranking: Inf values are replaced before ranking (no infinite scores)", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df <- data.frame(REBASED = c(0.01, Inf, 0.5), stringsAsFactors = FALSE)
  # Should not crash — Inf is replaced with 1E300 internally
  result <- SEMseeker:::sem_metrics_ranking("PVALUE", df, column_to_rank = "REBASED")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3L)
  expect_false(any(is.infinite(result$SCORE)))
})

# ---------------------------------------------------------------------------
# 3. assoc_model_performance: train+test path and overfitting detection
# (extends the basic tests in test-0-pure-functions.R)
# ---------------------------------------------------------------------------

test_that("assoc_model_performance: train+test path adds *_test columns", {
  train_fit <- c(1.0, 2.0, 3.0, 4.0)
  train_exp <- c(1.1, 1.9, 3.1, 3.9)
  test_fit  <- c(5.0, 6.0)
  test_exp  <- c(5.5, 6.5)
  res <- SEMseeker:::assoc_model_performance(train_fit, train_exp, test_fit, test_exp)
  expected_test_cols <- c("mse_test", "rmse_test", "mape_test",
                          "mae_test", "r_squared_test", "overfitting")
  expect_true(all(expected_test_cols %in% colnames(res)))
})

test_that("assoc_model_performance: rmse_test = sqrt(mse_test)", {
  train_fit <- c(1, 2, 3, 4, 5)
  train_exp <- c(1, 2, 3, 4, 5)
  test_fit  <- c(1.5, 2.5, 3.5)
  test_exp  <- c(1.0, 2.0, 3.0)
  res <- SEMseeker:::assoc_model_performance(train_fit, train_exp, test_fit, test_exp)
  expect_equal(res$rmse_test, sqrt(res$mse_test), tolerance = 1e-10)
})

test_that("assoc_model_performance: perfect train and test → overfitting = FALSE", {
  x <- c(1, 2, 3, 4, 5)
  res <- SEMseeker:::assoc_model_performance(x, x, x, x)
  expect_equal(res$mse,      0)
  expect_equal(res$mse_test, 0)
  expect_false(res$overfitting)
})

test_that("assoc_model_performance: terrible test performance → overfitting = TRUE", {
  # Perfect fit on train, reversed fit on test
  train_fit <- c(1, 2, 3, 4, 5)
  train_exp <- c(1, 2, 3, 4, 5)  # mse_train = 0
  test_fit  <- c(1, 2, 3, 4, 5)
  test_exp  <- c(5, 4, 3, 2, 1)  # mse_test >> 0
  res <- SEMseeker:::assoc_model_performance(train_fit, train_exp, test_fit, test_exp)
  expect_true(res$overfitting)
})

test_that("assoc_model_performance: r_squared_test in [0, 1] for near-perfect predictions", {
  train_fit <- c(1.0, 2.0, 3.0, 4.0)
  train_exp <- c(1.0, 2.0, 3.0, 4.0)
  test_fit  <- c(1.1, 1.9, 3.1, 3.9)
  test_exp  <- c(1.0, 2.0, 3.0, 4.0)
  res <- SEMseeker:::assoc_model_performance(train_fit, train_exp, test_fit, test_exp)
  expect_gte(res$r_squared_test, 0.9)
  expect_lte(res$r_squared_test, 1.0)
})

test_that("assoc_model_performance: mse_test reflects actual squared residuals on test set", {
  # Use varied training data to avoid 0/0 → NaN in r_squared computation
  train_fit <- c(1, 2, 3, 4)
  train_exp <- c(1, 2, 3, 4)  # perfect train, non-constant → r_squared = 1
  test_fit  <- c(2, 3, 4, 5)
  test_exp  <- c(1, 2, 3, 4)  # each prediction off by 1 → mse_test = 1
  res <- SEMseeker:::assoc_model_performance(train_fit, train_exp, test_fit, test_exp)
  expect_equal(res$mse_test, 1.0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 4. assoc_compute_quantreg_permutation
# ---------------------------------------------------------------------------

test_that("assoc_compute_quantreg_permutation: returns a numeric scalar", {
  skip_if_not_installed("lqmm")
  set.seed(42)
  n  <- 30
  df <- data.frame(
    burden = c(stats::rnorm(n/2, 1), stats::rnorm(n/2, 3)),
    age    = stats::rnorm(n, 50, 10)
  )
  f       <- stats::as.formula("burden ~ age")
  control <- list(loop_tol_ll = 1e-5, loop_max_iter = 10000, verbose = FALSE)
  result  <- SEMseeker:::assoc_compute_quantreg_permutation(f, df, tau = 0.5, lqm_control = control)
  expect_length(result, 1)
  expect_true(is.numeric(result))
})

test_that("assoc_compute_quantreg_permutation: shuffle changes the coefficient", {
  skip_if_not_installed("lqmm")
  set.seed(7)
  n  <- 40
  df <- data.frame(
    burden = 1:n + stats::rnorm(n, 0, 0.5),
    pheno  = 1:n
  )
  f       <- stats::as.formula("burden ~ pheno")
  control <- list(loop_tol_ll = 1e-5, loop_max_iter = 10000, verbose = FALSE)

  obs   <- SEMseeker:::assoc_compute_quantreg_permutation(f, df, tau = 0.5, lqm_control = control)
  perms <- replicate(20, SEMseeker:::assoc_compute_quantreg_permutation(f, df, tau = 0.5,
                                                                  lqm_control = control))
  # Shuffled permutations should not all equal the observed coefficient
  expect_false(all(perms == obs))
})
