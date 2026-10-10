# 2026-06-09: polars-native equivalent of
# `io_data_preparation()` for the lazy batch path. Pins down:
#
#   1. transformation_y (none/log/log2/log10/exp/scale) produces values
#      bit-equal (tol 1e-7) to R-side log() / log10() / etc.
#   2. universal degenerate-burden filter drops rows where every
#      sample column carries the same value.
#   3. Unsupported transformation_y (factor, quantile_<N>) emits a
#      warning, falls back to "none", and continues — the lazy batch
#      path keeps running (no hard stop, by design).
#   4. The function ALWAYS returns a polars_lazy_frame (never materialises).

test_that("io_data_preparation_lazy passes 'none' through unchanged when no degenerate rows", {
  skip_if_not_installed("polars")
  set.seed(11L)
  mat <- matrix(stats::rnorm(50L * 8L), nrow = 50L, ncol = 8L)
  df <- data.frame(AREA = paste0("p_", seq_len(50L)), mat,
                    check.names = FALSE)
  colnames(df)[-1] <- paste0("S", sprintf("%02d", 1:8))
  pivot <- polars::as_polars_df(df)$lazy()

  out <- SEMseeker:::io_data_preparation_lazy(
    pivot_lazy        = pivot,
    sample_cols       = paste0("S", sprintf("%02d", 1:8)),
    transformation_y  = "none",
    apply_degenerate_filter = FALSE
  )
  expect_s3_class(out, "polars_lazy_frame")
  collected <- as.data.frame(out$collect())
  expect_equal(nrow(collected), 50L)
  expect_equal(collected[, -1L], df[, -1L], tolerance = 1e-12)
})

test_that("degenerate-burden filter drops rows with var(Y)==0", {
  skip_if_not_installed("polars")
  set.seed(12L)
  n <- 30L
  mat <- matrix(stats::rnorm(n * 5L), nrow = n, ncol = 5L)
  mat[3L, ]  <- 1.5      # degenerate (constant across all 5 samples)
  mat[17L, ] <- -2.7     # degenerate
  df <- data.frame(AREA = paste0("p_", seq_len(n)), mat,
                    check.names = FALSE)
  colnames(df)[-1] <- paste0("S", sprintf("%02d", 1:5))
  pivot <- polars::as_polars_df(df)$lazy()

  out <- SEMseeker:::io_data_preparation_lazy(
    pivot, paste0("S", sprintf("%02d", 1:5)),
    transformation_y = "none",
    apply_degenerate_filter = TRUE
  )
  collected <- as.data.frame(out$collect())
  expect_equal(nrow(collected), n - 2L)
  expect_false("p_3"  %in% collected$AREA)
  expect_false("p_17" %in% collected$AREA)
})

test_that("transformation_y=log10 matches R log10 bit-equal (tol 1e-9)", {
  skip_if_not_installed("polars")
  set.seed(13L)
  mat <- abs(matrix(stats::rnorm(20L * 4L), nrow = 20L, ncol = 4L)) + 1
  df <- data.frame(AREA = paste0("p_", seq_len(20L)), mat,
                    check.names = FALSE)
  colnames(df)[-1] <- paste0("S", sprintf("%02d", 1:4))
  pivot <- polars::as_polars_df(df)$lazy()

  out <- SEMseeker:::io_data_preparation_lazy(
    pivot, paste0("S", sprintf("%02d", 1:4)),
    transformation_y = "log10",
    apply_degenerate_filter = FALSE
  )
  collected <- as.data.frame(out$collect())

  # R-side reference: log10(x + 1e-9) — matches the EPS used inside
  # io_data_preparation_lazy() to avoid log10(0) → -Inf.
  expected <- log10(mat + 1e-9)
  expect_equal(unname(as.matrix(collected[, -1L])), expected,
                tolerance = 1e-9)
})

test_that("transformation_y=scale produces mean~0, sd~1 per sample col", {
  skip_if_not_installed("polars")
  set.seed(14L)
  mat <- matrix(stats::rnorm(40L * 3L, mean = 5, sd = 2),
                 nrow = 40L, ncol = 3L)
  df <- data.frame(AREA = paste0("p_", seq_len(40L)), mat,
                    check.names = FALSE)
  colnames(df)[-1] <- paste0("S", sprintf("%02d", 1:3))
  pivot <- polars::as_polars_df(df)$lazy()

  out <- SEMseeker:::io_data_preparation_lazy(
    pivot, paste0("S", sprintf("%02d", 1:3)),
    transformation_y = "scale",
    apply_degenerate_filter = FALSE
  )
  collected <- as.data.frame(out$collect())
  for (s in paste0("S", sprintf("%02d", 1:3))) {
    expect_equal(mean(collected[[s]]), 0, tolerance = 1e-9)
    expect_equal(stats::sd(collected[[s]]), 1, tolerance = 1e-9)
  }
})

test_that("unsupported transformation_y (factor) → warning + fallback to none", {
  skip_if_not_installed("polars")
  set.seed(15L)
  mat <- matrix(stats::rnorm(30L * 4L), nrow = 30L, ncol = 4L)
  df <- data.frame(AREA = paste0("p_", seq_len(30L)), mat,
                    check.names = FALSE)
  colnames(df)[-1] <- paste0("S", sprintf("%02d", 1:4))
  pivot <- polars::as_polars_df(df)$lazy()

  expect_warning(
    out <- SEMseeker:::io_data_preparation_lazy(
      pivot, paste0("S", sprintf("%02d", 1:4)),
      transformation_y = "factor",
      apply_degenerate_filter = FALSE
    ),
    regexp = "factor.*NOT supported"
  )
  # Output is the input unchanged (transformation skipped, no degenerate filter)
  collected <- as.data.frame(out$collect())
  expect_equal(nrow(collected), 30L)
  expect_equal(collected[, -1L], df[, -1L], tolerance = 1e-12)
})

test_that("unsupported quantile_4 → warning + fallback to none", {
  skip_if_not_installed("polars")
  set.seed(16L)
  mat <- matrix(stats::rnorm(20L * 3L), nrow = 20L, ncol = 3L)
  df <- data.frame(AREA = paste0("p_", seq_len(20L)), mat,
                    check.names = FALSE)
  colnames(df)[-1] <- paste0("S", sprintf("%02d", 1:3))
  pivot <- polars::as_polars_df(df)$lazy()

  expect_warning(
    SEMseeker:::io_data_preparation_lazy(
      pivot, paste0("S", sprintf("%02d", 1:3)),
      transformation_y = "quantile_4"
    ),
    regexp = "quantile_4.*NOT supported"
  )
})

test_that("io_data_preparation_lazy refuses non-Polars input gracefully", {
  skip_if_not_installed("polars")
  expect_error(
    SEMseeker:::io_data_preparation_lazy(
      pivot_lazy = data.frame(S1 = 1:5),
      sample_cols = "S1"
    ),
    regexp = "polars LazyFrame or DataFrame"
  )
})

test_that("io_data_preparation_lazy refuses empty sample_cols", {
  skip_if_not_installed("polars")
  df <- data.frame(AREA = "a", S1 = 1)
  pivot <- polars::as_polars_df(df)$lazy()
  expect_error(
    SEMseeker:::io_data_preparation_lazy(
      pivot_lazy = pivot,
      sample_cols = character(0)
    ),
    regexp = "non-empty character vector"
  )
})
