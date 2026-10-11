# 2026-06-09: KNN memory gate.
#
# What this file pins down:
#   1. .sem_knn_memory_gate() does NOT block when the estimated memory fits
#      under the available × SEMSEEKER_KNN_MEM_FRACTION budget.
#   2. .sem_knn_memory_gate() throws an actionable error when the matrix is
#      too large for the budget. The error message includes the four
#      breakdown numbers (matrix, distance, buffers, total) and the
#      fall-back options (median/mean/upstream KNN).
#   3. SEMSEEKER_KNN_MEM_FRACTION (env var) overrides the default 0.6.
#   4. .sem_total_ram_GB() returns a positive number on this platform.

test_that(".sem_total_ram_GB returns a positive number on all platforms", {
  ram <- SEMseeker:::.sem_total_ram_GB()
  expect_true(is.numeric(ram))
  expect_false(is.na(ram))
  expect_gt(ram, 0)
})

test_that(".sem_knn_memory_gate passes silently when matrix fits in budget", {
  # 100 × 50 doubles = 40 kB matrix; trivially under 60% of any real RAM.
  mat <- matrix(stats::rnorm(100L * 50L), nrow = 100L, ncol = 50L)
  expect_silent(SEMseeker:::.sem_knn_memory_gate(mat))
})

test_that(".sem_knn_memory_gate fails when SEMSEEKER_KNN_MEM_FRACTION is set very low", {
  mat <- matrix(stats::rnorm(1000L * 100L), nrow = 1000L, ncol = 100L)
  withr::with_envvar(
    new = c(SEMSEEKER_KNN_MEM_FRACTION = "0.00000001"),  # effectively no budget
    code = {
      expect_error(SEMseeker:::.sem_knn_memory_gate(mat),
                   regexp = "KNN imputation needs")
      err <- tryCatch(SEMseeker:::.sem_knn_memory_gate(mat), error = function(e) e)
      msg <- conditionMessage(err)
      # Error message includes the four breakdown numbers
      expect_match(msg, "matrix=", fixed = TRUE)
      expect_match(msg, "distance=", fixed = TRUE)
      expect_match(msg, "buffers=", fixed = TRUE)
      # ...and the fall-back options
      expect_match(msg, "median", fixed = TRUE)
      expect_match(msg, "upstream", fixed = TRUE)
    }
  )
})

test_that("SEMSEEKER_KNN_MEM_FRACTION default falls back to 0.6 on invalid values", {
  # Negative, zero, NA-string → fallback to 0.6 (not block on a real fit).
  mat <- matrix(stats::rnorm(100L * 50L), nrow = 100L, ncol = 50L)
  for (val in c("-0.5", "0", "abc", "")) {
    withr::with_envvar(
      new = c(SEMSEEKER_KNN_MEM_FRACTION = val),
      code = expect_silent(SEMseeker:::.sem_knn_memory_gate(mat))
    )
  }
})
