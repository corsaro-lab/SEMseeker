# =============================================================================
# The thresholds are Tukey's fences, capped at [0, 1] for beta values only.
#
# From February 2024 the fences were clamped to the minimum and maximum the
# Reference observed. The intent, in that commit, was to keep the thresholds
# inside the range of the signal; the extremes of a handful of controls are
# much tighter than that range. With six controls the clamp was the effective
# threshold on every KCNQ1OT1 probe, and a healthy father barely below all six
# came out with epimutations at most imprinting regions.
#
# The expected values are computed here with stats::quantile on the input.
# =============================================================================

.fence_session <- function(beta, env = parent.frame()) {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE,
                            parallel_strategy = "sequential")
  ssEnv <- SEMseeker:::core_get_session_info()
  ssEnv$beta <- beta
  SEMseeker:::core_update_session_info(ssEnv)
  withr::defer({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, envir = env)
}

.fence_features <- function(probes) {
  data.frame(PROBE = probes, CHR = "1", START = seq_along(probes) * 100L,
             END = seq_along(probes) * 100L, stringsAsFactors = FALSE)
}

test_that("the fences are Q1 - 3 IQR and Q3 + 3 IQR, not the Reference extremes", {
  .fence_session(beta = TRUE)
  # one probe whose fences fall well inside (0, 1) and outside the observed range
  ref <- matrix(c(0.50, 0.51, 0.52, 0.53, 0.54, 0.55), nrow = 1,
                dimnames = list("cg_mid", paste0("R", 1:6)))
  th <- SEMseeker:::sem_signal_range_values(ref, "fences_mid", .fence_features("cg_mid"))

  q <- stats::quantile(ref[1, ], c(0.25, 0.75), names = FALSE)
  expect_equal(th$signal_inferior_thresholds, q[1] - 3 * (q[2] - q[1]))
  expect_equal(th$signal_superior_thresholds, q[2] + 3 * (q[2] - q[1]))
  # the clamp would have put them at 0.50 and 0.55
  expect_lt(th$signal_inferior_thresholds, min(ref))
  expect_gt(th$signal_superior_thresholds, max(ref))
})

test_that("for beta values the fences are capped at 0 and 1", {
  .fence_session(beta = TRUE)
  ref <- rbind(cg_low  = c(0.01, 0.02, 0.05, 0.10, 0.15, 0.20),
               cg_high = c(0.80, 0.85, 0.90, 0.95, 0.98, 0.99))
  colnames(ref) <- paste0("R", 1:6)
  th <- SEMseeker:::sem_signal_range_values(ref, "fences_cap", .fence_features(rownames(ref)))
  th <- th[match(rownames(ref), th$PROBE), ]

  # Tukey goes below 0 for the first probe and above 1 for the second
  expect_equal(th$signal_inferior_thresholds[1], 0)
  expect_equal(th$signal_superior_thresholds[2], 1)
  expect_true(all(th$signal_inferior_thresholds >= 0 & th$signal_superior_thresholds <= 1))
})

test_that("for M-values the fences are not capped", {
  .fence_session(beta = FALSE)
  ref <- matrix(c(-4, -3.5, -3, -2.5, -2, -1), nrow = 1,
                dimnames = list("cg_m", paste0("R", 1:6)))
  th <- SEMseeker:::sem_signal_range_values(ref, "fences_m", .fence_features("cg_m"))

  q <- stats::quantile(ref[1, ], c(0.25, 0.75), names = FALSE)
  expect_equal(th$signal_inferior_thresholds, q[1] - 3 * (q[2] - q[1]))
  expect_lt(th$signal_inferior_thresholds, 0)
})
