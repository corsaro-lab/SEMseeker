# Shared fixture for the burden / per-sample statistics tests.
# Lives in a helper because testthat isolates each test file's environment:
# test-8-burden-integration.R and test-0-sample-stats-sibling.R both use it.
#
# AI-254 — why every fixture is a PARAMETER and nothing is read from the
# surrounding scope.
#
# This function used to read `nsamples`, `probe_features`, `mySampleSheet` and
# `signal_data` as free variables, and it was the only helper in the suite doing
# so. That made it the only one whose result depended on HOW the tests were run:
#
#   devtools::test()  -> the closure of a helper is `package:SEMseeker`, whose
#                        parent chain is devtools_shims -> package:testthat ->
#                        ... -> base -> R_EmptyEnv. Neither R_GlobalEnv (where
#                        setup.R's `<<-` lands) nor setup.R's own source
#                        environment (where its `<-` lands) is on that chain, so
#                        NEITHER form of assignment reaches here and the function
#                        died on its first free variable.
#   R CMD check       -> the package is installed and attached, the chain does
#                        reach them, and the same code worked.
#
# The measured consequence was five tests that failed locally, passed in CI, and
# — worse — died before reaching a single assertion, hiding real defects behind
# an environment problem.
#
# helper-bedmethyl.R never had the issue because it takes everything it needs as
# arguments. That is the convention; this file had departed from it.
.burden_setup_signal_with_outliers <- function(seed = 12345L,
                                               n_samples,
                                               probe_features,
                                               sample_sheet,
                                               signal_data) {
  set.seed(seed)
  n_probes_b  <- 200L
  n_samples_b <- n_samples
  local_probes  <- probe_features[1:n_probes_b, ]
  local_samples <- sample_sheet

  # Background: mostly methylated (Beta(90, 10)).
  # IMPORTANT: outliers are detected per-probe across samples (IQR × 3 on
  # the inter-sample distribution). If we inject the SAME band across all
  # samples uniformly, the per-probe q1/q3 collapse to that uniform value
  # and IQR → 0 → no sample is flagged. So we inject DISJOINT, RANDOM
  # per-sample probe sets — each sample is an outlier at probes the other
  # samples are not, which is exactly the inter-sample dispersion the
  # threshold needs.
  local_sig <- matrix(
    stats::rbeta(n_probes_b * n_samples_b, 90L, 10L),
    nrow = n_probes_b, ncol = n_samples_b
  )
  for (s in seq_len(n_samples_b)) {
    # 10 random HYPO probes from positions 1:100 → values ≈ 0
    hypo_probes <- sample.int(100L, 10L)
    local_sig[hypo_probes, s] <- stats::rbeta(10L, 1L, 100L)
    # 10 random HYPER probes from positions 101:200 → values ≈ 1
    hyper_probes <- 100L + sample.int(100L, 10L)
    local_sig[hyper_probes, s] <- stats::rbeta(10L, 100L, 1L)
  }

  rownames(local_sig) <- local_probes$PROBE
  local_sig <- as.data.frame(local_sig)
  # signal_data has 10 unique columns; mySampleSheet has 16 rows (Reference reuse pattern)
  colnames(local_sig) <- colnames(signal_data)
  list(signal = local_sig, samples = local_samples, probes = local_probes)
}

.burden_required_markers <- c("MUTATIONS",
                              "DELTAP", "DELTAQ", "DELTARP", "DELTARQ",
                              "DELTAS", "DELTAR")
# AI-223: the burden lives in the statistics sibling under the SAMPLE scope,
# so the column names carry the scope prefix. Composed here the same way
# io_burden_colname() composes them on the package side.
# AI-248: the name carries the aggregation, and which one is produced by
# default depends on the class of the marker — a count is summed, a continuous
# deviation carries the descriptor set (of which MEAN is the historical value).
.burden_discrete_markers   <- c("MUTATIONS", "DELTAP", "DELTAQ", "DELTARP", "DELTARQ")
.burden_continuous_markers <- c("DELTAS", "DELTAR")
.burden_cols <- c(
  paste0("SAMPLE_", c(paste0(.burden_discrete_markers, "_HYPER"),
                      paste0(.burden_discrete_markers, "_HYPO")), "_SUM"),
  paste0("SAMPLE_", c(paste0(.burden_continuous_markers, "_HYPER"),
                      paste0(.burden_continuous_markers, "_HYPO")), "_MEAN"))
