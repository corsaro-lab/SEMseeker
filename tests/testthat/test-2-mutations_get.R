# ---------------------------------------------------------------------------
# Behavioural unit tests for sem_mutations_get() — A-10
#
# These tests verify:
#   1. The Polars inner join correctly computes mutations over the
#      intersection of input positions and beta-range positions.
#   2. Extra positions in values (not in thresholds) are silently dropped.
#   3. Extra positions in thresholds (not in values) are silently dropped.
#   4. Zero overlap returns an empty data.frame without crashing.
#   5. Both HYPO and HYPER directions work.
#
# All tests use synthetic, fully controlled data so expected outputs can be
# verified analytically without loading real methylation data.
# ---------------------------------------------------------------------------

# Helper: minimal values data.frame (CHR, START, END, VALUE)
.make_values_df_mut <- function(chr, starts, values_vec) {
  data.frame(
    CHR   = as.character(chr),
    START = as.integer(starts),
    END   = as.integer(starts),
    VALUE = as.numeric(values_vec),
    stringsAsFactors = FALSE
  )
}

# Helper: minimal thresholds data.frame
.make_thresholds_df_mut <- function(chr, starts, inf_thresh, sup_thresh) {
  data.frame(
    CHR                        = as.character(chr),
    START                      = as.integer(starts),
    END                        = as.integer(starts),
    signal_inferior_thresholds = as.numeric(inf_thresh),
    signal_superior_thresholds = as.numeric(sup_thresh),
    stringsAsFactors           = FALSE
  )
}

test_that("sem_mutations_get (A-10): HYPO full overlap — correct mutation count", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # 10 positions; threshold: inferior = 0.3, superior = 0.7
  # Values 1-5: 0.1 (below 0.3 → HYPO mutation)
  # Values 6-10: 0.5 (above inferior, below superior → no mutation)
  starts    <- seq(1000L, by = 1000L, length.out = 10L)
  values_df <- .make_values_df_mut("chr1", starts, c(rep(0.1, 5L), rep(0.5, 5L)))
  thresh_df <- .make_thresholds_df_mut("chr1", starts,
                                       inf_thresh = rep(0.3, 10L),
                                       sup_thresh = rep(0.7, 10L))

  result <- SEMseeker:::sem_mutations_get(
    values     = values_df,
    figure     = "HYPO",
    thresholds = thresh_df,
    sampleName = "test_hypo_full"
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("CHR", "START", "END", "MUTATIONS") %in% colnames(result)))
  expect_equal(nrow(result), 10L)          # all 10 positions covered
  expect_equal(sum(result$MUTATIONS), 5L)  # exactly 5 below-threshold positions
})

test_that("sem_mutations_get (A-10): HYPER full overlap — correct mutation count", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  starts    <- seq(1000L, by = 1000L, length.out = 10L)
  # Values 6-10: 0.9 (above 0.7 → HYPER mutation)
  values_df <- .make_values_df_mut("chr1", starts, c(rep(0.5, 5L), rep(0.9, 5L)))
  thresh_df <- .make_thresholds_df_mut("chr1", starts,
                                       inf_thresh = rep(0.3, 10L),
                                       sup_thresh = rep(0.7, 10L))

  result <- SEMseeker:::sem_mutations_get(
    values     = values_df,
    figure     = "HYPER",
    thresholds = thresh_df,
    sampleName = "test_hyper_full"
  )

  expect_equal(nrow(result), 10L)
  expect_equal(sum(result$MUTATIONS), 5L)
})

test_that("sem_mutations_get (A-10): values has extra positions → only intersection returned", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # Cross-run scenario: 15 value positions, 10 threshold positions
  # Only positions 1-10 are shared; positions 11-15 exist only in values.
  all_starts    <- seq(1000L, by = 1000L, length.out = 15L)
  shared_starts <- all_starts[1:10]

  values_df <- .make_values_df_mut("chr1", all_starts, rep(0.1, 15L))
  thresh_df <- .make_thresholds_df_mut("chr1", shared_starts,
                                       inf_thresh = rep(0.3, 10L),
                                       sup_thresh = rep(0.7, 10L))

  result <- SEMseeker:::sem_mutations_get(
    values     = values_df,
    figure     = "HYPO",
    thresholds = thresh_df,
    sampleName = "test_extra_values"
  )

  # Only the 10 shared positions should appear
  expect_equal(nrow(result), 10L)
  expect_equal(sum(result$MUTATIONS), 10L)
  # Positions 11-15 (starts 11000+) must NOT appear in output
  expect_true(all(result$START %in% shared_starts))
})

test_that("sem_mutations_get (A-10): thresholds has extra positions → only intersection returned", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # 10 value positions, 15 threshold positions — overlap = positions 1-10
  all_starts    <- seq(1000L, by = 1000L, length.out = 15L)
  shared_starts <- all_starts[1:10]

  values_df <- .make_values_df_mut("chr1", shared_starts, rep(0.9, 10L))
  thresh_df <- .make_thresholds_df_mut("chr1", all_starts,
                                       inf_thresh = rep(0.3, 15L),
                                       sup_thresh = rep(0.7, 15L))

  result <- SEMseeker:::sem_mutations_get(
    values     = values_df,
    figure     = "HYPER",
    thresholds = thresh_df,
    sampleName = "test_extra_thresholds"
  )

  expect_equal(nrow(result), 10L)
  expect_equal(sum(result$MUTATIONS), 10L)
  expect_true(all(result$START %in% shared_starts))
})

test_that("sem_mutations_get (A-10): zero overlap → empty data.frame, no crash", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # Values on chr1, thresholds on chr2 — no shared positions
  values_df <- .make_values_df_mut("chr1",
                                   seq(1000L, by = 1000L, length.out = 10L),
                                   rep(0.1, 10L))
  thresh_df <- .make_thresholds_df_mut("chr2",
                                       seq(1000L, by = 1000L, length.out = 10L),
                                       inf_thresh = rep(0.3, 10L),
                                       sup_thresh = rep(0.7, 10L))

  result <- SEMseeker:::sem_mutations_get(
    values     = values_df,
    figure     = "HYPO",
    thresholds = thresh_df,
    sampleName = "test_zero_overlap"
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0L)
  expect_true(all(c("CHR", "START", "END", "MUTATIONS") %in% colnames(result)))
})

# ---------------------------------------------------------------------------
# Original integration test (preserved unchanged)
# ---------------------------------------------------------------------------

test_that("sem_mutations_get", {

  tempFolder <- sem_test_folder()
  ssEnv <- SEMseeker:::core_init_env(result_folder = tempFolder, inpute = "median")

  if (!exists("signal_thresholds")) {
    signal_data <- SEMseeker:::sem_inpute_missing_values(signal_data)
    signal_thresholds <<- SEMseeker:::sem_signal_range_values(signal_data, batch_id)
  }
  probe_features <<- probe_features[probe_features$PROBE %in% rownames(signal_data), ]

  values_df <- data.frame(
    CHR   = probe_features$CHR[match(rownames(signal_data), probe_features$PROBE)],
    START = probe_features$START[match(rownames(signal_data), probe_features$PROBE)],
    END   = probe_features$END[match(rownames(signal_data), probe_features$PROBE)],
    VALUE = as.numeric(signal_data[, 1])
  )

  # ── HYPO: basic existence ──────────────────────────────────────────────────
  mutations_hypo <- SEMseeker:::sem_mutations_get(
    values     = values_df,
    figure     = "HYPO",
    thresholds = signal_thresholds,
    sampleName = mySampleSheet[1, "Sample_ID"]
  )

  testthat::expect_false(length(mutations_hypo) == 0)

  # ── HYPO: output structure ─────────────────────────────────────────────────
  testthat::expect_s3_class(mutations_hypo, "data.frame")
  testthat::expect_true(all(c("CHR", "START", "END", "MUTATIONS") %in% colnames(mutations_hypo)))

  # MUTATIONS column must be binary (0 / 1)
  testthat::expect_true(all(mutations_hypo$MUTATIONS %in% c(0, 1)))

  # row count matches sorted thresholds (no probes lost)
  testthat::expect_true(nrow(mutations_hypo) > 0)

  # ── HYPER: symmetric test ──────────────────────────────────────────────────
  mutations_hyper <- SEMseeker:::sem_mutations_get(
    values     = values_df,
    figure     = "HYPER",
    thresholds = signal_thresholds,
    sampleName = mySampleSheet[1, "Sample_ID"]
  )

  testthat::expect_s3_class(mutations_hyper, "data.frame")
  testthat::expect_true(all(c("CHR", "START", "END", "MUTATIONS") %in% colnames(mutations_hyper)))
  testthat::expect_true(all(mutations_hyper$MUTATIONS %in% c(0, 1)))
  # row count must be identical for HYPO and HYPER (same probes, different direction)
  testthat::expect_equal(nrow(mutations_hypo), nrow(mutations_hyper))

  # ── Boundary: threshold = -Inf  →  zero HYPO mutations ────────────────────
  thresholds_zero <- signal_thresholds
  thresholds_zero$signal_inferior_thresholds <- -Inf
  mutations_none <- SEMseeker:::sem_mutations_get(
    values     = values_df,
    figure     = "HYPO",
    thresholds = thresholds_zero,
    sampleName = mySampleSheet[1, "Sample_ID"]
  )
  testthat::expect_equal(sum(mutations_none$MUTATIONS), 0)

  # ── Boundary: threshold = +Inf  →  all HYPO mutations ─────────────────────
  thresholds_all <- signal_thresholds
  thresholds_all$signal_inferior_thresholds <- Inf
  mutations_all <- SEMseeker:::sem_mutations_get(
    values     = values_df,
    figure     = "HYPO",
    thresholds = thresholds_all,
    sampleName = mySampleSheet[1, "Sample_ID"]
  )
  testthat::expect_equal(sum(mutations_all$MUTATIONS), nrow(mutations_all))

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})
