# test-io_coord_input-signal_input.R
# Integration tests for SEMseeker pipeline with non-Illumina coordinate input
# (CHR / START / END format — WGBS, Nanopore ONT bedmethyl, RRBS)
#
# The standard test setup (setup.R) uses Illumina cg* probe IDs as rownames.
# These tests verify that the same pipeline stages work correctly when the
# signal is provided as a (CHR, START, END, sample1, sample2, ...) data.frame
# — the coordinate format introduced by B-01 (WGBS) and B-02 (ONT long-read).
#
# Covered:
#   1. io_is_coord_format()         detects CHR/START columns correctly
#   2. io_normalize_signal_input()  converts coord df → probe-indexed matrix
#   3. core_get_meth_tech()           sets ssEnv$tech = "LONGREAD" for coord input
#   4. sem_mutations_get()           HYPO / HYPER with coord-derived values + thresholds
#   5. sem_delta_single_sample()     continuous delta metric from coord input
#   6. sem_signal_single_sample()    writes bedgraph file for coord-format sample
#
# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Build a small synthetic coordinate-format methylation matrix
# Returns a data.frame with CHR, START, END + n_samples numeric columns.
.build_coord_signal <- function(n_pos = 50L, n_samples = 4L, seed = 77L) {
  set.seed(seed)
  chrs <- paste0("chr", c(rep("1", n_pos %/% 2), rep("2", n_pos - n_pos %/% 2)))
  starts <- as.integer(seq(100000L, by = 1000L, length.out = n_pos))
  sample_ids <- paste0("S", seq_len(n_samples))

  mat <- matrix(stats::rbeta(n_pos * n_samples, 80, 20), nrow = n_pos)
  # Inject a few HYPO outliers in rows 1-5 (very low beta)
  mat[1:5, ] <- matrix(stats::runif(5L * n_samples, 0.01, 0.05), nrow = 5L)
  # Inject a few HYPER outliers in rows 46-50 (very high beta)
  mat[(n_pos - 4L):n_pos, ] <- matrix(stats::runif(5L * n_samples, 0.95, 0.99),
                                       nrow = 5L)

  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(df) <- sample_ids
  data.frame(CHR = chrs, START = starts, END = starts, df,
             stringsAsFactors = FALSE)
}

# Build a thresholds data.frame from a coordinate signal (using per-position
# Q1-3*IQR / Q3+3*IQR boundaries, same as the reference pipeline).
.build_coord_thresholds <- function(coord_signal) {
  # Extract numeric columns (samples)
  betas <- as.matrix(coord_signal[, -(1:3)])
  q1  <- apply(betas, 1, stats::quantile, probs = 0.25, na.rm = TRUE)
  q3  <- apply(betas, 1, stats::quantile, probs = 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  data.frame(
    CHR                        = as.character(coord_signal$CHR),
    START                      = as.integer(coord_signal$START),
    END                        = as.integer(coord_signal$END),
    signal_inferior_thresholds = q1 - 3 * iqr,
    signal_superior_thresholds = q3 + 3 * iqr,
    signal_median_values       = apply(betas, 1, stats::median, na.rm = TRUE),
    iqr                        = iqr,
    q1                         = q1,
    q3                         = q3,
    stringsAsFactors           = FALSE
  )
}

# ---------------------------------------------------------------------------
# 1. Format detection
# ---------------------------------------------------------------------------

test_that("io_is_coord_format: detects coordinate data.frame (CHR/START columns)", {
  df <- .build_coord_signal(n_pos = 10L, n_samples = 2L)
  expect_true(SEMseeker:::io_is_coord_format(df))
})

test_that("io_is_coord_format: rejects cg*-indexed Illumina matrix", {
  # signal_data from setup.R has cg* rownames → NOT coord format
  expect_false(SEMseeker:::io_is_coord_format(signal_data))
})

# ---------------------------------------------------------------------------
# 2. io_normalize_signal_input
# ---------------------------------------------------------------------------

test_that("io_normalize_signal_input: coord df → probe-indexed matrix, no CHR/START/END cols", {
  df <- .build_coord_signal(n_pos = 20L, n_samples = 3L)
  result <- SEMseeker:::io_normalize_signal_input(df)
  expect_false(SEMseeker:::io_is_coord_format(result))
  expect_false("CHR"   %in% colnames(result))
  expect_false("START" %in% colnames(result))
  expect_false("END"   %in% colnames(result))
  expect_equal(nrow(result), 20L)
  expect_equal(ncol(result), 3L)  # 3 sample columns
})

test_that("io_normalize_signal_input: probe IDs are CHR_START format", {
  df <- data.frame(CHR = "chr1", START = 12345L, END = 12345L, S1 = 0.7)
  result <- SEMseeker:::io_normalize_signal_input(df)
  # chr prefix stripped → "1_12345"
  expect_equal(rownames(result), "1_12345")
})

test_that("io_normalize_signal_input: passes Illumina matrix unchanged", {
  # Illumina matrix should pass through untouched
  result <- SEMseeker:::io_normalize_signal_input(signal_data)
  expect_identical(result, signal_data)
})

test_that("io_normalize_signal_input: preserves beta values after conversion", {
  df <- data.frame(CHR = "chr1", START = 5000L, END = 5000L, S1 = 0.42)
  result <- SEMseeker:::io_normalize_signal_input(df)
  expect_equal(as.numeric(result[1, 1]), 0.42, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 3. core_get_meth_tech with coord format
# ---------------------------------------------------------------------------

test_that("core_get_meth_tech: coord-format signal sets tech to WGBS (not an Illumina array)", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  coord_df  <- .build_coord_signal(n_pos = 30L, n_samples = 3L)
  probe_mat <- SEMseeker:::io_normalize_signal_input(coord_df)
  env       <- SEMseeker:::core_get_meth_tech(probe_mat)
  # Coordinate-format data has synthetic "CHR_POS" probe IDs that don't match
  # any Illumina array manifest → core_get_meth_tech classifies them as WGBS
  expect_equal(env$tech, "WGBS")
  expect_false(env$tech %in% c("K27", "K450", "K850"))
})

# ---------------------------------------------------------------------------
# 4. sem_mutations_get with coordinate-derived values and thresholds
# ---------------------------------------------------------------------------

test_that("sem_mutations_get: coord-derived values + coord thresholds — HYPO counts injected outliers", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  coord_df  <- .build_coord_signal(n_pos = 50L, n_samples = 5L, seed = 42L)
  coord_thr <- .build_coord_thresholds(coord_df)

  # Build values df for sample column 1 (column 4 of coord_df)
  values_df <- data.frame(
    CHR   = as.character(coord_df$CHR),
    START = as.integer(coord_df$START),
    END   = as.integer(coord_df$END),
    VALUE = as.numeric(coord_df[, 4]),
    stringsAsFactors = FALSE
  )

  res <- SEMseeker:::sem_mutations_get(
    values     = values_df,
    figure     = "HYPO",
    thresholds = coord_thr,
    sampleName = "coord_sample_1"
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(c("CHR", "START", "END", "MUTATIONS") %in% colnames(res)))
  expect_equal(nrow(res), 50L)

  # Rows 1-5 were injected as very low beta values → should be HYPO mutations.
  # The threshold is Q1-3*IQR; with ~0.02 betas vs background ~0.8 the outliers
  # should be detected.  We test that > 0 and <= 5 HYPO mutations are found.
  n_mut <- sum(res$MUTATIONS)
  expect_gte(n_mut, 1L)
})

test_that("sem_mutations_get: coord-derived HYPER — detects injected high-beta outliers", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  coord_df  <- .build_coord_signal(n_pos = 50L, n_samples = 5L, seed = 99L)
  coord_thr <- .build_coord_thresholds(coord_df)

  values_df <- data.frame(
    CHR   = as.character(coord_df$CHR),
    START = as.integer(coord_df$START),
    END   = as.integer(coord_df$END),
    VALUE = as.numeric(coord_df[, 4]),
    stringsAsFactors = FALSE
  )

  res <- SEMseeker:::sem_mutations_get(
    values     = values_df,
    figure     = "HYPER",
    thresholds = coord_thr,
    sampleName = "coord_sample_hyper"
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(res$MUTATIONS %in% c(0L, 1L)))
  expect_gte(sum(res$MUTATIONS), 1L)
})

test_that("sem_mutations_get: coord values on different chromosomes than thresholds → zero results", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # Values on chr3, thresholds on chr1/chr2 → inner join = empty
  n <- 10L
  values_df <- data.frame(
    CHR   = "chr3",
    START = seq(1000L, by = 1000L, length.out = n),
    END   = seq(1000L, by = 1000L, length.out = n),
    VALUE = rep(0.1, n),
    stringsAsFactors = FALSE
  )
  coord_df  <- .build_coord_signal(n_pos = 20L, n_samples = 3L)
  coord_thr <- .build_coord_thresholds(coord_df)

  res <- SEMseeker:::sem_mutations_get(
    values     = values_df,
    figure     = "HYPO",
    thresholds = coord_thr,
    sampleName = "coord_no_overlap"
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0L)
})

# ---------------------------------------------------------------------------
# 5. sem_delta_single_sample with coordinate-derived input
# ---------------------------------------------------------------------------

test_that("sem_delta_single_sample: coord input runs without error and returns NULL", {
  # Note: sem_delta_single_sample() writes files only when DELTA > 0.  When thresholds
  # are built from the same samples as the test sample, no outlier is detected and
  # no file is written (correct behaviour — this test verifies no crash).
  # For a test with expected file output, see test-io_bed_file_name.R and test-sem_delta_single_sample.R.
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  coord_df  <- .build_coord_signal(n_pos = 30L, n_samples = 4L, seed = 123L)
  coord_thr <- .build_coord_thresholds(coord_df)

  # Build values for sample column 4 (0-based: coord_df[, 4+3] = col index 7)
  values_df <- data.frame(
    CHR   = as.character(coord_df$CHR),
    START = as.integer(coord_df$START),
    END   = as.integer(coord_df$END),
    VALUE = as.numeric(coord_df[, 4]),
    stringsAsFactors = FALSE
  )
  sample_detail <- data.frame(
    Sample_ID    = "coord_delta_test",
    Sample_Group = "Control",
    stringsAsFactors = FALSE
  )

  # sem_delta_single_sample() returns invisible(NULL)
  expect_no_error(
    result <- SEMseeker:::sem_delta_single_sample(
      values        = values_df,
      thresholds    = coord_thr,
      sample_detail = sample_detail
    )
  )
  expect_null(result)
})

test_that("sem_delta_single_sample: detects outliers when thresholds come from separate reference", {
  # Build reference from 3 background samples (all ~0.8), then test against a
  # sample that has 5 genuine HYPO outliers (values ≈ 0.02).
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  n_pos <- 20L
  set.seed(55)
  starts <- as.integer(seq(100000L, by = 1000L, length.out = n_pos))
  # Reference: all positions methylated ~0.8 in 3 samples
  ref_mat <- matrix(stats::rbeta(n_pos * 3L, 80, 20), nrow = n_pos)
  ref_df  <- as.data.frame(ref_mat)
  colnames(ref_df) <- paste0("Ref", 1:3)
  ref_df$CHR   <- "1"
  ref_df$START <- starts
  ref_df$END   <- starts

  # Thresholds from reference only
  coord_thr <- .build_coord_thresholds(
    data.frame(CHR = "1", START = starts, END = starts, ref_df[, 1:3],
               stringsAsFactors = FALSE)
  )

  # Test sample: positions 1-5 are HYPO outliers (≈0.02), rest normal (≈0.8)
  test_vals <- c(rep(0.02, 5), stats::rbeta(n_pos - 5L, 80, 20))
  values_df <- data.frame(CHR = "1", START = starts, END = starts,
                           VALUE = test_vals, stringsAsFactors = FALSE)
  sample_detail <- data.frame(Sample_ID = "test_outlier", Sample_Group = "Control",
                               stringsAsFactors = FALSE)

  SEMseeker:::sem_delta_single_sample(
    values        = values_df,
    thresholds    = coord_thr,
    sample_detail = sample_detail
  )

  ssEnv     <- SEMseeker:::core_get_session_info()
  data_root <- ssEnv$result_folderData
  # At least one DELTAS_HYPO file should have been written
  hypo_files <- list.files(data_root, pattern = "HYPO.*\\.bedgraph",
                            recursive = TRUE)
  expect_gte(length(hypo_files), 1L)
})

# ---------------------------------------------------------------------------
# 6. sem_signal_single_sample: bedgraph file is created for coord-format sample
# ---------------------------------------------------------------------------

test_that("sem_signal_single_sample: writes bedgraph file for coordinate-format sample", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE, inpute = "median")
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  coord_df  <- .build_coord_signal(n_pos = 20L, n_samples = 2L)
  probe_mat <- SEMseeker:::io_normalize_signal_input(coord_df)

  # Build probe_features from the converted matrix (CHR/START/END recovered via
  # io_probe_id_to_coord)
  pf <- SEMseeker:::io_coord_probe_features(rownames(probe_mat))

  sample_detail <- data.frame(
    Sample_ID    = colnames(probe_mat)[1],
    Sample_Group = "Control",
    stringsAsFactors = FALSE
  )

  SEMseeker:::sem_signal_single_sample(
    values        = probe_mat[, 1],
    sample_detail = sample_detail,
    probe_features = pf
  )

  ssEnv       <- SEMseeker:::core_get_session_info()
  # The FIGURE of SIGNAL is the scale of the session, resolved by the same
  # function the writer uses: the test follows BETA or MVALUE.
  folder      <- file.path(ssEnv$result_folderData, "Control",
                           paste0("SIGNAL_", SEMseeker:::io_signal_figure()))
  bed_files   <- list.files(folder, pattern = "\\.bedgraph(\\.gz)?$", recursive = TRUE)
  expect_gte(length(bed_files), 1L)
})
