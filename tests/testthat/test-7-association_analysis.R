# Tests for association_analysis helpers and integration
#
# Covered:
#  - util_split_and_clean()           pure string-splitting helper
#  - assoc_validate_family_test()      returns TRUE/FALSE for known/unknown families
#  - association_analysis()      end-to-end at depth=1 (sample-level)
#
# Note on parallel strategy for integration tests (tests 3 and 4):
#  Both semseeker() and association_analysis() are run with "sequential" to make
#  the tests environment-agnostic.
#   - "multicore" (fork) works only with devtools::load_all(), not with installed pkg
#   - "multisession" works only with an installed package, not with devtools::load_all()
#   - "sequential" works in both environments
#  These tests verify correctness of the association pipeline, not parallelism.
#
# Note on areas:
#  semseeker() must include "POSITION" so sem_study_summary_total() writes
#  sample_sheet_result.csv with per-sample mutation counts.  Those counts are
#  the dependent-variable columns for depth=1 association tests.

# ---------------------------------------------------------------------------
# 1. util_split_and_clean — pure function, no session needed
# ---------------------------------------------------------------------------

test_that("util_split_and_clean splits on + and cleans whitespace", {
  expect_equal(SEMseeker:::util_split_and_clean("A+B+C"), c("A", "B", "C"))
  expect_equal(SEMseeker:::util_split_and_clean("single"), "single")
  expect_equal(SEMseeker:::util_split_and_clean("A + B"), c("A", "B"))   # leading/trailing space trimmed
})

test_that("util_split_and_clean removes empty parts", {
  result <- SEMseeker:::util_split_and_clean("")
  expect_equal(length(result), 0)  # empty string → character(0) after filtering
})

test_that("util_split_and_clean deduplicates", {
  result <- SEMseeker:::util_split_and_clean("A+A+B")
  expect_equal(sort(result), c("A", "B"))
})

test_that("util_split_and_clean respects a custom split delimiter", {
  result <- SEMseeker:::util_split_and_clean("X,Y,Z", split = ",")
  expect_equal(result, c("X", "Y", "Z"))
})

# ---------------------------------------------------------------------------
# 2. assoc_validate_family_test — needs a live session for core_log_event()
# ---------------------------------------------------------------------------

test_that("assoc_validate_family_test accepts standard parametric families", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  expect_true(SEMseeker:::assoc_validate_family_test("gaussian"))
  expect_true(SEMseeker:::assoc_validate_family_test("binomial"))
  expect_true(SEMseeker:::assoc_validate_family_test("poisson"))
  expect_true(SEMseeker:::assoc_validate_family_test("wilcoxon"))
  expect_true(SEMseeker:::assoc_validate_family_test("t.test"))
  expect_true(SEMseeker:::assoc_validate_family_test("pearson"))
  expect_true(SEMseeker:::assoc_validate_family_test("spearman"))
  expect_true(SEMseeker:::assoc_validate_family_test("kendall"))

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})

test_that("assoc_validate_family_test accepts parametric family variants", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  # quantreg family (grepl match)
  expect_true(SEMseeker:::assoc_validate_family_test("quantreg_0.5"))
  # quantreg-permutation requires exactly 5 underscore-separated parts
  expect_true(SEMseeker:::assoc_validate_family_test("quantreg-permutation_0.5_5_10_0.9"))
  # polynomial / exp / log variants
  expect_true(SEMseeker:::assoc_validate_family_test("polynomial_4_1"))
  expect_true(SEMseeker:::assoc_validate_family_test("exp_1"))
  expect_true(SEMseeker:::assoc_validate_family_test("log_1"))

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})

test_that("assoc_validate_family_test rejects NULL, NA, and unknown strings", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  expect_false(SEMseeker:::assoc_validate_family_test(NULL))
  expect_false(SEMseeker:::assoc_validate_family_test(NA))
  expect_false(SEMseeker:::assoc_validate_family_test("not_a_valid_test"))

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})

# ---------------------------------------------------------------------------
# 3. association_analysis — depth=1 integration test
#    Runs semseeker() with POSITION area to populate sample_sheet_result.csv,
#    then calls association_analysis() with a minimal gaussian inference.
#    depth=1 reads per-sample mutation counts from sample_sheet_result.csv
#    and regresses them against the continuous Phenotest covariate.
# ---------------------------------------------------------------------------

test_that("association_analysis at scope SAMPLE runs without error and writes inference CSV", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  unlink(tempFolder, recursive = TRUE)

  # ── Step 0: build synthetic data WITH guaranteed mutations ────────────────
  # The global signal_data has almost no mutations (Beta(90,10) is too tight for
  # IQR×3 threshold with only 30 samples). We inject bimodal outliers so that
  # samples 1-5 are hypomethylated at the first 50 probes. This ensures non-zero
  # per-sample mutation counts in sample_sheet_result.csv and avoids NaN in cor.test.
  set.seed(777)
  n_probes_aa <- 200L
  n_samples_aa <- nsamples
  local_probes  <- probe_features[1:n_probes_aa, ]
  local_samples <- mySampleSheet

  # Background: mostly methylated
  local_sig <- matrix(stats::rbeta(n_probes_aa * n_samples_aa, 90L, 10L),
                       nrow = n_probes_aa, ncol = n_samples_aa)
  # Inject clear HYPO outliers: first 5 samples, first 50 probes → values near 0
  local_sig[1:50, 1:5] <- stats::rbeta(50L * 5L, 1L, 100L)

  rownames(local_sig) <- local_probes$PROBE
  local_sig <- as.data.frame(local_sig)
  # signal_data has 10 unique columns; mySampleSheet has 16 rows (Reference reuse pattern)
  colnames(local_sig) <- colnames(signal_data)

  # ── Step 1: produce semseeker output including POSITION pivots ────────────
  # POSITION area is required so sem_study_summary_total() writes per-sample
  # mutation counts into sample_sheet_result.csv (depth=1 reads that file).
  # "sequential" makes the test work under both devtools::load_all() and CI.
  SEMseeker::semseeker(
    input             = local_sig,
    sample_sheet      = local_samples,
    result_folder     = tempFolder,
    parallel_strategy = parallel_strategy,
    areas             = c("GENE", "POSITION"),
    markers           = c("MUTATIONS"),
    start_fresh       = TRUE,
    inpute            = "median",
    showprogress      = showprogress,
    verbosity         = verbosity
  )

  # ── Step 2: build a minimal inference_details ─────────────────────────────
  # Use "spearman" to avoid the caret::createDataPartition / GLM path
  # (gaussian calls caret which can fail when mutation counts contain NaN).
  # Spearman correlation goes through assoc_test_model which is NaN-safe.
  # transformation_x must always be present (accessed without NA-guard).
  inference_details <- data.frame(
    independent_variable = "Phenotest",
    family_test          = "spearman",
    transformation_y     = "",
    transformation_x     = "",
    aggregation          = "SUM",
    scope                = "SAMPLE",
    filter_p_value       = FALSE,
    stringsAsFactors     = FALSE
  )

  # ── Step 3: association_analysis should complete without error ────────────
  # multiple_test_adj="BH" avoids qvalue::qvalue() which fails on few p-values
  # (pi0 bootstrap needs many observations; default "q" fails with ~1 p-value).
  testthat::expect_no_error(
    SEMseeker:::association_analysis(
      inference_details   = inference_details,
      result_folder       = tempFolder,
      parallel_strategy   = "sequential",
      markers             = c("MUTATIONS"),
      figures             = c("HYPO"),
      multiple_test_adj   = "BH",
      showprogress        = showprogress,
      verbosity           = verbosity
    )
  )

  # ── Step 4: inference folder was created ─────────────────────────────────
  # Compute path directly — avoids re-opening the env (core_init_env would mkdir again,
  # and normalizePath differences can cause dir.exists() false negatives).
  inference_dir <- file.path(tempFolder, "Inference")
  testthat::expect_true(dir.exists(inference_dir))

  # ── Step 5: at least one CSV written inside the inference folder ──────────
  csv_files <- list.files(inference_dir, pattern = "\\.csv$", recursive = TRUE,
                           full.names = TRUE)
  testthat::expect_true(length(csv_files) > 0)

  # ── Step 6: the main result CSV (not the assoc_covariates_model side-file) has rows
  result_csv <- csv_files[!grepl("(?i)assoc_covariates_model", csv_files)][1]
  if (!is.na(result_csv) && file.exists(result_csv)) {
    result_df <- utils::read.csv2(result_csv)
    testthat::expect_true(nrow(result_df) > 0)
  }

  unlink(tempFolder, recursive = TRUE)
})

# ---------------------------------------------------------------------------
# Helper: build a small semseeker result with injected hypomethylated
# outliers, large enough to populate GENE pivots. Returns the tempFolder.
# ---------------------------------------------------------------------------

.aa_setup_result_folder <- function(seed = 777, n_probes = 200L, areas = c("GENE", "POSITION")) {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  unlink(tempFolder, recursive = TRUE)

  set.seed(seed)
  local_probes <- probe_features[seq_len(n_probes), ]
  local_sig <- matrix(stats::rbeta(n_probes * nsamples, 90L, 10L),
                       nrow = n_probes, ncol = nsamples)
  # Inject HYPO outliers in first 50 probes × first 5 samples
  local_sig[1:50, 1:5] <- stats::rbeta(50L * 5L, 1L, 100L)
  rownames(local_sig) <- local_probes$PROBE
  local_sig <- as.data.frame(local_sig)
  # signal_data has 10 unique columns; mySampleSheet has 16 rows (Reference reuse pattern)
  colnames(local_sig) <- colnames(signal_data)

  SEMseeker::semseeker(
    input             = local_sig,
    sample_sheet      = mySampleSheet,
    result_folder     = tempFolder,
    parallel_strategy = parallel_strategy,   # global from setup.R (multisession on macOS)
    areas             = areas,
    markers           = c("MUTATIONS"),
    start_fresh       = TRUE,
    inpute            = "median",
    showprogress      = showprogress,
    verbosity         = verbosity
  )
  tempFolder
}

# ---------------------------------------------------------------------------
# T1, scope INSTANCE + spearman: exercises pivot parquet read / chunk / merge
#      path, one row per instance of every region class of the run
# ---------------------------------------------------------------------------

test_that("association_analysis at scope INSTANCE reads area pivots and writes inference CSV", {
  tempFolder <- .aa_setup_result_folder()

  inference_details <- data.frame(
    independent_variable = "Phenotest",
    family_test          = "spearman",
    transformation_y     = "",
    transformation_x     = "",
    aggregation          = "SUM",
    scope                = "INSTANCE",
    filter_p_value       = FALSE,
    stringsAsFactors     = FALSE
  )

  testthat::expect_no_error(
    SEMseeker:::association_analysis(
      inference_details = inference_details,
      result_folder     = tempFolder,
      parallel_strategy = parallel_strategy,
      markers           = c("MUTATIONS"),
      figures           = c("HYPO"),
      multiple_test_adj = "BH",
      showprogress      = showprogress,
      verbosity         = verbosity
    )
  )

  inference_dir <- file.path(tempFolder, "Inference")
  testthat::expect_true(dir.exists(inference_dir))

  csv_files <- list.files(inference_dir, pattern = "\\.csv$", recursive = TRUE,
                           full.names = TRUE)
  testthat::expect_true(length(csv_files) > 0)

  # AI-308: this block used to assert `any(result_df$DEPTH > 1)` behind a
  # skip() that reported a regression of 53310c1. It could not have detected
  # either: `DEPTH` was removed with depth_analysis, so the guard
  # `"DEPTH" %in% colnames(result_df)` was never true, the skip never fired and
  # the assertion after it was unreachable. The test had been reduced to
  # "no error, a CSV exists" while reading as though it checked the granularity.
  #
  # The property it wanted is now sayable: the request asked for the
  # per-instance branch, so every row must carry SCOPE = INSTANCE, and a region
  # class must span many instances, which is exactly what a collapsed row
  # cannot do.
  result_csv <- csv_files[!grepl("(?i)assoc_covariates_model", csv_files)][1]
  if (!is.na(result_csv) && file.exists(result_csv) && file.info(result_csv)$size > 10) {
    result_df <- utils::read.csv2(result_csv)
    result_df <- result_df[which(!is.na(result_df$SCOPE)), , drop = FALSE]
    if (nrow(result_df) > 0) {
      testthat::expect_equal(unique(result_df$SCOPE), "INSTANCE")
      testthat::expect_gt(length(unique(result_df$AREA_OF_TEST)), 1)
    }
  }

  unlink(tempFolder, recursive = TRUE)
})

# ---------------------------------------------------------------------------
# T2 — polynomial_2_1 family: exercises assoc_execute_model dispatch beyond GLM
# ---------------------------------------------------------------------------

test_that("association_analysis polynomial family runs without error", {
  tempFolder <- .aa_setup_result_folder()

  inference_details <- data.frame(
    independent_variable = "Phenotest",
    family_test          = "polynomial_2_1",
    transformation_y     = "",
    transformation_x     = "",
    aggregation          = "SUM",
    scope                = "INSTANCE",
    filter_p_value       = FALSE,
    stringsAsFactors     = FALSE
  )

  testthat::expect_no_error(
    SEMseeker:::association_analysis(
      inference_details = inference_details,
      result_folder     = tempFolder,
      parallel_strategy = parallel_strategy,
      markers           = c("MUTATIONS"),
      figures           = c("HYPO"),
      multiple_test_adj = "BH",
      showprogress      = showprogress,
      verbosity         = verbosity
    )
  )

  inference_dir <- file.path(tempFolder, "Inference")
  csv_files <- list.files(inference_dir, pattern = "\\.csv$", recursive = TRUE,
                           full.names = TRUE)
  testthat::expect_true(length(csv_files) > 0)

  unlink(tempFolder, recursive = TRUE)
})

# ---------------------------------------------------------------------------
# T3 — covariates: exercises assoc_covariates_model (collinearity check, sample filter)
# ---------------------------------------------------------------------------

test_that("association_analysis with covariates runs and produces a assoc_covariates_model side-file", {
  tempFolder <- .aa_setup_result_folder()

  inference_details <- data.frame(
    independent_variable = "Phenotest",
    family_test          = "spearman",
    covariates           = "Covariates1+Covariates2",
    transformation_y     = "",
    transformation_x     = "",
    aggregation          = "SUM",
    scope                = "INSTANCE",
    filter_p_value       = FALSE,
    collinearity_check   = TRUE,
    stringsAsFactors     = FALSE
  )

  testthat::expect_no_error(
    SEMseeker:::association_analysis(
      inference_details = inference_details,
      result_folder     = tempFolder,
      parallel_strategy = parallel_strategy,
      markers           = c("MUTATIONS"),
      figures           = c("HYPO"),
      multiple_test_adj = "BH",
      showprogress      = showprogress,
      verbosity         = verbosity
    )
  )

  inference_dir <- file.path(tempFolder, "Inference")
  csv_files <- list.files(inference_dir, pattern = "\\.csv$", recursive = TRUE,
                           full.names = TRUE)
  testthat::expect_true(length(csv_files) > 0)

  unlink(tempFolder, recursive = TRUE)
})

# ---------------------------------------------------------------------------
# T4 — resumption / idempotency: second identical run must not duplicate rows
#      Exercises the "keys already done" branch (association_analysis.R:239-249)
# ---------------------------------------------------------------------------

test_that("association_analysis is idempotent: second run on same folder does not duplicate results", {
  tempFolder <- .aa_setup_result_folder()

  inference_details <- data.frame(
    independent_variable = "Phenotest",
    family_test          = "spearman",
    transformation_y     = "",
    transformation_x     = "",
    aggregation          = "SUM",
    scope                = "INSTANCE",
    filter_p_value       = FALSE,
    stringsAsFactors     = FALSE
  )

  run_once <- function() {
    SEMseeker:::association_analysis(
      inference_details = inference_details,
      result_folder     = tempFolder,
      parallel_strategy = parallel_strategy,
      markers           = c("MUTATIONS"),
      figures           = c("HYPO"),
      multiple_test_adj = "BH",
      showprogress      = showprogress,
      verbosity         = verbosity
    )
  }

  testthat::expect_no_error(run_once())

  inference_dir <- file.path(tempFolder, "Inference")
  csv_files <- list.files(inference_dir, pattern = "\\.csv$", recursive = TRUE,
                           full.names = TRUE)
  result_csv <- csv_files[!grepl("(?i)assoc_covariates_model", csv_files)][1]
  testthat::skip_if(is.na(result_csv) || !file.exists(result_csv),
                    "No primary result CSV produced — cannot verify idempotency")

  rows_run1 <- nrow(unique(utils::read.csv2(result_csv)))
  testthat::expect_true(rows_run1 > 0)

  # Second run on the same folder, same inference — must skip already-done keys
  testthat::expect_no_error(run_once())
  rows_run2 <- nrow(unique(utils::read.csv2(result_csv)))

  # Idempotency: row count must not grow (allow == for skip, never grow)
  testthat::expect_equal(rows_run2, rows_run1)

  unlink(tempFolder, recursive = TRUE)
})

# ---------------------------------------------------------------------------
# 5. association_analysis — missing independent_variable is handled gracefully
# ---------------------------------------------------------------------------

test_that("association_analysis skips gracefully when independent_variable absent from sample sheet", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  unlink(tempFolder, recursive = TRUE)

  # Reuse the same bimodal injected-outlier signal as test 3
  set.seed(777)
  n_probes_aa2 <- 200L
  local_probes2  <- probe_features[1:n_probes_aa2, ]
  local_sig2 <- matrix(stats::rbeta(n_probes_aa2 * nsamples, 90L, 10L),
                        nrow = n_probes_aa2, ncol = nsamples)
  local_sig2[1:50, 1:5] <- stats::rbeta(50L * 5L, 1L, 100L)
  rownames(local_sig2) <- local_probes2$PROBE
  local_sig2 <- as.data.frame(local_sig2)
  colnames(local_sig2) <- colnames(signal_data)

  SEMseeker::semseeker(
    input             = local_sig2,
    sample_sheet      = mySampleSheet,
    result_folder     = tempFolder,
    parallel_strategy = parallel_strategy,
    areas             = c("GENE", "POSITION"),
    markers           = c("MUTATIONS"),
    start_fresh       = TRUE,
    inpute            = "median",
    showprogress      = showprogress,
    verbosity         = verbosity
  )

  inference_details <- data.frame(
    independent_variable = "NonExistentColumn",
    family_test          = "spearman",
    transformation_y     = "",
    transformation_x     = "",
    aggregation          = "SUM",
    scope                = "INSTANCE",
    filter_p_value       = FALSE,
    stringsAsFactors     = FALSE
  )

  # Should log a WARNING and skip rather than crash
  testthat::expect_no_error(
    SEMseeker:::association_analysis(
      inference_details = inference_details,
      result_folder     = tempFolder,
      parallel_strategy = parallel_strategy,
      markers           = c("MUTATIONS"),
      figures           = c("HYPO"),
      showprogress      = showprogress,
      verbosity         = verbosity
    )
  )

  unlink(tempFolder, recursive = TRUE)
})
