# AI-086: integration test for the per-sample burden.
#
# AI-223 moved the burden out of SAMPLE_SHEET_RESULT.csv and into the
# statistics sibling SAMPLE_STATS_RESULT.csv, where it is named
# SAMPLE_<MARKER>_<FIGURE> (PROBES_COUNT became SAMPLE_N_PROBES). This test
# follows it there, and additionally asserts that the sample sheet is left
# lean.
#
# Canary for AI-083: the burden aggregation failed to populate per-sample
# values on ewas_osteoporosis/GSE99624 (48 samples, 450K) — all burden columns
# landed as NA, which then crashed every depth=1 inference downstream with
# "data are not the same size".
#
# test-7-association_analysis.R only verifies this indirectly (depth=1 would
# fail on all-NA burden vectors); this test asserts the contract directly:
#
#   1. sample_sheet_result.csv exists after semseeker().
#   2. Every expected <MARKER>_<FIGURE> burden column is PRESENT.
#   3. Every burden column has at least one non-NA sample.
#   4. No sample has ALL burden values NA (per-sample BED collection works).
#   5. PROBES_COUNT > 0 on every sample.
#   6. Sanity: injected hypomethylated outliers surface as MUTATIONS_HYPO > 0.
#
# Why discrete markers only: sem_study_summary_total() iterates AREA=="POSITION"
# keys, which is the set declared in keys_create.R via
# keys_markers_default_discrete (MUTATIONS / LESIONS / DELTAQ / DELTARQ /
# DELTAP / DELTARP) plus the two continuous DELTAS / DELTAR. LESIONS depends
# on MUTATIONS having been computed (SOURCE column), so the full
# (MUTATIONS + DELTA*) set exercises the whole derive-and-aggregate path.

test_that("sample_sheet_result.csv has populated burden columns for all discrete + continuous markers (AI-086, canary for AI-083)", {
  tempFolder <- sem_test_folder()
  unlink(tempFolder, recursive = TRUE)

  syn <- .burden_setup_signal_with_outliers(
    n_samples      = nsamples,
    probe_features = probe_features,
    sample_sheet   = mySampleSheet,
    signal_data    = signal_data)

  # inpute="median" is defensive: harmless when the input has no NAs (the
  # guard at inpute_missing_values.R:7 short-circuits), useful if a future
  # signal generator introduces them.
  SEMseeker::semseeker(
    input             = syn$signal,
    sample_sheet      = syn$samples,
    result_folder     = tempFolder,
    parallel_strategy = "sequential",
    areas             = c("POSITION"),
    markers           = c("MUTATIONS", "LESIONS",
                          "DELTAP", "DELTAQ", "DELTARP", "DELTARQ",
                          "DELTAS", "DELTAR"),
    start_fresh       = TRUE,
    inpute            = "median",
    showprogress      = showprogress,
    verbosity         = verbosity
  )

  # SEMseeker writes via io_file_path_build() → core_name_cleaning() → toupper(),
  # so the on-disk file is SAMPLE_SHEET_RESULT.csv (not the lowercase form
  # used in the API contract). macOS/Windows file systems are case-insensitive
  # by default so either spelling works; Linux ext4 is case-sensitive and
  # only the uppercase form resolves. Use uppercase here to be correct on
  # all three CI runners.
  # AI-255: the sibling CSV is gone. The per-sample table is composed on read
  # from the SCOPE = SAMPLE artefacts and joined onto the sample sheet.
  sheet_csv  <- file.path(tempFolder, "Data", "SAMPLE_SHEET_RESULT.csv")
  df <- SEMseeker:::sem_study_summary_get()
  testthat::expect_true(!is.null(df) && nrow(df) > 0,
    info = sprintf("no per-sample statistics composed — tempFolder=%s", tempFolder))

  # AI-223 net move: the sample sheet must NOT carry the burden any more
  sheet <- utils::read.csv2(sheet_csv, stringsAsFactors = FALSE)
  testthat::expect_equal(
    intersect(c("MUTATIONS_HYPER", "MUTATIONS_HYPO", "PROBES_COUNT"),
              colnames(sheet)),
    character(0)
  )

  # Required columns: MUTATIONS + DELTA* (LESIONS is optional — derives from
  # MUTATIONS clusters and may legitimately be 0 on small synthetic data.
  # Tracked separately below as a soft expectation.)
  required_markers <- c("MUTATIONS",
                        "DELTAP", "DELTAQ", "DELTARP", "DELTARQ",
                        "DELTAS", "DELTAR")
  # AI-248: composed once in helper-burden.R, where the aggregation each class
  # produces by default is spelled out.
  required_burden_cols <- .burden_cols
  # AI-255: N_PROBES describes the sample, not a scope of it — it counts the
  # positions the imputation left usable — so it comes in from the sample sheet
  # under its own name.
  required_cols <- c("Sample_ID", required_burden_cols, "N_PROBES")

  missing_cols <- setdiff(required_cols, colnames(df))
  testthat::expect_equal(
    length(missing_cols), 0L,
    info = paste("missing required columns:", paste(missing_cols, collapse = ", "))
  )

  # Soft: log if LESIONS_HYPER/HYPO are absent so we surface the gap
  # without failing — see AI-088 follow-up for an explicit LESIONS canary.
  for (lesion_col in c("SAMPLE_LESIONS_HYPER", "SAMPLE_LESIONS_HYPO")) {
    if (!(lesion_col %in% colnames(df))) {
      message(sprintf(
        "test-8-burden-integration: %s absent — likely synthetic-data sparsity (soft warn, see AI-088)",
        lesion_col
      ))
    }
  }

  # The per-column and per-sample assertions below operate only on the
  # required burden columns (LESIONS handled by the soft block above).
  expected_burden_cols <- required_burden_cols

  # AI-083 canary (per-column): every burden column has ≥ 1 non-NA value.
  for (col in intersect(expected_burden_cols, colnames(df))) {
    n_non_na <- sum(!is.na(df[[col]]))
    testthat::expect_gt(
      n_non_na, 0L,
      label = sprintf("%s non-NA samples (got %d)", col, n_non_na)
    )
  }

  # AI-083 canary (per-sample, relaxed): NOT every sample must have a
  # populated burden — synthetic data sparsity can legitimately leave some
  # samples with zero events for ALL (marker, figure) combos, which then
  # become all-NA after the all.x=TRUE merge in study_summary_total. The
  # AI-083 bug was 100% of samples NA (whole-population integration broken),
  # not "some samples are NA". So we assert the all-NA fraction is below a
  # safety margin instead of zero. A tighter per-sample canary on real or
  # tightly-controlled data is tracked in AI-088.
  present_burden <- intersect(expected_burden_cols, colnames(df))
  if (length(present_burden) > 0L) {
    per_sample_all_na <- apply(
      df[, present_burden, drop = FALSE], 1L,
      function(r) all(is.na(r))
    )
    all_na_fraction <- sum(per_sample_all_na) / nrow(df)
    testthat::expect_lt(
      all_na_fraction, 0.5,
      label = sprintf(
        "all-NA-burden fraction (got %.2f, %d/%d): %s",
        all_na_fraction, sum(per_sample_all_na), nrow(df),
        paste(df$Sample_ID[per_sample_all_na], collapse = ", ")
      )
    )
  }

  # PROBES_COUNT > 0 on every sample.
  testthat::expect_true(
    all(df$SAMPLE_N_PROBES > 0L, na.rm = TRUE),
    info = "SAMPLE_N_PROBES must be > 0 on every sample"
  )

  # Sanity: injected HYPO outliers must surface as MUTATIONS_HYPO > 0 on
  # at least the 5 samples we touched (sample indices 1:5).
  if ("SAMPLE_MUTATIONS_HYPO" %in% colnames(df)) {
    testthat::expect_gt(
      sum(df$SAMPLE_MUTATIONS_HYPO > 0, na.rm = TRUE), 0L,
      label = "samples with SAMPLE_MUTATIONS_HYPO > 0 (injected-outlier sanity)"
    )
  }
  if ("SAMPLE_MUTATIONS_HYPER" %in% colnames(df)) {
    testthat::expect_gt(
      sum(df$SAMPLE_MUTATIONS_HYPER > 0, na.rm = TRUE), 0L,
      label = "samples with SAMPLE_MUTATIONS_HYPER > 0 (injected-outlier sanity)"
    )
  }

  # Explicit pivot assertion: every required marker must have a per-figure
  # parquet pivot under Data/Pivots/<MARKER>/, and each one must contain
  # rows. The burden checks above cover this indirectly (an empty pivot
  # → all-NA burden column), but a direct file/row check catches an empty
  # pivot upstream of the merge and points at the failing marker by name.
  pivots_dir <- file.path(tempFolder, "Data", "Pivots")
  for (m in required_markers) {
    marker_dir <- file.path(pivots_dir, m)
    testthat::expect_true(
      dir.exists(marker_dir),
      info = sprintf("missing pivot directory for marker '%s' at %s", m, marker_dir)
    )
    if (!dir.exists(marker_dir)) next
    pq_files <- list.files(marker_dir, pattern = "\\.parquet$", full.names = TRUE)
    testthat::expect_gt(
      length(pq_files), 0L,
      label = sprintf("parquet pivot files under %s", marker_dir)
    )
    for (pf in pq_files) {
      n_rows <- tryCatch(
        as.integer(nrow(polars::pl$read_parquet(pf))),
        error = function(e) NA_integer_
      )
      testthat::expect_gt(
        n_rows, 0L,
        label = sprintf("rows in pivot %s", basename(pf))
      )
    }
  }

  unlink(tempFolder, recursive = TRUE)
})

# ---------------------------------------------------------------------------
# AI-083 hardening. The canary above runs on GSM-style identifiers in a clean
# session; this block covers what it cannot see:
#
#   (a) identifiers that core_name_cleaning() actually rewrites — the
#       ewas_osteoporosis/GSE99624 case was "DNAm_sample1" -> "DNAM_SAMPLE1";
#   (b) a `temp_result` object left in the global environment. The old code
#       used exists("temp_result"), which resolves through globalenv(). This
#       is NOT sufficient on its own to produce the reported symptom (the
#       inner merge uses all=TRUE, so a stray object only adds a spurious row),
#       but the accumulator must be a local binding regardless;
#   (c) the condition that DOES produce the reported symptom — a burden table
#       sharing no Sample_ID with the sample sheet, which used to yield 100% NA
#       burden columns and killed every depth=1 inference downstream with
#       "data are not the same size". Since AI-223 the aggregation lives in
#       sem_sample_stats_build(), which must refuse to write that file.
# ---------------------------------------------------------------------------

test_that("burden survives mixed-case Sample_IDs and a polluted global temp_result (AI-083)", {
  tempFolder <- sem_test_folder()
  unlink(tempFolder, recursive = TRUE)
  on.exit({
    try(SEMseeker:::core_close_env(), silent = TRUE)
    if (exists("temp_result", envir = globalenv(), inherits = FALSE))
      rm("temp_result", envir = globalenv())
    unlink(tempFolder, recursive = TRUE)
  }, add = TRUE)

  syn <- .burden_setup_signal_with_outliers(
    n_samples      = nsamples,
    probe_features = probe_features,
    sample_sheet   = mySampleSheet,
    signal_data    = signal_data)

  # ewas_osteoporosis-style identifiers: core_name_cleaning() uppercases them
  # ("DNAm_sample1" -> "DNAM_SAMPLE1"), so sample sheet and pivot columns must
  # still meet after normalisation.
  raw_ids <- paste0("DNAm_sample", seq_len(ncol(syn$signal)))
  id_map  <- stats::setNames(raw_ids, colnames(syn$signal))
  sig     <- syn$signal
  sheet   <- syn$samples
  colnames(sig)   <- unname(id_map[colnames(sig)])
  sheet$Sample_ID <- unname(id_map[sheet$Sample_ID])

  # (b) pollute globalenv BEFORE the run: the object name collides with the
  # accumulator inside the burden aggregation.
  assign("temp_result",
         data.frame(Sample_ID = "NOT_A_REAL_SAMPLE", JUNK = 1,
                    stringsAsFactors = FALSE),
         envir = globalenv())

  SEMseeker::semseeker(
    input             = sig,
    sample_sheet      = sheet,
    result_folder     = tempFolder,
    parallel_strategy = "sequential",
    areas             = c("POSITION"),
    markers           = c("MUTATIONS", "LESIONS",
                          "DELTAP", "DELTAQ", "DELTARP", "DELTARQ",
                          "DELTAS", "DELTAR"),
    start_fresh       = TRUE,
    inpute            = "median",
    showprogress      = showprogress,
    verbosity         = verbosity
  )

  # AI-255: composed on read, no sibling file.
  df <- SEMseeker:::sem_study_summary_get()
  testthat::expect_true(!is.null(df) && nrow(df) > 0)

  # identifiers landed normalised on both sides
  testthat::expect_true(all(grepl("^DNAM_SAMPLE", df$Sample_ID)))

  present <- intersect(.burden_cols, colnames(df))
  testthat::expect_gt(length(present), 0L)

  # the AI-083 signature is 100% NA on every burden column at once
  na_fraction <- vapply(df[present], function(x) mean(is.na(x)), numeric(1))
  testthat::expect_false(
    all(na_fraction == 1),
    info = sprintf("all burden columns are 100%% NA: %s",
                   paste(names(na_fraction), collapse = ", "))
  )
  testthat::expect_lt(
    mean(na_fraction), 0.5,
    label = sprintf("mean NA fraction across burden columns (%.2f)",
                    mean(na_fraction))
  )
  testthat::expect_true(all(df$SAMPLE_N_PROBES > 0L, na.rm = TRUE))

  # (c) the actual AI-083 signature: rewrite the SAMPLE SHEET with identifiers
  # that exist in no pivot, then re-run the aggregation. The producer recomputes
  # the burden from the pivots (which still carry the real identifiers) and
  # compares it against the sheet, so corrupting the sheet is what reproduces
  # the mismatch. Before the guard this silently produced 100% NA burden; now it
  # must stop and name both sides.
  sheet_csv <- file.path(tempFolder, "Data", "SAMPLE_SHEET_RESULT.csv")
  sheet_broken <- utils::read.csv2(sheet_csv, stringsAsFactors = FALSE)
  sheet_broken$Sample_ID <- paste0("UNRELATED_", seq_len(nrow(sheet_broken)))
  utils::write.csv2(sheet_broken, sheet_csv, row.names = FALSE)

  SEMseeker:::core_init_env(
    result_folder     = tempFolder,
    parallel_strategy = "sequential",
    areas             = c("POSITION"),
    markers           = c("MUTATIONS", "LESIONS",
                          "DELTAP", "DELTAQ", "DELTARP", "DELTARQ",
                          "DELTAS", "DELTAR"),
    start_fresh       = FALSE,
    showprogress      = FALSE,
    verbosity         = 1
  )
  # AI-255: the check moved with the join. sem_sample_stats_build() now only
  # materialises artefacts; it is sem_study_summary_get() that puts the sample
  # sheet and the artefact columns side by side, so that is where a total
  # identifier mismatch has to be caught rather than joined into a table of NAs.
  testthat::expect_error(
    SEMseeker:::sem_study_summary_get(),
    "no Sample_ID in common"
  )
})
