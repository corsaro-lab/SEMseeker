## sample identifier normalisation must be forced on BOTH sides
## (sample sheet + signal matrix) before any name-based column subsetting.
##
## Regression: SEMseeker cleaned `colnames(signal_data)` with
## core_name_cleaning() but compared/subset them with the RAW `Sample_ID`
## values. With identifiers such as "C3L-00001-06" (-> "C3L_00001_06") the
## intersection is empty: util_exploratory_analysis() wrote a "cleaned" parquet
## containing only the PROBE column, and the SEM path failed downstream with
## "undefined columns selected". Invisible on the bundled fixtures because
## GSM\d+ identifiers make core_name_cleaning() a no-op.
##

# ---------------------------------------------------------------------------
# core_normalize_sample_ids  (pure)
# ---------------------------------------------------------------------------

test_that("core_normalize_sample_ids normalises sheet and signal columns coherently", {
  ss  <- data.frame(Sample_ID = c("C3L-00001-06", "C3L-00002-06"),
                    Sample_Group = c("Case", "Reference"),
                    stringsAsFactors = FALSE)
  sig <- data.frame(`C3L-00001-06` = c(0.1, 0.2), `C3L-00002-06` = c(0.3, 0.4),
                    check.names = FALSE)

  out <- SEMseeker:::core_normalize_sample_ids(ss, sig)

  expect_equal(out$sample_sheet$Sample_ID, c("C3L_00001_06", "C3L_00002_06"))
  expect_equal(colnames(out$signal_data), c("C3L_00001_06", "C3L_00002_06"))
  # the whole point: the two sides match after normalisation
  expect_true(all(out$sample_sheet$Sample_ID %in% colnames(out$signal_data)))
  expect_length(out$unmatched_ids, 0)
  expect_length(out$unmatched_columns, 0)
})

test_that("core_normalize_sample_ids is idempotent", {
  ss  <- data.frame(Sample_ID = c("C3L-00001-06", "GSM123456"),
                    stringsAsFactors = FALSE)
  sig <- data.frame(`C3L-00001-06` = 1, `GSM123456` = 2, check.names = FALSE)

  once  <- SEMseeker:::core_normalize_sample_ids(ss, sig)
  twice <- SEMseeker:::core_normalize_sample_ids(once$sample_sheet, once$signal_data)

  expect_equal(twice$sample_sheet$Sample_ID, once$sample_sheet$Sample_ID)
  expect_equal(colnames(twice$signal_data), colnames(once$signal_data))
  # GSM identifiers are untouched — this is why the bug stayed hidden
  expect_true("GSM123456" %in% once$sample_sheet$Sample_ID)
})

test_that("core_normalize_sample_ids leaves structural columns untouched", {
  sig <- data.frame(PROBE = "cg0001", CHR = "chr1", START = 1, END = 2,
                    `C3L-00001-06` = 0.5, check.names = FALSE)

  out <- SEMseeker:::core_normalize_sample_ids(NULL, sig)

  expect_equal(colnames(out$signal_data),
               c("PROBE", "CHR", "START", "END", "C3L_00001_06"))
})

test_that("core_normalize_sample_ids fails when normalisation collapses distinct ids", {
  ss <- data.frame(Sample_ID = c("S-1", "S.1"), stringsAsFactors = FALSE)

  expect_error(SEMseeker:::core_normalize_sample_ids(ss),
               "collapses distinct")
})

test_that("core_normalize_sample_ids reports orphan identifiers instead of 'undefined columns selected'", {
  ss  <- data.frame(Sample_ID = c("C3L-00001-06", "C3L-00099-06"),
                    stringsAsFactors = FALSE)
  sig <- data.frame(`C3L-00001-06` = c(0.1, 0.2), check.names = FALSE)

  expect_error(
    SEMseeker:::core_normalize_sample_ids(ss, sig, require_all_ids = TRUE),
    "C3L_00099_06"
  )
  # non-strict mode reports them without failing
  quiet <- SEMseeker:::core_normalize_sample_ids(ss, sig)
  expect_equal(quiet$unmatched_ids, "C3L_00099_06")
})

test_that("core_normalize_sample_ids fails loudly when the id column is missing", {
  expect_error(
    SEMseeker:::core_normalize_sample_ids(data.frame(ID = "a", stringsAsFactors = FALSE)),
    "no 'Sample_ID' column"
  )
})

# ---------------------------------------------------------------------------
# util_exploratory_analysis  (end-to-end regression on dashed identifiers)
# ---------------------------------------------------------------------------

test_that("util_exploratory_analysis keeps sample columns with non-alphanumeric Sample_IDs", {
  tf <- sem_test_folder()
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE); unlink(tf, recursive = TRUE) },
          add = TRUE)

  # CPTAC-style identifiers: core_name_cleaning() turns "-" into "_"
  raw_ids <- paste0("C3L-", sprintf("%05d", seq_len(ncol(signal_data))), "-06")
  dashed_signal <- signal_data
  colnames(dashed_signal) <- raw_ids

  dashed_sheet <- mySampleSheet
  dashed_sheet$Sample_ID <- raw_ids[match(dashed_sheet$Sample_ID, colnames(signal_data))]
  dashed_sheet <- dashed_sheet[!is.na(dashed_sheet$Sample_ID), ]

  expected_ids <- unique(SEMseeker:::core_name_cleaning(dashed_sheet$Sample_ID))

  expect_no_error(
    SEMseeker:::util_exploratory_analysis(
      categorical_variables = c("Sample_Group"),
      numerical_variables   = c("Phenotest", "Covariates1"),
      sample_sheet          = dashed_sheet,
      signal_data           = dashed_signal,
      result_folder         = tf,
      parallel_strategy     = "sequential",
      start_fresh           = TRUE,
      showprogress          = FALSE,
      verbosity             = 1
    )
  )

  explo_dir <- file.path(tf, "Data", "Exploratory_0")
  expect_true(dir.exists(explo_dir))

  parquet_path <- list.files(explo_dir, pattern = "CLEANED_SIGNAL_DATA.*\\.parquet$",
                             full.names = TRUE, ignore.case = TRUE)
  expect_length(parquet_path, 1)

  cleaned <- as.data.frame(polars::pl$read_parquet(parquet_path[1]))
  sample_cols <- setdiff(colnames(cleaned), "PROBE")

  # THE regression: before the fix this parquet held the PROBE column only
  expect_gt(length(sample_cols), 0)
  expect_true(all(expected_ids %in% sample_cols))

  # ... and the cleaned sample sheet must address exactly those columns
  sheet_path <- list.files(explo_dir, pattern = "CLEANED_SAMPLE_SHEET.*\\.csv$",
                           full.names = TRUE, ignore.case = TRUE)
  expect_length(sheet_path, 1)
  cleaned_sheet <- utils::read.csv2(sheet_path[1], stringsAsFactors = FALSE)
  expect_setequal(unique(as.character(cleaned_sheet$Sample_ID)), sample_cols)
})
