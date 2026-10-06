# 2026-06-09: contract-level E2E tests for the
# bedmethyl parser (modkit / nanopolish output). bedmethyl is the
# canonical long-read methylation-call file format; SEMseeker accepts
# it through `io_bedmethyl_read()` which converts a directory of modkit
# files into the wide SEMseeker coordinate-indexed data frame.
#
# Modkit bedmethyl schema (tab-separated, no header):
#   1  chrom
#   2  start_position
#   3  end_position
#   4  modified_base_code         (e.g. "m" for 5mC)
#   5  score
#   6  strand
#   7  start_code
#   8  end_code
#   9  color
#  10  N_valid_cov                (coverage at the position)
#  11  fraction_modified          (PERCENT, 0–100)
#  12+ additional modkit columns  (ignored by SEMseeker)
#
# Tests below generate synthetic modkit files in tempdir() and assert
# the contract that `io_bedmethyl_read()` must satisfy.

# ---- helper: write a synthetic modkit bedmethyl file ----------------------

.write_synthetic_bedmethyl <- function(path, rows) {
  # `rows` is a data.frame with columns
  # chrom, start, end, N_valid_cov, fraction_modified_percent.
  # We pad the remaining modkit columns with sensible placeholders.
  out <- data.frame(
    chrom              = rows$chrom,
    start              = rows$start,
    end                = rows$end,
    mod                = "m",
    score              = 0L,
    strand             = "+",
    start_code         = rows$start,
    end_code           = rows$end,
    color              = "0,0,0",
    N_valid_cov        = rows$N_valid_cov,
    fraction_modified  = rows$fraction_modified_percent,
    stringsAsFactors   = FALSE
  )
  utils::write.table(out, path, sep = "\t",
                     row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# ---- single-file read: schema + value-range contract --------------------

test_that("io_bedmethyl_read parses a single modkit file and returns SEMseeker shape", {
  skip_on_cran()

  td <- tempfile("ai109_bedm_single_")
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  rows <- data.frame(
    chrom                       = c("chr1", "chr1", "chr2"),
    start                       = c(100L, 200L, 300L),
    end                         = c(101L, 201L, 301L),
    N_valid_cov                 = c(10L, 8L, 12L),
    fraction_modified_percent   = c(75, 33, 100),
    stringsAsFactors            = FALSE
  )
  fpath <- file.path(td, "SAMPLE_A.bed")
  .write_synthetic_bedmethyl(fpath, rows)

  out <- SEMseeker:::io_bedmethyl_read(fpath)

  # Shape: CHR, START, END + 1 sample column
  expect_equal(colnames(out), c("CHR", "START", "END", "SAMPLE_A"))
  expect_equal(nrow(out), 3L)

  # Coordinates carried through verbatim
  expect_equal(out$CHR,   rows$chrom)
  expect_equal(out$START, rows$start)
  expect_equal(out$END,   rows$end)

  # Values normalised to [0, 1]
  expect_equal(out$SAMPLE_A, c(0.75, 0.33, 1.00))
  expect_true(all(out$SAMPLE_A >= 0 & out$SAMPLE_A <= 1))
})

# ---- coverage filter ----------------------------------------------------

test_that("io_bedmethyl_read drops positions below min_coverage", {
  skip_on_cran()

  td <- tempfile("ai109_bedm_cov_")
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  rows <- data.frame(
    chrom                       = c("chr1", "chr1", "chr1"),
    start                       = c(100L, 200L, 300L),
    end                         = c(101L, 201L, 301L),
    N_valid_cov                 = c(10L,  2L, 100L),    # middle row is below default 5
    fraction_modified_percent   = c(50,   90, 25),
    stringsAsFactors            = FALSE
  )
  fpath <- file.path(td, "S001.bed")
  .write_synthetic_bedmethyl(fpath, rows)

  # Default min_coverage = 5: drops the middle row
  out_default <- SEMseeker:::io_bedmethyl_read(fpath)
  expect_equal(nrow(out_default), 2L)
  expect_equal(out_default$START, c(100L, 300L))

  # min_coverage = 10: keeps only rows with coverage >= 10
  out_strict <- SEMseeker:::io_bedmethyl_read(fpath, min_coverage = 10L)
  expect_equal(nrow(out_strict), 2L)
  expect_equal(out_strict$START, c(100L, 300L))

  # min_coverage = 11: only the row with coverage 100 survives
  out_strict2 <- SEMseeker:::io_bedmethyl_read(fpath, min_coverage = 11L)
  expect_equal(nrow(out_strict2), 1L)
  expect_equal(out_strict2$START, 300L)
})

# ---- multi-sample outer join --------------------------------------------

test_that("io_bedmethyl_read outer-joins multiple files on (CHR, START, END)", {
  skip_on_cran()

  td <- tempfile("ai109_bedm_multi_")
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  rowsA <- data.frame(
    chrom = "chr1", start = c(100L, 200L), end = c(101L, 201L),
    N_valid_cov = 10L, fraction_modified_percent = c(50, 75),
    stringsAsFactors = FALSE
  )
  rowsB <- data.frame(
    chrom = "chr1", start = c(200L, 300L), end = c(201L, 301L),
    N_valid_cov = 10L, fraction_modified_percent = c(80, 25),
    stringsAsFactors = FALSE
  )
  fA <- file.path(td, "SAMP_A.bed")
  fB <- file.path(td, "SAMP_B.bed")
  .write_synthetic_bedmethyl(fA, rowsA)
  .write_synthetic_bedmethyl(fB, rowsB)

  out <- SEMseeker:::io_bedmethyl_read(c(fA, fB))

  # Outer-join: 3 unique positions (100 only in A, 200 in both, 300 only in B)
  expect_equal(colnames(out), c("CHR", "START", "END", "SAMP_A", "SAMP_B"))
  expect_equal(nrow(out), 3L)
  expect_equal(out$START, c(100L, 200L, 300L))

  # Position 100: only in SAMP_A → SAMP_B should be NA
  expect_equal(out[out$START == 100L, "SAMP_A"], 0.50)
  expect_true(is.na(out[out$START == 100L, "SAMP_B"]))

  # Position 200: in both
  expect_equal(out[out$START == 200L, "SAMP_A"], 0.75)
  expect_equal(out[out$START == 200L, "SAMP_B"], 0.80)

  # Position 300: only in SAMP_B → SAMP_A should be NA
  expect_true(is.na(out[out$START == 300L, "SAMP_A"]))
  expect_equal(out[out$START == 300L, "SAMP_B"], 0.25)
})

# ---- custom sample IDs override basename --------------------------------

test_that("io_bedmethyl_read uses caller-supplied sample_ids over filename", {
  skip_on_cran()

  td <- tempfile("ai109_bedm_ids_")
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  rows <- data.frame(
    chrom = "chr1", start = 100L, end = 101L,
    N_valid_cov = 10L, fraction_modified_percent = 60,
    stringsAsFactors = FALSE
  )
  fA <- file.path(td, "weird-name with spaces.bed")
  .write_synthetic_bedmethyl(fA, rows)

  out <- SEMseeker:::io_bedmethyl_read(fA, sample_ids = "TCGA_AB_001")
  expect_equal(colnames(out), c("CHR", "START", "END", "TCGA_AB_001"))
  expect_equal(out$TCGA_AB_001, 0.60)
})

# ---- error handling: missing file ---------------------------------------

test_that("io_bedmethyl_read errors clearly on missing file path", {
  expect_error(
    SEMseeker:::io_bedmethyl_read(c("/tmp/this/does/not/exist.bed")),
    "file\\(s\\) not found"
  )
})

test_that("io_bedmethyl_read errors on sample_ids length mismatch", {
  skip_on_cran()

  td <- tempfile("ai109_bedm_lenmis_")
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  rows <- data.frame(
    chrom = "chr1", start = 100L, end = 101L,
    N_valid_cov = 10L, fraction_modified_percent = 60,
    stringsAsFactors = FALSE
  )
  fA <- file.path(td, "A.bed")
  .write_synthetic_bedmethyl(fA, rows)

  expect_error(
    SEMseeker:::io_bedmethyl_read(fA, sample_ids = c("X", "Y")),
    "sample_ids length must match"
  )
})

# ---- stable ordering: sort by CHR then START ----------------------------

test_that("io_bedmethyl_read output is sorted by CHR then START", {
  skip_on_cran()

  td <- tempfile("ai109_bedm_sort_")
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  # Write rows OUT of order on purpose
  rows <- data.frame(
    chrom                       = c("chr2", "chr1", "chr2", "chr1"),
    start                       = c(100L, 500L, 200L, 100L),
    end                         = c(101L, 501L, 201L, 101L),
    N_valid_cov                 = c(10L, 10L, 10L, 10L),
    fraction_modified_percent   = c(50, 25, 75, 100),
    stringsAsFactors            = FALSE
  )
  fpath <- file.path(td, "S.bed")
  .write_synthetic_bedmethyl(fpath, rows)

  out <- SEMseeker:::io_bedmethyl_read(fpath)

  # Sorted: chr1:100, chr1:500, chr2:100, chr2:200
  expect_equal(out$CHR,   c("chr1", "chr1", "chr2", "chr2"))
  expect_equal(out$START, c(100L,   500L,   100L,   200L))
})

# ---- empty file_paths fails fast ----------------------------------------

test_that("io_bedmethyl_read errors when no file_paths are provided", {
  expect_error(
    SEMseeker:::io_bedmethyl_read(character(0)),
    "no file_paths provided"
  )
})

# ---- integration: bedmethyl → io_coord_probe_features round-trip -----------
#
# bedmethyl gives us CHR/START/END. Downstream SEMseeker code expects
# PROBE IDs in "{CHR}_{START}" form (see R/coord_input.R). Verify the
# bedmethyl output can be converted to that form and that
# io_coord_probe_features() round-trips coordinates back unchanged.

test_that("bedmethyl output is compatible with io_coord_probe_features", {
  skip_on_cran()

  td <- tempfile("ai109_bedm_pf_")
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  rows <- data.frame(
    chrom                       = c("1", "1", "X"),
    start                       = c(10000L, 20000L, 99999L),
    end                         = c(10001L, 20001L, 100000L),
    N_valid_cov                 = c(10L, 10L, 10L),
    fraction_modified_percent   = c(50, 75, 25),
    stringsAsFactors            = FALSE
  )
  fpath <- file.path(td, "SAMP.bed")
  .write_synthetic_bedmethyl(fpath, rows)

  out <- SEMseeker:::io_bedmethyl_read(fpath)

  # Build synthetic probe IDs in SEMseeker {CHR}_{START} format
  probe_ids <- paste(out$CHR, out$START, sep = "_")
  expect_equal(probe_ids, c("1_10000", "1_20000", "X_99999"))

  pf <- SEMseeker:::io_coord_probe_features(probe_ids)
  # Round-trip: io_coord_probe_features recovers the same CHR/START/END
  expect_equal(pf$CHR,   out$CHR)
  expect_equal(pf$START, out$START)
  expect_equal(pf$END,   out$START + 1L)
})
