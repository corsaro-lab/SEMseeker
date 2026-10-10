# =============================================================================
# test-anno_normalize_chr.R
#
# Tests for SEMseeker:::anno_normalize_chr() — the single entry point for chromosome name
# normalization between internal (bare: "1", "X") and output (UCSC: "chr1",
# "chrX") conventions.
# =============================================================================

test_that("internal direction strips chr prefix", {
  expect_equal(SEMseeker:::anno_normalize_chr(c("chr1", "chr22", "chrX"), "internal"),
               c("1", "22", "X"))
})

test_that("internal direction is case-insensitive", {
  expect_equal(SEMseeker:::anno_normalize_chr(c("Chr1", "CHR2", "CHRX"), "internal"),
               c("1", "2", "X"))
})

test_that("internal direction is idempotent on bare names", {
  expect_equal(SEMseeker:::anno_normalize_chr(c("1", "X", "22"), "internal"),
               c("1", "X", "22"))
})

test_that("output direction adds chr prefix", {
  expect_equal(SEMseeker:::anno_normalize_chr(c("1", "22", "X"), "output"),
               c("chr1", "chr22", "chrX"))
})

test_that("output direction is idempotent on prefixed names", {
  expect_equal(SEMseeker:::anno_normalize_chr(c("chr1", "chrX"), "output"),
               c("chr1", "chrX"))
})

test_that("output direction preserves existing case of chr prefix", {
  expect_equal(SEMseeker:::anno_normalize_chr(c("Chr1", "CHR2"), "output"),
               c("Chr1", "CHR2"))
})

test_that("handles factor input", {
  f <- factor(c("chr1", "chr2", "chrX"))
  expect_equal(SEMseeker:::anno_normalize_chr(f, "internal"), c("1", "2", "X"))
})

test_that("handles integer input", {
  expect_equal(SEMseeker:::anno_normalize_chr(c(1L, 22L), "internal"), c("1", "22"))
  expect_equal(SEMseeker:::anno_normalize_chr(c(1L, 22L), "output"), c("chr1", "chr22"))
})

test_that("handles empty vector", {
  expect_equal(SEMseeker:::anno_normalize_chr(character(0), "internal"), character(0))
  expect_equal(SEMseeker:::anno_normalize_chr(character(0), "output"), character(0))
})

test_that("invalid direction raises error", {
  expect_error(SEMseeker:::anno_normalize_chr("chr1", "invalid"))
})

test_that("mixed input: some with chr, some without", {
  expect_equal(SEMseeker:::anno_normalize_chr(c("chr1", "2", "chrX"), "internal"),
               c("1", "2", "X"))
  expect_equal(SEMseeker:::anno_normalize_chr(c("chr1", "2", "chrX"), "output"),
               c("chr1", "chr2", "chrX"))
})

test_that("always returns character, never numeric", {
  expect_type(SEMseeker:::anno_normalize_chr(c(1, 22), "internal"), "character")
  expect_type(SEMseeker:::anno_normalize_chr(c(1, 22), "output"), "character")
  expect_type(SEMseeker:::anno_normalize_chr("chr1", "internal"), "character")
})

test_that("round-trip: internal → output → internal is identity", {
  bare <- c("1", "22", "X", "Y", "MT")
  expect_equal(SEMseeker:::anno_normalize_chr(SEMseeker:::anno_normalize_chr(bare, "output"), "internal"), bare)
})

test_that("round-trip: output → internal → output is identity", {
  prefixed <- c("chr1", "chr22", "chrX", "chrY", "chrMT")
  expect_equal(SEMseeker:::anno_normalize_chr(SEMseeker:::anno_normalize_chr(prefixed, "internal"), "output"),
               prefixed)
})

# ---------------------------------------------------------------------------
# Regression test: BED write → read round-trip must preserve joinability
#
# This is the exact scenario that caused covered_by_inner_join = 0:
# signal_thresholds has CHR = "1" (bare, from probe_features join),
# bedgraph files written by io_dump_sample_as_bed_file have CHR = "chr1".
# Reading back without SEMseeker:::anno_normalize_chr("internal") breaks the join.
# ---------------------------------------------------------------------------
test_that("BED round-trip: bare CHR survives write+read cycle", {
  skip_if_not_installed("readr")

  # Simulate internal data with bare CHR (as probe_features produces)
  internal_df <- data.frame(
    CHR   = c("1", "1", "X"),
    START = c(15865, 18827, 100000),
    END   = c(15865, 18827, 100000),
    VALUE = c(0.88, 0.56, 0.42),
    stringsAsFactors = FALSE
  )

  # Write as BED (adds chr prefix, as io_dump_sample_as_bed_file does)
  tmp <- tempfile(fileext = ".bed")
  on.exit(unlink(tmp), add = TRUE)
  bed_out <- internal_df
  bed_out$CHR <- SEMseeker:::anno_normalize_chr(bed_out$CHR, "output")
  readr::write_tsv(bed_out, tmp, col_names = FALSE)

  # Read back and normalize (as fixed sem_analyze_population does)
  bed_in <- utils::read.delim(tmp, header = FALSE, sep = "\t",
                               stringsAsFactors = FALSE)
  colnames(bed_in) <- c("CHR", "START", "END", "VALUE")
  bed_in$CHR <- SEMseeker:::anno_normalize_chr(bed_in$CHR, "internal")

  # The CHR values must match the original internal format
  expect_equal(bed_in$CHR, internal_df$CHR)

  # Inner join on CHR+START must find all rows (the original bug: 0 matches)
  joined <- merge(bed_in[, c("CHR", "START", "END")],
                  internal_df[, c("CHR", "START", "END")],
                  by = c("CHR", "START", "END"))
  expect_equal(nrow(joined), nrow(internal_df))
})
