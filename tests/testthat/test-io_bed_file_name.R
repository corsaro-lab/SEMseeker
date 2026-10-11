# test-io_bed_file_name.R
# Tests for BED / bedgraph file generation functions
#
# Covered:
#   - io_dump_sample_as_bed_file()  writes a tab-separated BED file;
#                                prepends "chr" if absent; sorts by CHR/START/END
#
# The companion chart test (test-plot_box_plot.R) already covers plot_box_plot() PNG
# generation. The claim that used to stand here, that plot_manhattan_plot_per_area()
# was tested end to end by test-6 and test-7, was false: nothing in the package
# called it, so no pipeline reached it, and it could not have drawn anything if one
# had. It is now sem_marker_value_per_area_plot(), with its own tests.
# The tests here focus on the low-level I/O helper used by sem_analyze_single_sample().

# ---------------------------------------------------------------------------
# Helper: build a minimal BED-like data.frame
# ---------------------------------------------------------------------------
.make_bed_df <- function(n = 5L, add_chr_prefix = TRUE) {
  chrs <- if (add_chr_prefix) paste0("chr", 1:n) else as.character(1:n)
  data.frame(
    CHR   = chrs,
    START = as.numeric(seq(1000L, by = 1000L, length.out = n)),
    END   = as.numeric(seq(1001L, by = 1000L, length.out = n)),
    VALUE = round(stats::runif(n), 4),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# 1. Basic file creation
# ---------------------------------------------------------------------------

test_that("io_dump_sample_as_bed_file: creates the output file", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  set.seed(1)
  df   <- .make_bed_df(5L)
  path <- file.path(tf, "test_output.bed")
  SEMseeker:::io_dump_sample_as_bed_file(data_to_dump = df, fileName = path)
  expect_true(file.exists(path))
})

test_that("io_dump_sample_as_bed_file: output has correct number of rows", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  set.seed(2)
  n    <- 8L
  df   <- .make_bed_df(n)
  path <- file.path(tf, "test_rows.bed")
  SEMseeker:::io_dump_sample_as_bed_file(data_to_dump = df, fileName = path)

  written <- read.table(path, sep = "\t", header = FALSE)
  expect_equal(nrow(written), n)
})

test_that("io_dump_sample_as_bed_file: prepends 'chr' when chromosome has no prefix", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df   <- .make_bed_df(3L, add_chr_prefix = FALSE)  # CHR = "1", "2", "3"
  path <- file.path(tf, "test_chr_prefix.bed")
  SEMseeker:::io_dump_sample_as_bed_file(data_to_dump = df, fileName = path)

  written <- read.table(path, sep = "\t", header = FALSE)
  expect_true(all(startsWith(as.character(written[[1]]), "chr")))
})

test_that("io_dump_sample_as_bed_file: does NOT double-prepend 'chr'", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df   <- .make_bed_df(3L, add_chr_prefix = TRUE)  # CHR = "chr1", "chr2", "chr3"
  path <- file.path(tf, "test_no_double_chr.bed")
  SEMseeker:::io_dump_sample_as_bed_file(data_to_dump = df, fileName = path)

  written <- read.table(path, sep = "\t", header = FALSE)
  expect_false(any(startsWith(as.character(written[[1]]), "chrchr")))
})

test_that("io_dump_sample_as_bed_file: output is sorted by CHR then START", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # Deliberately unsorted input
  df <- data.frame(
    CHR   = c("chr2", "chr1", "chr2", "chr1"),
    START = as.numeric(c(3000, 2000, 1000, 1000)),
    END   = as.numeric(c(3001, 2001, 1001, 1001)),
    VALUE = c(0.1, 0.2, 0.3, 0.4),
    stringsAsFactors = FALSE
  )
  path <- file.path(tf, "test_sorted.bed")
  SEMseeker:::io_dump_sample_as_bed_file(data_to_dump = df, fileName = path)

  written <- read.table(path, sep = "\t", header = FALSE,
                         colClasses = c("character", "numeric", "numeric", "numeric"))
  # After sorting: chr1/1000, chr1/2000, chr2/1000, chr2/3000
  expect_equal(as.character(written[[1]]), c("chr1", "chr1", "chr2", "chr2"))
  expect_equal(written[[2]], c(1000, 2000, 1000, 3000))
})

test_that("io_dump_sample_as_bed_file: empty data.frame does not create file", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df_empty <- data.frame(CHR = character(0), START = numeric(0),
                          END = numeric(0), VALUE = numeric(0))
  path <- file.path(tf, "test_empty.bed")
  SEMseeker:::io_dump_sample_as_bed_file(data_to_dump = df_empty, fileName = path)
  # Empty data.frame → plyr::empty returns TRUE → file should NOT be written
  expect_false(file.exists(path))
})

test_that("io_dump_sample_as_bed_file: rows with NA START are dropped", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df <- data.frame(
    CHR   = c("chr1", "chr1", "chr1"),
    START = c(1000, NA_real_, 3000),
    END   = c(1001, NA_real_, 3001),
    VALUE = c(0.5, 0.6, 0.7),
    stringsAsFactors = FALSE
  )
  path <- file.path(tf, "test_na_start.bed")
  SEMseeker:::io_dump_sample_as_bed_file(data_to_dump = df, fileName = path)
  written <- read.table(path, sep = "\t", header = FALSE)
  expect_equal(nrow(written), 2L)  # NA row dropped
})

test_that("io_dump_sample_as_bed_file: gz extension — readr writes without error", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df   <- .make_bed_df(4L)
  path <- file.path(tf, "test_output.bed")
  # Even without .gz extension, the function writes via readr::write_tsv which
  # auto-detects gz from the path if needed. Without .gz, it writes plain text.
  expect_no_error(
    SEMseeker:::io_dump_sample_as_bed_file(data_to_dump = df, fileName = path)
  )
  expect_true(file.exists(path))
})
