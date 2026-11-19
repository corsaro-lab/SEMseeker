# io_read_pivot() dispatch on tiny fixtures: io_pivot_ensure() makes the file
# available, io_read_pivot() scans it with the data engine, and the polars
# bridge reads the same file for the callers not yet moved.

test_that("io_read_pivot() returns NULL when no storage exists", {

  skip_on_cran()

  tempFolder <- tempfile("read_pivot_test_")
  dir.create(file.path(tempFolder, "Data"), recursive = TRUE)
  ssEnv <- SEMseeker:::core_init_env(tempFolder,
                                parallel_strategy = "sequential",
                                iqrTimes = 3, verbosity = 1)
  on.exit({ SEMseeker:::core_close_env(); unlink(tempFolder, recursive = TRUE) },
          add = TRUE)

  testthat::expect_null(
    SEMseeker:::io_read_pivot("MUTATIONS", "HYPER")
  )
})

test_that("io_read_pivot() prefers cached parquet over bed files", {

  skip_on_cran()

  tempFolder <- tempfile("read_pivot_test_")
  dir.create(file.path(tempFolder, "Data"), recursive = TRUE)
  ssEnv <- SEMseeker:::core_init_env(tempFolder,
                                parallel_strategy = "sequential",
                                iqrTimes = 3, verbosity = 1)
  on.exit({ SEMseeker:::core_close_env(); unlink(tempFolder, recursive = TRUE) },
          add = TRUE)

  # Write a 1-row stub parquet at the expected pivot location
  pivot_path <- SEMseeker:::io_pivot_file_name_parquet("MUTATIONS", "HYPER",
                                                   "POSITION", "WHOLE")
  dir.create(dirname(pivot_path), recursive = TRUE, showWarnings = FALSE)
  polars::as_polars_df(data.frame(CHR = "1", START = 1L, END = 2L,
                                  SAMPLE_X = 0.0))$write_parquet(pivot_path)

  testthat::expect_identical(SEMseeker:::io_pivot_ensure("MUTATIONS", "HYPER"), pivot_path)

  result <- SEMseeker:::io_read_pivot("MUTATIONS", "HYPER")
  testthat::expect_s3_class(result, "tbl_lazy")
  testthat::expect_true("SAMPLE_X" %in% colnames(dplyr::collect(result)))

  bridged <- SEMseeker:::io_read_pivot_polars("MUTATIONS", "HYPER")
  testthat::expect_s3_class(bridged, "polars_lazy_frame")
  testthat::expect_identical(as.data.frame(bridged$collect()), as.data.frame(dplyr::collect(result)))
})

test_that("io_pivot_ensure() merges the bed files once and the merged file is the cache", {

  skip_on_cran()

  tempFolder <- tempfile("read_pivot_test_")
  ssEnv <- SEMseeker:::core_init_env(tempFolder,
                                parallel_strategy = "sequential",
                                iqrTimes = 3, verbosity = 1)
  on.exit({ SEMseeker:::core_close_env(); unlink(tempFolder, recursive = TRUE) },
          add = TRUE)
  bed_dir <- file.path(ssEnv$result_folderData, "Case", "MUTATIONS_HYPER")
  dir.create(bed_dir, recursive = TRUE)
  writeLines(c("chr1\t100\t101\t0.5", "chr1\t200\t201\t0.7"), file.path(bed_dir, "S001_MUTATIONS_HYPER.bedgraph"))
  writeLines(c("chr1\t100\t101\t0.3", "chr1\t300\t301\t0.9"), file.path(bed_dir, "S002_MUTATIONS_HYPER.bedgraph"))

  pivot_path <- SEMseeker:::io_pivot_ensure("MUTATIONS", "HYPER")
  testthat::expect_identical(pivot_path, SEMseeker:::io_pivot_file_name_parquet("MUTATIONS", "HYPER", "POSITION", "WHOLE"))
  testthat::expect_true(file.exists(pivot_path))

  merged <- as.data.frame(dplyr::collect(SEMseeker:::io_read_pivot("MUTATIONS", "HYPER")))
  testthat::expect_equal(nrow(merged), 3L)
  testthat::expect_true(all(c("CHR", "START", "END", "S001", "S002") %in% colnames(merged)))

  # The bed files are no longer needed: a second call reads the cached pivot.
  unlink(bed_dir, recursive = TRUE)
  testthat::expect_identical(SEMseeker:::io_pivot_ensure("MUTATIONS", "HYPER"), pivot_path)
})

test_that("io_stream_merge_bed() writes the wide pivot, samples in the order of the files", {

  skip_on_cran()
  skip_if_not_installed("nanoparquet")

  tempFolder <- tempfile("read_pivot_test_")
  bed_dir <- file.path(tempFolder, "Data", "Case", "MUTATIONS_HYPER")
  dir.create(bed_dir, recursive = TRUE)

  # S002 is given first: its column must come first, not in alphabetical order.
  bg2 <- file.path(bed_dir, "S002_MUTATIONS_HYPER.bedgraph")
  bg1 <- file.path(bed_dir, "S001_MUTATIONS_HYPER.bedgraph")
  writeLines(c("chr1\t100\t101\t0.3", "chr1\t300\t301\t0.9"), bg2)
  writeLines(c("chr1\t100\t101\t0.5", "chr2\t200\t201\t0.7"), bg1)
  pivot_file <- file.path(tempFolder, "merged.parquet")

  SEMseeker:::io_stream_merge_bed(c(bg2, bg1), "MUTATIONS", "HYPER", pivot_file)
  out <- as.data.frame(nanoparquet::read_parquet(pivot_file))

  testthat::expect_identical(colnames(out), c("CHR", "START", "END", "S002", "S001"))
  # One row per position, sorted on CHR then START, "chr" removed.
  testthat::expect_identical(out$CHR, c("1", "1", "2"))
  testthat::expect_equal(out$START, c(100, 300, 200))
  testthat::expect_equal(out$S002, c(0.3, 0.9, NA))
  testthat::expect_equal(out$S001, c(0.5, NA, 0.7))
})

test_that("io_stream_merge_bed() refuses a position repeated within one sample", {

  skip_on_cran()

  tempFolder <- tempfile("read_pivot_test_")
  bed_dir <- file.path(tempFolder, "Data", "Case", "MUTATIONS_HYPER")
  dir.create(bed_dir, recursive = TRUE)
  bg1 <- file.path(bed_dir, "S001_MUTATIONS_HYPER.bedgraph")
  # The same position twice, once with and once without the "chr" prefix.
  writeLines(c("chr1\t100\t101\t0.5", "1\t100\t101\t0.6"), bg1)

  testthat::expect_error(
    SEMseeker:::io_stream_merge_bed(bg1, "MUTATIONS", "HYPER", file.path(tempFolder, "merged.parquet")),
    "repeated within one sample")
})
