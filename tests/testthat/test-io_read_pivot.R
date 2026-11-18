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

test_that("io_stream_merge_bed() builds a lazy frame from minimal bed fixtures", {

  skip_on_cran()

  tempFolder <- tempfile("read_pivot_test_")
  bed_dir <- file.path(tempFolder, "Data", "Case", "MUTATIONS_HYPER")
  dir.create(bed_dir, recursive = TRUE)

  # Two tiny bedgraph files, sample S001 and S002, overlapping coordinates
  bg1 <- file.path(bed_dir, "S001_MUTATIONS_HYPER.bedgraph")
  bg2 <- file.path(bed_dir, "S002_MUTATIONS_HYPER.bedgraph")
  writeLines(c("chr1\t100\t101\t0.5", "chr1\t200\t201\t0.7"), bg1)
  writeLines(c("chr1\t100\t101\t0.3", "chr1\t300\t301\t0.9"), bg2)

  merged <- SEMseeker:::io_stream_merge_bed(c(bg1, bg2), "MUTATIONS", "HYPER")
  testthat::expect_s3_class(merged, "polars_lazy_frame")

  out <- merged$collect()
  testthat::expect_true(all(c("CHR","START","END","S001","S002") %in% colnames(out)))
  testthat::expect_equal(nrow(out), 3L)   # 100, 200, 300 = 3 distinct rows
})
