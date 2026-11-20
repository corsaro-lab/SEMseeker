# The position pivot grows as new samples arrive: a second call merges the new
# samples' bed files into the pivot the first call wrote. The values below are
# written by hand, so the expected table is known without running the code
# under test.

#' Write one per-sample bed file where the pipeline would have written it.
#'
#' @param sample_id the sample, used by io_bed_file_name() for the path
#' @param sample_group its group, which is part of the path
#' @param positions data.frame with CHR, START, END and value
.write_sample_bed <- function(sample_id, sample_group, positions) {
  bed_path <- SEMseeker:::io_bed_file_name(sample_id, sample_group, "MUTATIONS", "HYPER")
  bed_connection <- if (grepl("\\.gz$", bed_path)) gzfile(bed_path, "w") else file(bed_path, "w")
  utils::write.table(positions, bed_connection, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  close(bed_connection)
  bed_path
}

test_that("a second call merges new samples into the existing position pivot", {
  skip_if_not_installed("nanoparquet")
  temp_folder <- sem_test_folder()
  unlink(temp_folder, recursive = TRUE, force = TRUE)
  SEMseeker:::core_init_env(temp_folder, parallel_strategy = "sequential", start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(temp_folder, recursive = TRUE) }, add = TRUE)
  keys <- data.frame(MARKER = "MUTATIONS", FIGURE = "HYPER")

  .write_sample_bed("S1", "Case", data.frame(CHR = c("chr1", "chr2"), START = c(100, 50), END = c(101, 51), value = c(1, 1)))
  # S3 found nothing and wrote no bed file: it must still get its column.
  SEMseeker:::anno_create_position_pivots(data.frame(Sample_ID = c("S1", "S3"), Sample_Group = "Case"), keys)

  pivot_path <- SEMseeker:::io_pivot_file_name_parquet("MUTATIONS", "HYPER", "POSITION", "WHOLE")
  first_pivot <- as.data.frame(nanoparquet::read_parquet(pivot_path))
  expect_identical(colnames(first_pivot), c("CHR", "START", "END", "S1", "S3"))
  expect_identical(first_pivot$CHR, c("1", "2"))
  expect_equal(first_pivot$S3, c(0, 0))

  # S2 shares chr1:100 with S1 and adds chr1:200, a position the pivot does
  # not have yet; S4 found nothing.
  .write_sample_bed("S2", "Case", data.frame(CHR = c("chr1", "chr1"), START = c(100, 200), END = c(101, 201), value = c(2, 3)))
  SEMseeker:::anno_create_position_pivots(data.frame(Sample_ID = c("S1", "S2", "S3", "S4"), Sample_Group = "Case"), keys)

  second_pivot <- as.data.frame(nanoparquet::read_parquet(pivot_path))
  expect_identical(colnames(second_pivot), c("CHR", "START", "END", "S1", "S3", "S2", "S4"))
  # the position found only by S2 keeps its coordinates
  expect_identical(second_pivot$CHR, c("1", "1", "2"))
  expect_equal(as.numeric(second_pivot$START), c(100, 200, 50))
  expect_equal(as.numeric(second_pivot$END), c(101, 201, 51))
  expect_equal(second_pivot$S1, c(1, NA, 1))
  expect_equal(second_pivot$S3, c(0, NA, 0))
  expect_equal(second_pivot$S2, c(2, 3, NA))
  expect_equal(second_pivot$S4, c(0, 0, 0))
  expect_false(file.exists(paste0(pivot_path, ".tmp")))
})
