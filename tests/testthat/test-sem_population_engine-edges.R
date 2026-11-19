# Edge cases of the population markers computed in the data engine: a missing
# signal, a constant reference, and values that fall exactly on a bin edge.
# The expected values are worked out by hand from the formulas in the header
# of sem_analyze_population_bulk() and sem_deltaX_get(), not by running them.

#' Open a session in a fresh folder for one test.
#'
#' @param folder_prefix prefix of the temporary result folder
#' @return the folder, removed with the session when the calling test ends
.population_edges_session <- function(folder_prefix) {
  temp_folder <- tempfile(folder_prefix)
  dir.create(temp_folder, recursive = TRUE, showWarnings = FALSE)
  SEMseeker:::core_init_env(temp_folder, parallel_strategy = "sequential", bulk_population = TRUE,
                            LESIONS_BP = 5000L, bonferroni_threshold = 0.05)
  withr::defer({ SEMseeker:::core_close_env(); unlink(temp_folder, recursive = TRUE) }, envir = parent.frame())
  temp_folder
}

#' Read a POSITION pivot with a reader that is not the engine under test.
#'
#' @param marker,figure the pivot to read
#' @return a data.frame sorted on CHR and START
.read_position_pivot <- function(marker, figure) {
  pivot <- as.data.frame(nanoparquet::read_parquet(SEMseeker:::io_pivot_file_name_parquet(marker, figure, "POSITION", "WHOLE")))
  pivot[order(pivot$CHR, pivot$START), , drop = FALSE]
}

test_that("a missing signal stays missing in DELTAS, MUTATIONS and DELTAR", {
  skip_if_not_installed("nanoparquet")
  .population_edges_session("population_missing_")

  # S1 is above the high threshold, missing, and below the low one; the
  # third position has a constant reference (high == low).
  signal_pivot <- data.frame(CHR = "1", START = c(100L, 200L, 300L), END = c(101L, 201L, 301L),
                             S1 = c(0.9, NA, 0.1), stringsAsFactors = FALSE)
  signal_path <- SEMseeker:::io_pivot_file_name_parquet("SIGNAL", SEMseeker:::io_signal_figure(), "POSITION", "WHOLE")
  dir.create(dirname(signal_path), recursive = TRUE, showWarnings = FALSE)
  nanoparquet::write_parquet(signal_pivot, signal_path)
  thresholds <- data.frame(CHR = "1", START = signal_pivot$START, END = signal_pivot$END,
                           signal_superior_thresholds = c(0.8, 0.8, 0.5),
                           signal_inferior_thresholds = c(0.2, 0.2, 0.5),
                           stringsAsFactors = FALSE)

  # LESIONS is out of scope here, and sem_lesions_get_bulk() stops on a
  # missing MUTATIONS value (the pipeline imputes the signal before this
  # step): existing LESIONS pivots make the function skip it.
  for (figure in c("HYPER", "HYPO")) {
    lesions_path <- SEMseeker:::io_pivot_file_name_parquet("LESIONS", figure, "POSITION", "WHOLE")
    dir.create(dirname(lesions_path), recursive = TRUE, showWarnings = FALSE)
    nanoparquet::write_parquet(signal_pivot[, c("CHR", "START", "END")], lesions_path)
  }

  SEMseeker:::sem_analyze_population_bulk(signal_data = NULL, sample_sheet = NULL,
                                          signal_thresholds = thresholds, probe_features = NULL)

  # DELTAS_HYPER = max(S1 - high, 0): 0.1, NA, 0
  expect_equal(.read_position_pivot("DELTAS", "HYPER")$S1, c(0.9 - 0.8, NA, 0))
  # DELTAS_HYPO = max(low - S1, 0): 0, NA, 0.4
  expect_equal(.read_position_pivot("DELTAS", "HYPO")$S1, c(0, NA, 0.5 - 0.1))
  expect_identical(.read_position_pivot("MUTATIONS", "HYPER")$S1, c(1L, NA, 0L))
  expect_identical(.read_position_pivot("MUTATIONS", "HYPO")$S1, c(0L, NA, 1L))
  # DELTAR = DELTAS / max(high - low, 1e-9): the constant reference divides
  # by 1e-9 rather than by zero
  expect_equal(.read_position_pivot("DELTAR", "HYPER")$S1, c((0.9 - 0.8) / 0.6, NA, 0))
  expect_equal(.read_position_pivot("DELTAR", "HYPO")$S1, c(0, NA, (0.5 - 0.1) / 1e-9))
})

test_that("DELTAP bins are closed on the right and leave zero and missing values empty", {
  skip_if_not_installed("nanoparquet")
  .population_edges_session("population_cut_")

  # Non-zero values over HYPER and HYPO: 1, 2, 3, 5. With Q = 4 the breaks
  # are 1, 2, 3, 4, 5 and the bins (-inf, 2], (2, 3], (3, 4], (4, inf).
  coordinates <- data.frame(CHR = "1", START = seq(100L, 500L, by = 100L), END = seq(101L, 501L, by = 100L),
                            stringsAsFactors = FALSE)
  for (figure in c("HYPER", "HYPO")) {
    deltas <- coordinates
    deltas$S1 <- if (figure == "HYPER") c(1, 2, 0, NA, 3) else c(5, 0, 0, 0, 0)
    deltas_path <- SEMseeker:::io_pivot_file_name_parquet("DELTAS", figure, "POSITION", "WHOLE")
    dir.create(dirname(deltas_path), recursive = TRUE, showWarnings = FALSE)
    nanoparquet::write_parquet(deltas, deltas_path)
  }

  SEMseeker:::sem_deltaX_get(markers = "DELTAP")

  # 1 and 2 (on the edge) fall in the first bin, 3 (on the edge) in the
  # second, 5 in the last; 0 and NA give no label
  expect_identical(.read_position_pivot("DELTAP", "HYPER")$S1, c(1L, 1L, NA, NA, 2L))
  expect_identical(.read_position_pivot("DELTAP", "HYPO")$S1, c(4L, NA, NA, NA, NA))
})
