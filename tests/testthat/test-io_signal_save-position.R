# The POSITION pivot of the signal is written by one sort, or, when that sort
# runs out of memory, one chromosome at a time. Both must give the same file,
# in natural chromosome order. The out-of-memory error is simulated on the
# single query only: the fallback and the dispatcher run for real. A real
# memory limit was not used because the window in which the single sort
# fails and the per-chromosome one succeeds is narrow and depends on the
# width of the file (measured in the engine migration notes).

#' Six probes on four chromosomes, out of order, with two samples.
#'
#' @return list with signal_data (probes in the rows, ids as row names),
#'   sample_sheet and probe_features (PROBE, CHR, START, END)
.signal_position_fixture <- function() {
  probe_ids <- sprintf("p%d", 1:6)
  probe_features <- data.frame(PROBE = probe_ids, CHR = c("10", "2", "X", "2", "1", "10"),
                               START = c(500L, 300L, 100L, 100L, 700L, 200L), stringsAsFactors = FALSE)
  probe_features$END <- probe_features$START + 1L
  signal_data <- data.frame(S1 = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), S2 = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4),
                            row.names = probe_ids)
  sample_sheet <- data.frame(Sample_ID = c("S1", "S2"), Sample_Group = "Case", stringsAsFactors = FALSE)
  list(signal_data = signal_data, sample_sheet = sample_sheet, probe_features = probe_features)
}

#' Save the fixture in a fresh session and read the POSITION pivot back.
#'
#' @return the POSITION pivot as a data.frame, in file order, read with a
#'   reader that is not the engine under test
.save_signal_position <- function() {
  temp_folder <- tempfile("signal_position_")
  dir.create(temp_folder, recursive = TRUE, showWarnings = FALSE)
  SEMseeker:::core_init_env(temp_folder, parallel_strategy = "sequential", tech = "K850")
  on.exit({ SEMseeker:::core_close_env(); unlink(temp_folder, recursive = TRUE) }, add = TRUE)
  fixture <- .signal_position_fixture()
  SEMseeker:::io_signal_save(fixture$signal_data, fixture$sample_sheet, batch_id = 1L, probe_features = fixture$probe_features)
  position_path <- SEMseeker:::io_pivot_file_name_parquet("SIGNAL", SEMseeker:::io_signal_figure(), "POSITION", "WHOLE")
  as.data.frame(nanoparquet::read_parquet(position_path))
}

# Natural order: 1:700 (p5), 2:100 (p4), 2:300 (p2), 10:200 (p6), 10:500 (p1), X:100 (p3).
.expected_position <- data.frame(CHR = c("1", "2", "2", "10", "10", "X"),
                                 START = c(700L, 100L, 300L, 200L, 500L, 100L),
                                 S1 = c(0.5, 0.4, 0.2, 0.6, 0.1, 0.3),
                                 S2 = c(0.5, 0.6, 0.8, 0.4, 0.9, 0.7), stringsAsFactors = FALSE)

test_that("the single sort writes the POSITION pivot in natural chromosome order", {
  skip_if_not_installed("nanoparquet")
  position <- .save_signal_position()
  expect_identical(names(position), c("CHR", "START", "END", "S1", "S2"))
  expect_identical(as.character(position$CHR), .expected_position$CHR)
  expect_identical(as.integer(position$START), .expected_position$START)
  expect_identical(as.integer(position$END), .expected_position$START + 1L)
  expect_identical(position$S1, .expected_position$S1)
  expect_identical(position$S2, .expected_position$S2)
})

test_that("an out-of-memory sort falls back on one chromosome at a time with the same result", {
  skip_if_not_installed("nanoparquet")
  # The failed query leaves a partial file behind, as a real COPY can.
  local_mocked_bindings(io_signal_position_single = function(connection, coordinates_view, probe_path, position_path) {
    writeLines("partial", position_path)
    stop("Out of Memory Error: could not allocate block of size 256.0 KiB (simulated)")
  })
  position <- .save_signal_position()
  expect_identical(names(position), c("CHR", "START", "END", "S1", "S2"))
  expect_identical(as.character(position$CHR), .expected_position$CHR)
  expect_identical(as.integer(position$START), .expected_position$START)
  expect_identical(position$S1, .expected_position$S1)
  expect_identical(position$S2, .expected_position$S2)
})

test_that("an error other than out of memory stops without the fallback", {
  fallback_called <- FALSE
  local_mocked_bindings(
    io_signal_position_single = function(connection, coordinates_view, probe_path, position_path) stop("IO Error: disk full (simulated)"),
    io_signal_position_per_chromosome = function(...) fallback_called <<- TRUE
  )
  expect_error(.save_signal_position(), "disk full")
  expect_false(fallback_called)
})

test_that("a failed fallback leaves no partial POSITION pivot behind", {
  # A partial file would make the next run believe the signal is saved.
  local_mocked_bindings(
    io_signal_position_single = function(connection, coordinates_view, probe_path, position_path) {
      writeLines("partial", position_path)
      stop("Out of Memory Error: could not allocate block of size 256.0 KiB (simulated)")
    },
    io_signal_position_per_chromosome = function(...) stop("Out of Memory Error: the fallback failed too (simulated)")
  )
  temp_folder <- tempfile("signal_position_")
  dir.create(temp_folder, recursive = TRUE, showWarnings = FALSE)
  SEMseeker:::core_init_env(temp_folder, parallel_strategy = "sequential", tech = "K850")
  withr::defer({ SEMseeker:::core_close_env(); unlink(temp_folder, recursive = TRUE) })
  fixture <- .signal_position_fixture()
  expect_error(SEMseeker:::io_signal_save(fixture$signal_data, fixture$sample_sheet, batch_id = 1L, probe_features = fixture$probe_features),
               "the fallback failed too")
  expect_false(file.exists(SEMseeker:::io_pivot_file_name_parquet("SIGNAL", SEMseeker:::io_signal_figure(), "POSITION", "WHOLE")))
})
