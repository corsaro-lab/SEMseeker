# Results are indexed by coordinates, never by the probe ids of a manufacturer.
#
# Through 0.99.6 the marker pivots at PROBE level were keyed by CHR_START while
# the SIGNAL pivot and the thresholds kept the Illumina probe id: a probe chosen
# on one could not be looked up on the other. Everything is now keyed by
# position, and a probe id given to a chart is refused.
#
# One real run on the probes of two imprinting DMRs of the bundled GSE133774
# fixture, area PROBE, shared by the tests of this file. Expected values are
# the shape of a position id, the equality of two files the run wrote, and the
# absence of probe ids; none of them is computed by the code under test.

.positions_run <- local({
  result_folder <- NULL
  function() {
    if (!is.null(result_folder) && dir.exists(result_folder))
      return(result_folder)
    utils::data("test_signal_gse133774", "test_samplesheet_gse133774", "test_master_features",
                package = "SEMseeker", envir = environment())
    demo_probes <- subset(test_master_features, DMR_LABEL %in% c("KCNQ1OT1:TSS-DMR", "H19_IGF2:IG-DMR") &
                            PROBE %in% rownames(test_signal_gse133774))
    result_folder <<- sem_test_folder()
    SEMseeker::semseeker(input = test_signal_gse133774[demo_probes$PROBE, , drop = FALSE],
                         sample_sheet = test_samplesheet_gse133774, result_folder = result_folder,
                         parallel_strategy = "sequential", areas = "PROBE")
    result_folder
  }
})

.position_pattern <- "^[^_]+_[0-9]+$"

.probe_level_pivots <- function(result_folder) {
  list.files(file.path(result_folder, "Data", "Pivots"), pattern = "_INSTANCE_PROBE_WHOLE_VALUE_.*\\.parquet$",
             recursive = TRUE, full.names = TRUE)
}

test_that("every PROBE pivot, SIGNAL included, is keyed by the same positions", {
  skip_on_cran()
  skip_if_not_installed("nanoparquet")
  pivots <- .probe_level_pivots(.positions_run())
  expect_true(any(grepl("/SIGNAL/", pivots)))
  expect_gt(length(pivots), 2L)

  areas <- lapply(pivots, function(pivot) sort(as.character(nanoparquet::read_parquet(pivot)$AREA)))
  for (i in seq_along(areas))
    expect_true(all(grepl(.position_pattern, areas[[i]])), info = basename(pivots[i]))
  for (i in seq_along(areas))
    expect_identical(areas[[i]], areas[[1]], info = basename(pivots[i]))
})

test_that("the SIGNAL PROBE pivot holds the values of the POSITION pivot", {
  skip_on_cran()
  skip_if_not_installed("nanoparquet")
  signal_folder <- file.path(.positions_run(), "Data", "Pivots", "SIGNAL")
  by_probe <- as.data.frame(nanoparquet::read_parquet(list.files(signal_folder, "_INSTANCE_PROBE_WHOLE_VALUE_.*\\.parquet$", full.names = TRUE)))
  by_position <- as.data.frame(nanoparquet::read_parquet(list.files(signal_folder, "_INSTANCE_POSITION_WHOLE_VALUE_.*\\.parquet$", full.names = TRUE)))
  rownames(by_position) <- paste0(by_position$CHR, "_", by_position$START)
  sample_columns <- setdiff(colnames(by_position), c("CHR", "START", "END"))

  expect_setequal(by_probe$AREA, rownames(by_position))
  expect_equal(as.matrix(by_probe[, sample_columns]), as.matrix(by_position[by_probe$AREA, sample_columns]), ignore_attr = TRUE)
})

test_that("a probe chosen on a marker pivot finds its thresholds", {
  skip_on_cran()
  result_folder <- .positions_run()
  on.exit(SEMseeker:::core_close_env(), add = TRUE)
  sheet <- utils::read.csv2(file.path(result_folder, "Data", "1_SAMPLE_SHEET_ORIGINAL.csv"))

  picked <- SEMseeker::sem_probe_select_by_statistic("DELTAS", "HYPER", "MAX", result_folder = result_folder)
  written <- SEMseeker::sem_signal_threshold_per_sample_plot(picked$PROBE, sheet, result_folder = result_folder)

  expect_true(grepl(.position_pattern, picked$PROBE))
  expect_false(is.null(written))
  expect_true(file.exists(written))
})

test_that("a manufacturer's probe id is refused before anything is read", {
  sheet <- data.frame(Sample_ID = "CASE_1", Sample_Group = "Case")
  expect_error(SEMseeker::sem_marker_value_per_sample_plot("DELTAS", "HYPER", "cg11680158", sheet), "not a position: cg11680158")
  expect_error(SEMseeker::sem_signal_threshold_per_sample_plot("cg11680158", sheet), "not a position: cg11680158")
  expect_error(SEMseeker::sem_marker_value_per_probe_plot("DELTAS", "HYPER", "CASE_1", probes_selection = c("11_2721619", "cg00000029_TC21")),
               "not a position: cg00000029_TC21")
})

test_that("nothing the run writes under Data carries a probe id", {
  skip_on_cran()
  skip_if_not_installed("nanoparquet")
  data_folder <- file.path(.positions_run(), "Data")
  probe_id_pattern <- "^cg[0-9]{8}"

  for (parquet in list.files(data_folder, "\\.parquet$", recursive = TRUE, full.names = TRUE)) {
    frame <- as.data.frame(nanoparquet::read_parquet(parquet))
    text_values <- unlist(lapply(frame[vapply(frame, is.character, logical(1))], unique), use.names = FALSE)
    expect_false(any(grepl(probe_id_pattern, c(text_values, colnames(frame)))), info = basename(parquet))
  }
  for (text_file in list.files(data_folder, "\\.(csv|json)$", recursive = TRUE, full.names = TRUE))
    expect_false(any(grepl("cg[0-9]{8}", readLines(text_file, warn = FALSE))), info = basename(text_file))
})
