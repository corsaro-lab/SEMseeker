# Reopening the results of a previous run with result_folder.
#
# The charts and readers used to take only the session left in memory by
# semseeker(), so in a new R session they stopped with "called without result
# folder!". Each case below closes the session that wrote the fixture, which
# empties it, before calling the function: that is the new R session.
#
# Expected values come from the fixture: where the files must land (the folder
# named, never the one the run was written in), the plot format the run chose,
# and the probe totals written by hand.

.results_fixture <- function(plot_format = "pdf") {
  folder <- sem_test_folder()
  dir.create(file.path(folder, "Data"), recursive = TRUE, showWarnings = FALSE)
  SEMseeker:::core_init_env(folder, parallel_strategy = "sequential", iqrTimes = 3,
                            verbosity = 1, plot_format = plot_format)
  write_pivot <- function(marker, figure, area, frame, ...) {
    path <- SEMseeker:::io_pivot_file_name_parquet(marker, figure, area, "WHOLE", ...)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    SEMseeker:::io_write_frame_parquet(frame, path)
  }
  # Probe totals of DELTAS BOTH: CG1 = 5 + 9 + 0 = 14, CG2 = 0 + 1 + 0 = 1.
  write_pivot("DELTAS", "BOTH", "PROBE", data.frame(AREA = c("CG1", "CG2"), CASE_1 = c(5, 0), CASE_2 = c(9, 1), REF_1 = c(0, 0)))
  write_pivot("DELTAS", "HYPER", "PROBE", data.frame(AREA = c("CG1", "CG2"), CASE_1 = c(5, 0), CASE_2 = c(0, 1), REF_1 = c(0, 0)))
  write_pivot("DELTAS", "HYPO", "PROBE", data.frame(AREA = c("CG1", "CG2"), CASE_1 = c(0, 0), CASE_2 = c(9, 0), REF_1 = c(0, 0)))
  write_pivot("MUTATIONS", "HYPER", "GENE", data.frame(AREA = c("BRCA1", "TP53"), CASE_1 = c(2, 0), CASE_2 = c(1, 3), REF_1 = c(0, 1)),
              aggregation = "SUM")
  write_pivot("SIGNAL", SEMseeker:::io_signal_figure(), "PROBE", data.frame(AREA = c("CG1", "CG2"), CASE_1 = c(0.95, 0.5), CASE_2 = c(0.1, 0.5), REF_1 = c(0.52, 0.5)))
  SEMseeker:::io_write_frame_parquet(data.frame(PROBE = "CG1", signal_inferior_thresholds = 0.2, signal_superior_thresholds = 0.8,
                                                signal_median_values = 0.5, iqr = 0.2, q1 = 0.4, q3 = 0.6),
                                     SEMseeker:::io_file_path_build(file.path(folder, "Data"), "1_signal_thresholds", "parquet"))
  SEMseeker:::core_close_env()
  folder
}

.results_sheet <- function() data.frame(Sample_ID = c("CASE_1", "CASE_2", "REF_1"), Sample_Group = c("Case", "Case", "Reference"),
                                        stringsAsFactors = FALSE)

.results_inference <- function() data.frame(independent_variable = "Sample_Group", family_test = "wilcoxon", transformation_y = "none",
                                            aggregation = "MEAN", scope = "INSTANCE", stringsAsFactors = FALSE)

# The five charts that draw from the fixture, each called on a named folder.
.results_charts <- list(
  heatmap    = function(folder) SEMseeker::sem_marker_heatmap_plot("MUTATIONS", "HYPER", "GENE", "WHOLE", "Sample_Group", .results_sheet(),
                                                                   aggregation = "SUM", result_folder = folder),
  per_area   = function(folder) SEMseeker::sem_marker_value_per_area_plot("MUTATIONS", "HYPER", "GENE", "WHOLE", "Sample_Group", .results_sheet(),
                                                                          aggregation = "SUM", result_folder = folder),
  per_probe  = function(folder) SEMseeker::sem_marker_value_per_probe_plot("DELTAS", "BOTH", "CASE_1", result_folder = folder),
  per_sample = function(folder) SEMseeker::sem_marker_value_per_sample_plot("DELTAS", "BOTH", "CG1", .results_sheet(), result_folder = folder),
  threshold  = function(folder) SEMseeker::sem_signal_threshold_per_sample_plot("CG1", .results_sheet(), result_folder = folder)
)

# All eight, called on a folder that is not a result folder.
.results_all <- c(.results_charts, list(
  probe_select = function(folder) SEMseeker::sem_probe_select_by_statistic("DELTAS", "BOTH", result_folder = folder),
  lollipop     = function(folder) SEMseeker::enrich_term_lollipop_plot("pathfindR", .results_inference(), result_folder = folder),
  circos       = function(folder) SEMseeker::enrich_circos_plot("pathfindR", .results_inference(), result_folder = folder)
))

.chart_files <- function(folder) list.files(file.path(folder, "Chart"), recursive = TRUE)

for (chart_name in names(.results_charts)) {
  test_that(paste0(chart_name, ": a new session draws into the folder named, in the run's format"), {
    skip_on_cran()
    original <- .results_fixture(plot_format = "pdf")
    copy_parent <- sem_test_folder()
    on.exit({ SEMseeker:::core_close_env(); unlink(c(original, copy_parent), recursive = TRUE) }, add = TRUE)
    expect_true(file.copy(original, copy_parent, recursive = TRUE))
    copy <- file.path(copy_parent, basename(original))

    written <- .results_charts[[chart_name]](copy)

    expect_true(file.exists(written))
    expect_true(startsWith(normalizePath(written), normalizePath(file.path(copy, "Chart"))))
    expect_identical(tools::file_ext(written), "pdf")
    expect_length(.chart_files(original), 0L)
    # The session opened for the chart is closed again.
    expect_error(SEMseeker:::core_get_session_info(), "without result folder")
  })
}

test_that("probe_select: a new session reads the pivots of the folder named", {
  skip_on_cran()
  folder <- .results_fixture()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) }, add = TRUE)

  picked <- SEMseeker::sem_probe_select_by_statistic("DELTAS", "BOTH", c("MIN", "MAX"), result_folder = folder)

  expect_identical(picked$PROBE, c("CG2", "CG1"))
  expect_equal(picked$VALUE, c(1, 14))
})

test_that("the enrichment charts open the folder named and find no report there", {
  skip_on_cran()
  folder <- .results_fixture()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) }, add = TRUE)

  expect_length(SEMseeker::enrich_term_lollipop_plot("pathfindR", .results_inference(), result_folder = folder), 0L)
  expect_length(SEMseeker::enrich_circos_plot("pathfindR", .results_inference(), result_folder = folder), 0L)
  expect_error(SEMseeker:::core_get_session_info(), "without result folder")
})

for (function_name in names(.results_all)) {
  test_that(paste0(function_name, ": a folder SEMseeker did not write is refused, and left as it was"), {
    skip_on_cran()
    not_results <- sem_test_folder()
    on.exit({ SEMseeker:::core_close_env(); unlink(not_results, recursive = TRUE) }, add = TRUE)
    dir.create(not_results, recursive = TRUE, showWarnings = FALSE)

    expect_error(.results_all[[function_name]](not_results), "not a SEMseeker result folder")
    expect_length(list.files(not_results, recursive = TRUE, include.dirs = TRUE), 0L)
  })
}

test_that("without result_folder the open session is used, and naming it does not reopen it", {
  skip_on_cran()
  folder <- .results_fixture(plot_format = "png")
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) }, add = TRUE)
  SEMseeker:::core_init_env(folder, parallel_strategy = "sequential", verbosity = 1)
  session_id <- SEMseeker:::core_get_session_info()$session_id

  without_argument <- SEMseeker::sem_marker_value_per_sample_plot("DELTAS", "BOTH", "CG1", .results_sheet())
  expect_true(startsWith(normalizePath(without_argument), normalizePath(file.path(folder, "Chart"))))
  expect_identical(tools::file_ext(without_argument), "png")

  SEMseeker::sem_marker_value_per_sample_plot("DELTAS", "BOTH", "CG1", .results_sheet(), overwrite = TRUE, result_folder = folder)
  # Still the same session, still open.
  expect_identical(SEMseeker:::core_get_session_info()$session_id, session_id)
})
