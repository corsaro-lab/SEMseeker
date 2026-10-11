# sem_marker_value_per_probe_plot(), which through 0.99.6 was
# sem_manhattan_plot_marker_per_sample(): exported, documented, and unable to
# run. Every test here writes real parquet pivots and calls the exported
# function; nothing on the path is stubbed.

.probe_session <- function() {
  folder <- sem_test_folder()
  dir.create(file.path(folder, "Data"), recursive = TRUE, showWarnings = FALSE)
  SEMseeker:::core_init_env(folder, parallel_strategy = "sequential",
                            iqrTimes = 3, verbosity = 1)
  folder
}

.write_probe_pivot <- function(figure, probes, ...) {
  path <- SEMseeker:::io_pivot_file_name_parquet("DELTAS", figure,
                                               "PROBE", "WHOLE")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  polars::as_polars_df(data.frame(AREA = probes, ..., stringsAsFactors = FALSE)
                       )$write_parquet(path)
  path
}

# One sample hyper at cg1, hypo at cg3, nothing at cg2.
.write_three <- function() {
  .write_probe_pivot("BOTH",  c("cg1", "cg2", "cg3"), CASE_1 = c(4, 0, 7),
                     REF_1 = c(0, 0, 0))
  .write_probe_pivot("HYPER", c("cg1", "cg2", "cg3"), CASE_1 = c(4, 0, 0),
                     REF_1 = c(0, 0, 0))
  .write_probe_pivot("HYPO",  c("cg1", "cg2", "cg3"), CASE_1 = c(0, 0, 7),
                     REF_1 = c(0, 0, 0))
}

test_that("the chart is written, and the coordinates are in its name", {
  skip_on_cran()
  folder <- .probe_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)
  .write_three()

  written <- SEMseeker::sem_marker_value_per_probe_plot("DELTAS", "BOTH", "CASE_1")

  expect_false(is.null(written))
  expect_true(file.exists(written))
  expect_gt(file.size(written), 1000)
  expect_match(basename(written), "DELTAS")
  expect_match(basename(written), "BOTH")
  expect_match(basename(written), "CASE_1")
  # it lands under its own folder, not mixed with the per-area charts.
  expect_match(dirname(written), "MARKER_VALUE_PER_PROBE")
})

test_that("the direction of each probe comes from HYPER and HYPO", {
  skip_on_cran()
  folder <- .probe_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)
  .write_three()

  # Reach the classification through the helper the door uses, so the assertion
  # is on the data that reaches the drawing rather than on a PNG.
  hyper <- SEMseeker:::.sem_probe_column_get("DELTAS", "HYPER", "CASE_1")
  hypo  <- SEMseeker:::.sem_probe_column_get("DELTAS", "HYPO",  "CASE_1")
  expect_equal(hyper$AREA[hyper$VALUE != 0], "cg1")
  expect_equal(hypo$AREA[hypo$VALUE != 0], "cg3")

  both <- SEMseeker:::.sem_probe_column_get("DELTAS", "BOTH", "CASE_1")
  expect_equal(both$VALUE, c(4, 0, 7))
  # the row order of the pivot is the order of the axis, so it is preserved.
  expect_equal(both$AREA, c("cg1", "cg2", "cg3"))
})

test_that("a sample with no events still produces a chart", {
  skip_on_cran()
  folder <- .probe_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)
  .write_three()

  # REF_1 is zero everywhere. The extent examined is still the whole axis: a
  # chart of a sample with no events is information, not an absence of data.
  written <- SEMseeker::sem_marker_value_per_probe_plot("DELTAS", "BOTH", "REF_1")
  expect_true(file.exists(written))
})

test_that("probes_selection names the probes instead of indexing a row order", {
  skip_on_cran()
  folder <- .probe_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)
  .write_three()

  kept <- SEMseeker:::.sem_probe_column_get("DELTAS", "BOTH", "CASE_1",
                                           probes_selection = c("cg1", "cg3"))
  expect_equal(kept$AREA, c("cg1", "cg3"))
  expect_equal(kept$VALUE, c(4, 7))

  # a selection naming nothing present is nothing to draw, not an error.
  expect_null(SEMseeker:::.sem_probe_column_get("DELTAS", "BOTH", "CASE_1",
                                               probes_selection = "cg_absent"))
})

test_that("a missing pivot and an unknown sample are both nothing to draw", {
  skip_on_cran()
  folder <- .probe_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  # no pivot at all
  expect_null(SEMseeker::sem_marker_value_per_probe_plot("DELTAS", "BOTH", "CASE_1"))

  .write_three()
  # a sample that is not a column of the pivot
  expect_null(SEMseeker::sem_marker_value_per_probe_plot("DELTAS", "BOTH", "NOT_A_SAMPLE"))
})

test_that("an existing chart is not redrawn unless asked", {
  skip_on_cran()
  folder <- .probe_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)
  .write_three()

  first <- SEMseeker::sem_marker_value_per_probe_plot("DELTAS", "BOTH", "CASE_1")
  stamp <- file.mtime(first)
  again <- SEMseeker::sem_marker_value_per_probe_plot("DELTAS", "BOTH", "CASE_1")
  expect_equal(again, first)
  expect_equal(file.mtime(again), stamp)

  redrawn <- SEMseeker::sem_marker_value_per_probe_plot("DELTAS", "BOTH", "CASE_1",
                                            overwrite = TRUE)
  expect_equal(redrawn, first)
  expect_gte(as.numeric(file.mtime(redrawn)), as.numeric(stamp))
})

test_that("the outlier palette names exactly the four classes the drawing knows", {
  skip_on_cran()
  folder <- .probe_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  palette <- SEMseeker:::.sem_outlier_palette(SEMseeker:::core_get_session_info())
  expect_setequal(names(palette),
                  c("Hyper", "Hypo", "Non-outlier", "Reference"))
  expect_true(all(nzchar(palette)))
  # the colours come from the session, so a study's charts agree with each other.
  expect_equal(unname(palette[["Hyper"]]),
               SEMseeker:::core_get_session_info()$color_palette_darker[1])
})
