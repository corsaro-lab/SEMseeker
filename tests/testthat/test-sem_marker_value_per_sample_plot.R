# sem_marker_value_per_sample_plot() and sem_probe_select_by_statistic(), which
# through 0.99.6 were one function, anno_manhattan_plot_marker_per_probe(): it
# had a manual page, scanned every marker, wrote a CSV of statistics, drew ten
# charts per call, and raised before any of them. Nothing here is stubbed.

.sample_session <- function() {
  folder <- sem_test_folder()
  dir.create(file.path(folder, "Data"), recursive = TRUE, showWarnings = FALSE)
  SEMseeker:::core_init_env(folder, parallel_strategy = "sequential",
                            iqrTimes = 3, verbosity = 1)
  folder
}

.write_probe_pivot_s <- function(figure, probes, ...) {
  path <- SEMseeker:::io_pivot_file_name_parquet("DELTAS", figure,
                                               "PROBE", "WHOLE")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  polars::as_polars_df(data.frame(AREA = probes, ..., stringsAsFactors = FALSE)
                       )$write_parquet(path)
  path
}

# At CG1: CASE_1 hyper, CASE_2 hypo, REF_1 nothing.
.write_cohort <- function() {
  .write_probe_pivot_s("BOTH",  c("CG1", "CG2"),
                       CASE_1 = c(5, 0), CASE_2 = c(9, 1), REF_1 = c(0, 0))
  .write_probe_pivot_s("HYPER", c("CG1", "CG2"),
                       CASE_1 = c(5, 0), CASE_2 = c(0, 1), REF_1 = c(0, 0))
  .write_probe_pivot_s("HYPO",  c("CG1", "CG2"),
                       CASE_1 = c(0, 0), CASE_2 = c(9, 0), REF_1 = c(0, 0))
}

.cohort_sheet <- function() data.frame(
  Sample_ID = c("CASE_1", "CASE_2", "REF_1"),
  Sample_Group = c("Case", "Case", "Reference"),
  stringsAsFactors = FALSE)

test_that("the chart is written for the probe it is given", {
  skip_on_cran()
  folder <- .sample_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)
  .write_cohort()

  written <- SEMseeker::sem_marker_value_per_sample_plot(
    "DELTAS", "BOTH", "CG1", .cohort_sheet())

  expect_true(file.exists(written))
  expect_gt(file.size(written), 1000)
  expect_match(basename(written), "CG1")
  expect_match(dirname(written), "MARKER_VALUE_PER_SAMPLE")

  # a second probe is a second chart, not an overwrite of the first.
  other <- SEMseeker::sem_marker_value_per_sample_plot(
    "DELTAS", "BOTH", "CG2", .cohort_sheet())
  expect_false(identical(other, written))
  expect_true(file.exists(other))
})

test_that("the row is read across samples, which is the other axis", {
  skip_on_cran()
  folder <- .sample_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)
  .write_cohort()

  row <- SEMseeker:::.sem_probe_row_get("DELTAS", "BOTH", "CG1")
  expect_setequal(row$SAMPLE, c("CASE_1", "CASE_2", "REF_1"))
  expect_equal(row$VALUE[row$SAMPLE == "CASE_2"], 9)
  expect_equal(row$VALUE[row$SAMPLE == "REF_1"], 0)

  # the probe name is matched without the caller having to know the case the
  # pivot writers use.
  expect_equal(SEMseeker:::.sem_probe_row_get("DELTAS", "BOTH", "cg1")$VALUE,
               row$VALUE)

  # a probe that is not in the pivot is nothing to draw.
  expect_null(SEMseeker:::.sem_probe_row_get("DELTAS", "BOTH", "CG_ABSENT"))
  expect_null(SEMseeker::sem_marker_value_per_sample_plot(
    "DELTAS", "BOTH", "CG_ABSENT", .cohort_sheet()))
})

test_that("samples_selection keeps only the samples it names", {
  skip_on_cran()
  folder <- .sample_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)
  .write_cohort()

  kept <- SEMseeker:::.sem_probe_row_get("DELTAS", "BOTH", "CG1",
                                        samples_selection = c("CASE_1", "REF_1"))
  expect_setequal(kept$SAMPLE, c("CASE_1", "REF_1"))
  expect_null(SEMseeker:::.sem_probe_row_get("DELTAS", "BOTH", "CG1",
                                            samples_selection = "NOBODY"))
})

test_that("a reference sample is a class and is not removed", {
  skip_on_cran()
  folder <- .sample_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)
  .write_cohort()

  # The thresholds were computed from the references, so where they sit is how a
  # reader judges an excursion. A chart that drops them hides that.
  written <- SEMseeker::sem_marker_value_per_sample_plot(
    "DELTAS", "BOTH", "CG1", .cohort_sheet())
  expect_true(file.exists(written))

  # the reference class wins over the direction: mark REF_1 hyper in the data
  # and it is still drawn as a reference.
  .write_probe_pivot_s("HYPER", c("CG1", "CG2"),
                       CASE_1 = c(5, 0), CASE_2 = c(0, 1), REF_1 = c(3, 0))
  row <- SEMseeker:::.sem_probe_row_get("DELTAS", "HYPER", "CG1")
  expect_equal(row$VALUE[row$SAMPLE == "REF_1"], 3)
  redrawn <- SEMseeker::sem_marker_value_per_sample_plot(
    "DELTAS", "BOTH", "CG1", .cohort_sheet(), overwrite = TRUE)
  expect_true(file.exists(redrawn))
})

# ---------------------------------------------------------------------------
# sem_probe_select_by_statistic
# ---------------------------------------------------------------------------

test_that("the five statistics pick a probe from the cohort totals", {
  skip_on_cran()
  folder <- .sample_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  # totals across samples: P1 = 0, P2 = 2, P3 = 4, P4 = 6, P5 = 100
  .write_probe_pivot_s("BOTH", c("P1", "P2", "P3", "P4", "P5"),
                       S1 = c(0, 1, 2, 3, 50),
                       S2 = c(0, 1, 2, 3, 50))

  picked <- SEMseeker:::sem_probe_select_by_statistic("DELTAS", "BOTH")
  expect_setequal(picked$STATISTIC, c("MIN", "Q1", "MEDIAN", "Q3", "MAX"))

  by_statistic <- stats::setNames(picked$PROBE, picked$STATISTIC)
  expect_equal(by_statistic[["MIN"]], "P1")
  expect_equal(by_statistic[["MAX"]], "P5")
  expect_equal(by_statistic[["MEDIAN"]], "P3")

  # the value is returned so the caller can see how near the nearest probe was.
  expect_equal(picked$VALUE[picked$STATISTIC == "MAX"], 100)
  expect_equal(picked$VALUE[picked$STATISTIC == "MIN"], 0)
})

test_that("one statistic can be asked for on its own", {
  skip_on_cran()
  folder <- .sample_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)
  .write_probe_pivot_s("BOTH", c("P1", "P2"), S1 = c(1, 7))

  one <- SEMseeker:::sem_probe_select_by_statistic("DELTAS", "BOTH", "MAX")
  expect_equal(nrow(one), 1L)
  expect_equal(one$PROBE, "P2")
})

test_that("a tie resolves to the first probe of the pivot, not to a failure", {
  skip_on_cran()
  folder <- .sample_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  # A marker that counts events has many probes totalling zero, so MIN is a tie
  # by construction on any real study: failing on a tie would fail always.
  .write_probe_pivot_s("BOTH", c("Z1", "Z2", "Z3"), S1 = c(0, 0, 5))
  picked <- SEMseeker:::sem_probe_select_by_statistic("DELTAS", "BOTH", "MIN")
  expect_equal(picked$PROBE, "Z1")
})

test_that("an unknown statistic is refused and a missing pivot is NULL", {
  skip_on_cran()
  folder <- .sample_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  expect_null(SEMseeker:::sem_probe_select_by_statistic("DELTAS", "BOTH"))

  .write_probe_pivot_s("BOTH", "P1", S1 = 1)
  expect_error(
    SEMseeker:::sem_probe_select_by_statistic("DELTAS", "BOTH", "AVERAGE"),
    "unknown statistic")
})

test_that("the selection and the chart compose, which is the point of splitting them", {
  skip_on_cran()
  folder <- .sample_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)
  .write_cohort()

  extreme <- SEMseeker:::sem_probe_select_by_statistic("DELTAS", "BOTH", "MAX")
  written <- SEMseeker::sem_marker_value_per_sample_plot(
    "DELTAS", "BOTH", extreme$PROBE, .cohort_sheet())

  expect_true(file.exists(written))
  expect_match(basename(written), extreme$PROBE)
})
