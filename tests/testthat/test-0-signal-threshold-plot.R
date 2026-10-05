# sem_signal_threshold_per_sample_plot(). The chart this replaces had a delta
# label built with a two-argument ifelse: it raised nothing, returned NA, and
# drew nothing, so the only quantitative content of the chart had never
# appeared. These assertions go through ggplot_build(), which renders.

.threshold_session <- function() {
  folder <- sem_test_folder()
  dir.create(file.path(folder, "Data"), recursive = TRUE, showWarnings = FALSE)
  SEMseeker:::core_init_env(folder, parallel_strategy = "sequential",
                            iqrTimes = 3, verbosity = 1)
  folder
}

.envelope <- function() list(superior = 0.8, inferior = 0.2, median = 0.5,
                             q1 = 0.4, q3 = 0.6)

# CASE_1 above the envelope, CASE_2 below it, MID inside, REF_1 a reference.
.signal <- function() data.frame(
  SAMPLE = c("CASE_1", "CASE_2", "MID", "REF_1"),
  VALUE = c(0.95, 0.05, 0.5, 0.52),
  stringsAsFactors = FALSE)

.sheet4 <- function() data.frame(
  Sample_ID = c("CASE_1", "CASE_2", "MID", "REF_1"),
  Sample_Group = c("Case", "Case", "Case", "Reference"),
  stringsAsFactors = FALSE)

.build <- function(...) SEMseeker:::.sem_signal_threshold_plot_build(...)

.palette <- c(Hyper = "blue", Hypo = "red",
              "Non-outlier" = "grey", Reference = "cyan")

test_that("the delta labels are drawn, one per excursion", {
  built <- ggplot2::ggplot_build(
    .build(.signal(), .envelope(), .sheet4(), "CG1", .palette))

  labels <- unlist(lapply(built$data, function(layer) layer$label))
  labels <- labels[!is.na(labels)]
  deltas <- labels[grepl("δ", labels)]

  # two excursions, two deltas. MID is inside the envelope and REF_1 is a
  # reference, so neither is annotated.
  expect_equal(length(deltas), 2L)
  expect_true(all(nzchar(deltas)))
})

test_that("the delta is the distance outside the envelope, in the right direction", {
  built <- ggplot2::ggplot_build(
    .build(.signal(), .envelope(), .sheet4(), "CG1", .palette))
  labels <- unlist(lapply(built$data, function(layer) layer$label))
  deltas <- labels[!is.na(labels) & grepl("δ", labels)]

  # CASE_1 is 0.95 against an upper limit of 0.8, so 0.15.
  # CASE_2 is 0.05 against a lower limit of 0.2, so 0.15 as well: a distance
  # outside the envelope, not a signed difference, so both are positive.
  # The number is parsed out rather than matched as a string: the exponent of
  # 1.5e-01 carries a minus of its own, so looking for "-" proves nothing.
  numbers <- as.numeric(sub("^\\S+\\s+", "", deltas))
  expect_equal(sort(numbers), c(0.15, 0.15))
  expect_true(all(numbers > 0))
})

test_that("the segment runs from the crossed threshold, so its length is the excursion", {
  built <- ggplot2::ggplot_build(
    .build(.signal(), .envelope(), .sheet4(), "CG1", .palette))
  segments <- built$data[[1]]

  # the hyper sample starts at the upper limit, the hypo one at the lower.
  expect_true(0.8 %in% segments$y)
  expect_true(0.2 %in% segments$y)
  # a sample inside the envelope gets a zero-length segment rather than a line
  # down to the axis: it has no excursion to show.
  inside <- segments[segments$y == segments$yend, , drop = FALSE]
  expect_gte(nrow(inside), 2L)
  expect_false(0 %in% segments$y)
})

test_that("the five envelope lines render, and their labels can be turned off", {
  with_labels <- ggplot2::ggplot_build(
    .build(.signal(), .envelope(), .sheet4(), "CG1", .palette,
           show_labels = TRUE))
  yintercepts <- unlist(lapply(with_labels$data,
                               function(layer) layer$yintercept))
  expect_setequal(yintercepts[!is.na(yintercepts)],
                  c(0.5, 0.8, 0.2, 0.4, 0.6))

  labels <- unlist(lapply(with_labels$data, function(layer) layer$label))
  labels <- labels[!is.na(labels)]
  expect_true("Upper Limit" %in% labels)
  expect_true("Q3" %in% labels)

  # without them the lines stay and the names go, and the deltas stay: they are
  # the content of the chart, not a decoration.
  without <- ggplot2::ggplot_build(
    .build(.signal(), .envelope(), .sheet4(), "CG1", .palette,
           show_labels = FALSE))
  bare_labels <- unlist(lapply(without$data, function(layer) layer$label))
  bare_labels <- bare_labels[!is.na(bare_labels)]
  expect_false("Upper Limit" %in% bare_labels)
  expect_equal(sum(grepl("δ", bare_labels)), 2L)
  bare_lines <- unlist(lapply(without$data, function(layer) layer$yintercept))
  expect_equal(length(bare_lines[!is.na(bare_lines)]), 5L)
})

test_that("a reference sample is marked as one even outside the envelope it defined", {
  signal <- .signal()
  signal$VALUE[signal$SAMPLE == "REF_1"] <- 0.99

  built <- ggplot2::ggplot_build(
    .build(signal, .envelope(), .sheet4(), "CG1", .palette))
  fills <- built$data[[2]]$fill

  # three points would be blue if the direction won; the reference colour is
  # there instead, and only two points are hyper-coloured.
  expect_true("cyan" %in% fills)
  expect_equal(sum(fills == "blue"), 1L)
  # and it carries no delta: a reference is not reported as an excursion.
  labels <- unlist(lapply(built$data, function(layer) layer$label))
  expect_equal(sum(grepl("δ", labels[!is.na(labels)])), 2L)
})

test_that("a non-finite threshold drops its line instead of breaking the chart", {
  envelope <- .envelope()
  envelope$q1 <- NA_real_
  envelope$inferior <- -Inf

  built <- ggplot2::ggplot_build(
    .build(.signal(), envelope, .sheet4(), "CG1", .palette))
  yintercepts <- unlist(lapply(built$data, function(layer) layer$yintercept))
  yintercepts <- yintercepts[!is.na(yintercepts)]

  # an M-value envelope is routinely infinite on one side; the finite lines are
  # still worth drawing.
  expect_setequal(yintercepts, c(0.5, 0.8, 0.6))
})

test_that("the thresholds are read and not recomputed, and a missing probe is NULL", {
  skip_on_cran()
  folder <- .threshold_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  # no file at all
  expect_null(SEMseeker:::.sem_probe_thresholds_get("CG1"))

  thresholds_path <- SEMseeker:::io_file_path_build(
    SEMseeker:::core_get_session_info()$result_folderData,
    "1_signal_thresholds", "parquet")
  polars::as_polars_df(data.frame(
    PROBE = c("CG1", "CG2"),
    signal_inferior_thresholds = c(0.2, 0.1),
    signal_superior_thresholds = c(0.8, 0.9),
    signal_median_values = c(0.5, 0.5),
    iqr = c(0.2, 0.3), q1 = c(0.4, 0.35), q3 = c(0.6, 0.65),
    stringsAsFactors = FALSE))$write_parquet(thresholds_path)

  envelope <- SEMseeker:::.sem_probe_thresholds_get("CG1")
  expect_equal(envelope$superior, 0.8)
  expect_equal(envelope$inferior, 0.2)
  expect_equal(envelope$q1, 0.4)
  # the case of the name is not the caller's problem.
  expect_equal(SEMseeker:::.sem_probe_thresholds_get("cg2")$superior, 0.9)
  # a probe with no thresholds is nothing to draw.
  expect_null(SEMseeker:::.sem_probe_thresholds_get("CG_ABSENT"))
})

test_that("the door writes a file and refuses to invent an envelope", {
  skip_on_cran()
  folder <- .threshold_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  signal_path <- SEMseeker:::io_pivot_file_name_parquet(
    "SIGNAL", SEMseeker:::io_signal_figure(), "PROBE", "WHOLE")
  dir.create(dirname(signal_path), recursive = TRUE, showWarnings = FALSE)
  polars::as_polars_df(data.frame(AREA = c("CG1", "CG2"),
                                  CASE_1 = c(0.95, 0.5), REF_1 = c(0.52, 0.5),
                                  stringsAsFactors = FALSE)
                       )$write_parquet(signal_path)

  sheet <- data.frame(Sample_ID = c("CASE_1", "REF_1"),
                      Sample_Group = c("Case", "Reference"),
                      stringsAsFactors = FALSE)

  # the signal is there and the thresholds are not: nothing is drawn, because a
  # threshold recomputed from the samples at hand is a different threshold.
  expect_null(SEMseeker::sem_signal_threshold_per_sample_plot("CG1", sheet))

  thresholds_path <- SEMseeker:::io_file_path_build(
    SEMseeker:::core_get_session_info()$result_folderData,
    "1_signal_thresholds", "parquet")
  polars::as_polars_df(data.frame(
    PROBE = "CG1", signal_inferior_thresholds = 0.2,
    signal_superior_thresholds = 0.8, signal_median_values = 0.5,
    iqr = 0.2, q1 = 0.4, q3 = 0.6))$write_parquet(thresholds_path)

  written <- SEMseeker::sem_signal_threshold_per_sample_plot("CG1", sheet)
  expect_true(file.exists(written))
  expect_gt(file.size(written), 1000)
  expect_match(dirname(written), "SIGNAL_THRESHOLD_PER_SAMPLE")
})
