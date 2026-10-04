# sem_marker_value_per_area_plot(), which through 0.99.5 was
# plot_manhattan_plot_per_area() and could never draw anything.
#
# The drawing is asserted on the OBJECT, not on the PNG: a saved file can only be
# checked to exist, where the object says which values it mapped and which scale
# it chose, and those are the two things that were wrong.

.long <- function(phenotype = c("case", "ctrl")) {
  data.frame(
    AREA      = rep(c("BRCA1", "TP53", "EGFR"), each = 2),
    VALUE     = c(1, 2, 3, 4, 5, 6),
    phenotype = rep(phenotype, 3),
    stringsAsFactors = FALSE)
}
.build <- function(...) SEMseeker:::.sem_area_value_plot_build(...)
.pal   <- c("#1f78b4", "#e31a1c", "#33a02c")

# ---- the two defects that stopped it drawing ---------------------------------

test_that("the region-class instance comes from AREA, which is the column the data has", {
  # It used to read pivot_data[, area] with area = "GENE", and the long frame has
  # AREA / VALUE / phenotype: an undefined column, every time.
  d <- .long()
  expect_false("GENE" %in% colnames(d))
  g <- .build(d, "MUTATIONS", "HYPER", "GENE", "WHOLE", "Sample_Group", .pal)
  built <- ggplot2::ggplot_build(g)

  # The x axis is the instances of the region class, taken from AREA.
  expect_setequal(built$layout$panel_params[[1]]$x$get_labels(),
                  c("BRCA1", "EGFR", "TP53"))
  # and every row of the long frame is a point: nothing dropped on the way
  expect_equal(nrow(built$data[[1]]), nrow(d))
})

test_that("a categorical phenotype gets a discrete scale, which is the common case", {
  # scale_fill_gradient() is continuous, so a phenotype of groups was
  # "Discrete value supplied to a continuous scale" and nothing was drawn.
  g <- .build(.long(), "MUTATIONS", "HYPER", "GENE", "WHOLE", "Sample_Group", .pal)
  expect_silent(invisible(ggplot2::ggplot_build(g)))
  scales <- vapply(g$scales$scales, function(s) class(s)[1], "")
  expect_true(any(grepl("ScaleDiscrete", scales)))
})

test_that("a numeric phenotype still gets a gradient", {
  g <- .build(.long(phenotype = c(55, 71)), "DELTARP", "HYPO", "ISLAND", "WHOLE",
              "Age", .pal)
  expect_silent(invisible(ggplot2::ggplot_build(g)))
  scales <- vapply(g$scales$scales, function(s) class(s)[1], "")
  expect_true(any(grepl("ScaleContinuous", scales)))
})

# ---- the flooring that hid half of SIGNAL on the M scale --------------------

test_that("negative values are drawn where they are, not at zero", {
  # The old body mapped ifelse(VALUE < 0, 0, VALUE). For counts and deltas that
  # is dead, because they are non-negative by construction; for SIGNAL on the
  # MVALUE figure it drew about half the range at zero.
  d <- .long()
  d$VALUE <- c(-3.2, -1.1, 0, 0.5, 2.0, 4.1)
  g <- .build(d, "SIGNAL", "MVALUE", "PROBE", "WHOLE", "Sample_Group", .pal)
  y <- ggplot2::ggplot_build(g)$data[[1]]$y

  expect_equal(sort(y), sort(d$VALUE))
  expect_lt(min(y), 0)
  expect_equal(sum(y == 0), 1L)   # the one true zero, not three
})

# ---- the door ---------------------------------------------------------------

test_that("the file name carries the coordinates, so two markers cannot collide", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  orig <- SEMseeker:::io_pivot_to_long_format
  unlockBinding("io_pivot_to_long_format", asNamespace("SEMseeker"))
  assign("io_pivot_to_long_format", function(...) .long(), envir = asNamespace("SEMseeker"))
  on.exit({
    assign("io_pivot_to_long_format", orig, envir = asNamespace("SEMseeker"))
    lockBinding("io_pivot_to_long_format", asNamespace("SEMseeker"))
  }, add = TRUE)

  a <- SEMseeker:::sem_marker_value_per_area_plot("MUTATIONS", "HYPER", "GENE",
                                                  "WHOLE", "Sample_Group", NULL)
  b <- SEMseeker:::sem_marker_value_per_area_plot("DELTARP", "HYPER", "GENE",
                                                  "WHOLE", "Sample_Group", NULL)
  expect_false(identical(a, b))
  expect_true(file.exists(a)); expect_true(file.exists(b))
  expect_match(basename(a), "MUTATIONS")
  expect_match(basename(b), "DELTARP")
})

test_that("nothing to draw is a warning and a NULL, not an empty chart", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  orig <- SEMseeker:::io_pivot_to_long_format
  unlockBinding("io_pivot_to_long_format", asNamespace("SEMseeker"))
  assign("io_pivot_to_long_format", function(...) .long()[0, ], envir = asNamespace("SEMseeker"))
  on.exit({
    assign("io_pivot_to_long_format", orig, envir = asNamespace("SEMseeker"))
    lockBinding("io_pivot_to_long_format", asNamespace("SEMseeker"))
  }, add = TRUE)

  out <- SEMseeker:::sem_marker_value_per_area_plot("MUTATIONS", "HYPER", "GENE",
                                                     "WHOLE", "Sample_Group", NULL)
  expect_null(out)
})

# ---------------------------------------------------------------------------
# End to end, through the REAL reader. The tests above stub
# io_pivot_to_long_format(), which is what let its defects through: the chart
# was repaired while the only path to a real pivot stayed broken. This one
# writes a parquet pivot and asserts on the file that comes out.
# ---------------------------------------------------------------------------

test_that("the chart is drawn from a real pivot on disk", {
  skip_on_cran()

  folder <- sem_test_folder()
  dir.create(file.path(folder, "Data"), recursive = TRUE, showWarnings = FALSE)
  SEMseeker:::core_init_env(folder, parallel_strategy = "sequential",
                            iqrTimes = 3, verbosity = 1)
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  pivot_path <- SEMseeker:::io_pivot_file_name_parquet("MUTATIONS", "HYPER",
                                                      "GENE", "WHOLE",
                                                      aggregation = "SUM")
  dir.create(dirname(pivot_path), recursive = TRUE, showWarnings = FALSE)
  polars::as_polars_df(data.frame(AREA = c("BRCA1", "TP53"),
                                  CASE_1 = c(3, -2), REF_1 = c(0, 1),
                                  stringsAsFactors = FALSE)
                       )$write_parquet(pivot_path)

  sheet <- data.frame(Sample_ID = c("CASE_1", "REF_1"),
                      Sample_Group = c("Case", "Reference"),
                      stringsAsFactors = FALSE)

  written <- SEMseeker:::sem_marker_value_per_area_plot(
    "MUTATIONS", "HYPER", "GENE", "WHOLE", "Sample_Group", sheet,
    aggregation = "SUM")

  expect_true(!is.null(written))
  expect_true(file.exists(written))
  expect_gt(file.size(written), 1000)
  # the sixth coordinate is in the name, so SUM and MEAN of the same class do
  # not land on the same file.
  expect_match(basename(written), "SUM")

  # asking again without overwrite does not redraw, and says which file it is.
  again <- SEMseeker:::sem_marker_value_per_area_plot(
    "MUTATIONS", "HYPER", "GENE", "WHOLE", "Sample_Group", sheet,
    aggregation = "SUM")
  expect_equal(again, written)
})
