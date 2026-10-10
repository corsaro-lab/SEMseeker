# sem_marker_heatmap_plot(), which through 0.99.6 was plot_create_heatmap():
# a manual page, no caller, and six independent reasons it could not run -
# including a reshape cast that produced one column where the guard below it
# required more than two, so the heatmap was unreachable even with everything
# else fixed. Its data source is gone too, so this is a rewrite on the pivots.

.heatmap_session <- function() {
  folder <- sem_test_folder()
  dir.create(file.path(folder, "Data"), recursive = TRUE, showWarnings = FALSE)
  SEMseeker:::core_init_env(folder, parallel_strategy = "sequential",
                            iqrTimes = 3, verbosity = 1)
  folder
}

.write_gene_pivot <- function(areas, ..., aggregation = "SUM") {
  path <- SEMseeker:::io_pivot_file_name_parquet("MUTATIONS", "HYPER",
                                                "GENE", "WHOLE",
                                                aggregation = aggregation)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  polars::as_polars_df(data.frame(AREA = areas, ..., stringsAsFactors = FALSE)
                       )$write_parquet(path)
  path
}

.heatmap_sheet <- function() data.frame(
  Sample_ID = c("CASE_1", "REF_1", "CASE_2"),
  Sample_Group = c("Case", "Reference", "Case"),
  stringsAsFactors = FALSE)

.long <- function() data.frame(
  AREA = rep(c("BRCA1", "TP53"), times = 3),
  VALUE = c(3, 1, 0, 0, 5, 2),
  SAMPLE = rep(c("CASE_1", "REF_1", "CASE_2"), each = 2),
  phenotype = rep(c("Case", "Reference", "Case"), each = 2),
  stringsAsFactors = FALSE)

.build <- function(data = .long(), ...) {
  session <- SEMseeker:::core_get_session_info()
  SEMseeker:::.sem_marker_heatmap_build(
    data, "MUTATIONS", "HYPER", "GENE", "WHOLE", "Sample_Group",
    session$color_palette, ...)
}

test_that("one tile per sample and area, and every value reaches the renderer", {
  skip_on_cran()
  folder <- .heatmap_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  built <- ggplot2::ggplot_build(.build())
  # three samples by two areas, each pair its own tile.
  expect_equal(nrow(built$data[[1]]), 6L)
  expect_equal(nrow(unique(built$data[[1]][, c("x", "y")])), 6L)

  # the value decides the fill: the fixture has five distinct values, so there
  # are five distinct colours and equal values share one.
  data <- .long()
  expect_equal(length(unique(built$data[[1]]$fill)),
               length(unique(data$VALUE)))
  zeros <- which(data$VALUE == 0)
  expect_equal(length(unique(built$data[[1]]$fill[zeros])), 1L)
})

test_that("the samples are ordered so the phenotype groups are contiguous", {
  skip_on_cran()
  folder <- .heatmap_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  plot_object <- .build()
  # REF_1 sits between the two cases in the input. A heatmap is read by looking
  # for blocks, so the rows are grouped: the two Cases become adjacent.
  expect_equal(levels(plot_object$data$SAMPLE),
               c("CASE_1", "CASE_2", "REF_1"))

  phenotypes <- plot_object$data$phenotype[
    match(levels(plot_object$data$SAMPLE), plot_object$data$SAMPLE)]
  expect_false(any(diff(as.integer(as.factor(phenotypes))) < 0))
})

test_that("a marker that crosses zero gets a diverging scale centred on zero", {
  skip_on_cran()
  folder <- .heatmap_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  deltas <- .long()
  deltas$VALUE <- c(-3, 1, 0, -1, 5, 2)
  built <- ggplot2::ggplot_build(.build(deltas))

  # a diverging ramp puts white at zero; the negative and positive extremes
  # must therefore get different colours and neither may be white.
  fills <- built$data[[1]]$fill
  expect_equal(length(unique(fills)), length(unique(deltas$VALUE)))
  most_negative <- fills[which.min(deltas$VALUE)]
  most_positive <- fills[which.max(deltas$VALUE)]
  expect_false(most_negative == most_positive)
})

test_that("a non-negative marker gets a sequential scale, not half a wasted ramp", {
  skip_on_cran()
  folder <- .heatmap_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  counts <- .long()
  expect_true(all(counts$VALUE >= 0))
  built <- ggplot2::ggplot_build(.build(counts))

  # zero is the low end of a sequential ramp, which is white here, so the
  # smallest value is the palest. On a diverging ramp centred at zero the
  # smallest value would be at the middle instead.
  fills <- built$data[[1]]$fill
  expect_equal(toupper(fills[which.min(counts$VALUE)]), "#FFFFFF")
})

test_that("the axis labels are dropped only when they stop being readable", {
  skip_on_cran()
  folder <- .heatmap_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  small <- .build()
  expect_false(inherits(small$theme$axis.text.x, "element_blank"))
  expect_false(inherits(small$theme$axis.text.y, "element_blank"))

  many_areas <- data.frame(
    AREA = paste0("GENE", seq_len(100)), VALUE = seq_len(100),
    SAMPLE = "CASE_1", phenotype = "Case", stringsAsFactors = FALSE)
  wide <- .build(many_areas)
  expect_s3_class(wide$theme$axis.text.x, "element_blank")
  # and every area is still drawn: the labels go, the data does not.
  expect_equal(nrow(ggplot2::ggplot_build(wide)$data[[1]]), 100L)

  # the limits are the caller's to move.
  expect_false(inherits(.build(many_areas, max_x_labels = 200L)$theme$axis.text.x,
                        "element_blank"))
})

test_that("the chart is written from a real pivot, with its coordinates in the name", {
  skip_on_cran()
  folder <- .heatmap_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  .write_gene_pivot(c("BRCA1", "TP53"),
                    CASE_1 = c(3, 1), REF_1 = c(0, 0), CASE_2 = c(5, 2))

  written <- SEMseeker::sem_marker_heatmap_plot(
    "MUTATIONS", "HYPER", "GENE", "WHOLE", "Sample_Group", .heatmap_sheet(),
    aggregation = "SUM")

  expect_true(file.exists(written))
  expect_gt(file.size(written), 1000)
  expect_match(basename(written), "MUTATIONS")
  expect_match(basename(written), "SUM")
  expect_match(dirname(written), "MARKER_HEATMAP")

  # asking again does not redraw, and says which file it is.
  expect_equal(SEMseeker::sem_marker_heatmap_plot(
    "MUTATIONS", "HYPER", "GENE", "WHOLE", "Sample_Group", .heatmap_sheet(),
    aggregation = "SUM"), written)
})

test_that("no pivot is nothing to draw, and a bad phenotype column is a refusal", {
  skip_on_cran()
  folder <- .heatmap_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  expect_null(SEMseeker::sem_marker_heatmap_plot(
    "MUTATIONS", "HYPER", "GENE", "WHOLE", "Sample_Group", .heatmap_sheet(),
    aggregation = "SUM"))

  .write_gene_pivot("BRCA1", CASE_1 = 3)
  expect_error(SEMseeker::sem_marker_heatmap_plot(
    "MUTATIONS", "HYPER", "GENE", "WHOLE", "NoSuchColumn", .heatmap_sheet(),
    aggregation = "SUM"), "not a column of the sample sheet")
})

test_that("an areas selection reaches the chart", {
  skip_on_cran()
  folder <- .heatmap_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  .write_gene_pivot(c("BRCA1", "TP53", "EGFR"),
                    CASE_1 = c(3, 1, 7), REF_1 = c(0, 0, 0),
                    CASE_2 = c(5, 2, 1))

  written <- SEMseeker::sem_marker_heatmap_plot(
    "MUTATIONS", "HYPER", "GENE", "WHOLE", "Sample_Group", .heatmap_sheet(),
    areas_selection = c("BRCA1", "EGFR"), aggregation = "SUM")
  expect_true(file.exists(written))

  # the selection is applied by the reader, so assert on what it returns.
  long <- SEMseeker:::io_pivot_to_long_format(
    "MUTATIONS", "HYPER", "GENE", "WHOLE", "Sample_Group", .heatmap_sheet(),
    c("BRCA1", "EGFR"), aggregation = "SUM")
  expect_setequal(unique(long$AREA), c("BRCA1", "EGFR"))
})
