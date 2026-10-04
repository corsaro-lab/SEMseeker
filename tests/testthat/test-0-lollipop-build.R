# .plot_lollipop_build(), the layout shared by the marker and the threshold
# charts. Every assertion goes through ggplot2::ggplot_build(), which renders
# the layers: an object that merely exists proves nothing, and the charts this
# replaces all built fine and failed when something tried to draw them.

.lolli <- function(...) SEMseeker:::.plot_lollipop_build(...)

.classes <- c(Hyper = "blue", Hypo = "orange",
              "Non-outlier" = "grey", Reference = "cyan")

.data3 <- function() data.frame(
  X = c("S3", "S1", "S2"),
  VALUE = c(5, -2, 0),
  CLASS = c("Hyper", "Hypo", "Non-outlier"),
  BASELINE = 0,
  SEGMENT_COLOUR = c("blue", "orange", "white"),
  stringsAsFactors = FALSE)

test_that("the layout renders, and renders every row", {
  built <- ggplot2::ggplot_build(
    .lolli(.data3(), "Sample", "DELTAS", .classes))

  # two layers at minimum: the segments and the points.
  expect_gte(length(built$data), 2L)
  expect_equal(nrow(built$data[[1]]), 3L)
  expect_equal(nrow(built$data[[2]]), 3L)
  # the y values reach the renderer unchanged, negatives included.
  expect_setequal(built$data[[2]]$y, c(5, -2, 0))
})

test_that("the x axis carries the names, in the order the caller gave", {
  # a character X keeps the row order rather than being sorted.
  built <- ggplot2::ggplot_build(
    .lolli(.data3(), "Sample", "DELTAS", .classes))
  expect_equal(built$layout$panel_params[[1]]$x$get_labels(),
               c("S3", "S1", "S2"))

  # and a factor's level order wins, which is how a declared order arrives.
  d <- .data3()
  d$X <- factor(d$X, levels = c("S1", "S2", "S3"))
  built_ordered <- ggplot2::ggplot_build(
    .lolli(d, "Sample", "DELTAS", .classes))
  expect_equal(built_ordered$layout$panel_params[[1]]$x$get_labels(),
               c("S1", "S2", "S3"))
})

test_that("segment colour is per row and point fill is per class", {
  built <- ggplot2::ggplot_build(
    .lolli(.data3(), "Sample", "DELTAS", .classes))

  # the segment layer takes its colour from the data, identity-scaled.
  expect_setequal(built$data[[1]]$colour, c("blue", "orange", "white"))
  # the point layer takes its fill from the class table, not from the data.
  fills <- built$data[[2]]$fill
  expect_true("blue" %in% fills)
  expect_true("grey" %in% fills)
  # the two channels are independent: a row whose segment is white still gets
  # the colour of its class on the point.
  white_segment <- which(built$data[[1]]$colour == "white")
  expect_equal(fills[white_segment], "grey")
})

test_that("the segment runs from the baseline to the value", {
  d <- .data3()
  d$BASELINE <- c(1, -1, 0)
  built <- ggplot2::ggplot_build(.lolli(d, "Sample", "DELTAS", .classes))
  expect_setequal(built$data[[1]]$y, c(1, -1, 0))
  expect_setequal(built$data[[1]]$yend, c(5, -2, 0))
})

test_that("reference lines render at their values and carry their labels", {
  lines <- data.frame(
    VALUE = c(0.2, 0.8, 0.5),
    LABEL = c("Lower Limit", "Upper Limit", ""),
    COLOUR = c("red", "red", "black"),
    LINETYPE = "dashed",
    stringsAsFactors = FALSE)

  built <- ggplot2::ggplot_build(
    .lolli(.data3(), "Sample", "DELTAS", .classes,
           reference_lines = lines))

  yintercepts <- unlist(lapply(built$data, function(layer) layer$yintercept))
  expect_setequal(yintercepts[!is.na(yintercepts)], c(0.2, 0.8, 0.5))

  # the unlabelled line is drawn and not labelled: two labels for three lines.
  labels <- unlist(lapply(built$data, function(layer) layer$label))
  expect_setequal(labels[!is.na(labels)], c("Lower Limit", "Upper Limit"))
})

test_that("no reference lines means no extra layers", {
  bare <- ggplot2::ggplot_build(.lolli(.data3(), "S", "V", .classes))
  lined <- ggplot2::ggplot_build(
    .lolli(.data3(), "S", "V", .classes,
           reference_lines = data.frame(VALUE = 1, LABEL = "x",
                                        COLOUR = "red", LINETYPE = "dashed",
                                        stringsAsFactors = FALSE)))
  expect_lt(length(bare$data), length(lined$data))
})

test_that("annotations render one label per point", {
  notes <- data.frame(X = c("S3", "S1"), VALUE = c(5, -2),
                      LABEL = c("d 1.2e-03", "d 4.0e-04"),
                      COLOUR = c("blue", "orange"),
                      stringsAsFactors = FALSE)
  built <- ggplot2::ggplot_build(
    .lolli(.data3(), "Sample", "DELTAS", .classes, annotations = notes))

  labels <- unlist(lapply(built$data, function(layer) layer$label))
  labels <- labels[!is.na(labels)]
  # the defect this replaces: ifelse(labels, "") coerced the vector to NA and
  # geom_text drew nothing, with no error anywhere.
  expect_equal(length(labels), 2L)
  expect_true(all(nzchar(labels)))
})

test_that("a missing column is a refusal, named", {
  d <- .data3()
  d$BASELINE <- NULL
  expect_error(.lolli(d, "S", "V", .classes), "BASELINE")

  d2 <- .data3()
  d2$CLASS <- NULL
  d2$SEGMENT_COLOUR <- NULL
  expect_error(.lolli(d2, "S", "V", .classes), "CLASS")
})

test_that("a class with no colour is refused, because ggplot2 would not refuse it", {
  d <- .data3()
  d$CLASS[1] <- "Unexpected"
  # Measured: scale_fill_manual() draws an unknown class in grey50, which is
  # exactly what a legitimate Non-outlier looks like here. So the drawing checks.
  expect_error(.lolli(d, "S", "V", .classes), "Unexpected")

  # and the message says which classes it does know.
  expect_error(.lolli(d, "S", "V", .classes), "Non-outlier")
})

test_that("one row draws, and the axis titles are the ones given", {
  one <- data.frame(X = "ONLY", VALUE = 3, CLASS = "Hyper",
                    BASELINE = 0, SEGMENT_COLOUR = "blue",
                    stringsAsFactors = FALSE)
  plot_object <- .lolli(one, "Probe", "MUTATIONS", .classes, title = "T")
  built <- ggplot2::ggplot_build(plot_object)
  expect_equal(nrow(built$data[[2]]), 1L)
  expect_equal(plot_object$labels$x, "Probe")
  expect_equal(plot_object$labels$y, "MUTATIONS")
  expect_equal(plot_object$labels$title, "T")
})

test_that("the x axis text is dropped once the categories outnumber the limit", {
  many <- data.frame(X = paste0("cg", sprintf("%05d", 1:80)),
                     VALUE = seq_len(80), CLASS = "Hyper",
                     BASELINE = 0, SEGMENT_COLOUR = "blue",
                     stringsAsFactors = FALSE)

  # 80 names rotated 90 degrees are an unreadable band, not a label.
  wide <- .lolli(many, "Probe", "MUTATIONS", .classes)
  expect_s3_class(wide$theme$axis.text.x, "element_blank")

  # but the categories themselves are all still drawn.
  expect_equal(nrow(ggplot2::ggplot_build(wide)$data[[2]]), 80L)

  # under the limit the names stay.
  narrow <- .lolli(many[1:10, ], "Probe", "MUTATIONS", .classes)
  expect_false(inherits(narrow$theme$axis.text.x, "element_blank"))

  # and the limit is the caller's to move.
  raised <- .lolli(many, "Probe", "MUTATIONS", .classes, max_x_labels = 100L)
  expect_false(inherits(raised$theme$axis.text.x, "element_blank"))
})
