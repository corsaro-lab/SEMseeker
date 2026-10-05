#' The lollipop layout, once, for every domain that needs it
#'
#' A point per observation over a categorical x axis, each joined to a baseline
#' by a segment, with the point filled by the class the observation belongs to.
#' Optional horizontal reference lines and per-point annotations turn the same
#' layout into the threshold chart without a second copy of it.
#'
#' It lives under the transversal `plot_` prefix and is internal. A layout
#' belongs to no domain: the same construction serves a marker across samples, a
#' marker across probes and a signal against its reference envelope. In front of
#' it stand thin doors, one per question, each with its own name and manual page,
#' because it is the name that says which question a chart answers.
#'
#' @section Why the two colour channels do not collide:
#' Segments carry their colour as data and are drawn through
#' `scale_colour_identity()`; points carry their class and are drawn through
#' `scale_fill_manual()` on a `shape = 21` point. Two different aesthetics, so a
#' per-observation segment colour and a per-class fill can coexist. The previous
#' code passed a vector to `colour =` outside `aes()`, which works only while the
#' vector stays in exactly the row order of the data and breaks silently when it
#' does not.
#'
#' @section The x axis carries names, not indices:
#' `x` is used as a factor, so the axis shows the categories and their order is
#' the order of the levels, which the caller sets. The previous code wrote
#' `as.numeric(factor(samples))`, which produces a continuous axis labelled 1..n
#' while keeping the 90-degree rotated text that only makes sense for names: the
#' charts had no way to tell which sample or which probe a point belonged to.
#'
#' @param plot_data data.frame with `X` (the category), `VALUE` (y), `CLASS` (a
#'   name of `class_colours`), `BASELINE` (where the segment starts) and
#'   `SEGMENT_COLOUR` (one colour per row). `X` may be a factor, and then its
#'   level order is the order of the axis.
#' @param x_label,y_label character. Axis titles.
#' @param class_colours named character vector, one colour per class. Every class
#'   present in `CLASS` must have one. Measured: `scale_fill_manual()` does not
#'   refuse a class it has no colour for, it draws it in `grey50`, which is
#'   indistinguishable from a legitimate grey class - a category silently
#'   recoloured as another. So the check is here instead.
#' @param reference_lines data.frame with `VALUE`, `LABEL`, `COLOUR` and
#'   `LINETYPE`, or `NULL`. One horizontal line each, labelled at the left edge.
#'   An empty `LABEL` draws the line without the label.
#' @param annotations data.frame with `X`, `VALUE`, `LABEL` and `COLOUR`, or
#'   `NULL`. One rotated label per point, overlap-checked.
#' @param segment_alpha numeric. Opacity of the segments.
#' @param label_font_size numeric. Size of the annotation labels.
#' @param max_x_labels integer. Above this many distinct categories the x axis
#'   text is dropped. A Manhattan-style chart legitimately carries thousands of
#'   probes, and thousands of names rotated 90 degrees are an unreadable black
#'   band that also takes most of the rendering time. The axis title still says
#'   what the dimension is.
#' @param title character. Plot title, empty by default.
#'
#' @return A `ggplot` object. Nothing is saved: the door that called this decides
#'   where the file goes, and a function that returns the object can be asserted
#'   on, where a saved PNG can only be asserted to exist.
#'
#' @keywords internal
#' @noRd
.plot_lollipop_build <- function(plot_data, x_label, y_label, class_colours,
                                 reference_lines = NULL, annotations = NULL,
                                 segment_alpha = 0.5, label_font_size = 3,
                                 max_x_labels = 60L, title = "") {

  required <- c("X", "VALUE", "CLASS", "BASELINE", "SEGMENT_COLOUR")
  missing_columns <- setdiff(required, colnames(plot_data))
  if (length(missing_columns) > 0L)
    stop("plot_data is missing: ", paste(missing_columns, collapse = ", "))

  # ggplot2 would draw an unknown class in grey50 rather than refuse it, and
  # grey is already what a non-outlier looks like. A chart that recolours a
  # category as another category is worse than no chart.
  uncoloured <- setdiff(unique(as.character(plot_data$CLASS)),
                        names(class_colours))
  if (length(uncoloured) > 0L)
    stop("class_colours has no colour for: ",
         paste(uncoloured, collapse = ", "),
         ". Known classes: ", paste(names(class_colours), collapse = ", "))

  # A factor keeps the caller's order; a character vector is ordered by the
  # order it appears in, which for a pivot is the order of the rows on disk.
  if (!is.factor(plot_data$X))
    plot_data$X <- factor(plot_data$X, levels = unique(as.character(plot_data$X)))

  plot_object <- ggplot2::ggplot(plot_data) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$X, xend = .data$X,
                   y = .data$BASELINE, yend = .data$VALUE,
                   colour = .data$SEGMENT_COLOUR),
      alpha = segment_alpha) +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$X, y = .data$VALUE, fill = .data$CLASS),
      shape = 21, alpha = 0.7, stroke = 0.1) +
    ggplot2::scale_fill_manual(values = class_colours) +
    ggplot2::labs(x = x_label, y = y_label, title = title, fill = NULL) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = if (nlevels(plot_data$X) > max_x_labels)
        ggplot2::element_blank()
      else
        ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.title = ggplot2::element_text(hjust = 0.5))

  if (!is.null(reference_lines) && nrow(reference_lines) > 0L) {
    plot_object <- plot_object +
      ggplot2::geom_hline(
        data = reference_lines,
        ggplot2::aes(yintercept = .data$VALUE, colour = .data$COLOUR,
                     linetype = .data$LINETYPE),
        linewidth = 0.2) +
      ggplot2::scale_linetype_identity()

    labelled <- reference_lines[nzchar(reference_lines$LABEL), , drop = FALSE]
    if (nrow(labelled) > 0L)
      plot_object <- plot_object +
        ggplot2::geom_text(
          data = labelled,
          ggplot2::aes(x = 1, y = .data$VALUE, label = .data$LABEL,
                       colour = .data$COLOUR),
          hjust = -0.1, vjust = -0.4, size = label_font_size)
  }

  if (!is.null(annotations) && nrow(annotations) > 0L)
    plot_object <- plot_object +
      ggplot2::geom_text(
        data = annotations,
        ggplot2::aes(x = .data$X, y = .data$VALUE, label = .data$LABEL,
                     colour = .data$COLOUR),
        angle = 90, check_overlap = TRUE, size = label_font_size,
        hjust = 1, vjust = -0.45)

  plot_object
}
