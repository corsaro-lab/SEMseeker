#' Marker values per instance of a region class, coloured by phenotype
#'
#' One point per sample per instance of the region class - one column of points
#' per gene, per island, per cytoband - with the points jittered so overlapping
#' samples stay visible and filled by the phenotype. It answers "how is this
#' marker distributed across the instances of this class, and do the phenotypes
#' sit differently within them", which is a description and not a test: no
#' p-value is computed or drawn.
#'
#' @section It is not a Manhattan plot, and it used to be called one:
#' Through 0.99.5 this was `plot_manhattan_plot_per_area()`. A Manhattan plot
#' puts genomic position on the x axis and a significance on the y; this puts a
#' categorical region-class instance on the x and a marker value on the y, which
#' is a different chart answering a different question. The two Manhattan
#' drawings the package does have are
#' `anno_manhattan_plot_marker_per_probe()` and
#' `sem_manhattan_plot_marker_per_sample()`.
#'
#' @section What it could not do before:
#' The old function could never draw anything, for reasons that compounded:
#' \itemize{
#'   \item it called `io_pivot_to_long_format()` with four arguments where six
#'     are required - `phenotype_column` and `sample_sheet` have no defaults and
#'     are used - so it stopped there, before any drawing;
#'   \item it read `pivot_data[, area]`, but the long frame names that column
#'     `AREA`, so with `area = "GENE"` the subscript was an undefined column;
#'   \item `scale_fill_gradient()` is continuous while a phenotype is usually
#'     categorical, which is `Discrete value supplied to a continuous scale`;
#'   \item `color = as.numeric(pivot_data[, area])` coerced gene symbols, so the
#'     colour was a vector of `NA`;
#'   \item its documentation promised a saved PNG and the body saved nothing, it
#'     returned the plot object - so the error surfaced only when something tried
#'     to render it, which is why calling it appeared to work;
#'   \item `family` and `adjust_method` were in the signature and never used, and
#'     `only_significant_areas` had an empty body.
#' }
#'
#' @section Two deliberate changes of behaviour:
#' The scale of the fill now follows the TYPE of the phenotype: a discrete scale
#' for a factor or a character, a gradient for a number. A phenotype in this
#' package is usually a group, so the continuous scale was wrong for the common
#' case rather than the rare one.
#'
#' And negative values are drawn where they are. The old body floored them with
#' `ifelse(VALUE < 0, 0, VALUE)`, which is dead for every count and delta marker
#' because those are non-negative by construction - but `SIGNAL` on the
#' `MVALUE` figure is negative for about half its range, and flooring drew all of
#' that at zero.
#'
#' `only_significant_areas` is gone rather than implemented. Restricting to the
#' significant instances means reading the inference file and deciding which
#' adjustment family defines significance, which is a different piece of work;
#' to look at a chosen set of instances, name them in `areas_selection`.
#'
#' @param marker character. SEM marker to plot, e.g. `"MUTATIONS"`, `"DELTARP"`,
#'   `"SIGNAL"`.
#' @param figure character. The figure of that marker: `"HYPER"` or `"HYPO"` for
#'   the epimutation markers, `"BETA"` or `"MVALUE"` for `SIGNAL`.
#' @param area,subarea character. The region class, e.g. `"GENE"` and
#'   `"TSS1500"`.
#' @param phenotype_column character. Sample sheet column the fill is taken
#'   from, e.g. `"Sample_Group"`.
#' @param sample_sheet data.frame. The sample sheet the phenotype is read from.
#' @param areas_selection character vector. Instances to keep, `NULL` for all.
#' @param overwrite logical. Redraw when the file exists. `FALSE` by default,
#'   like the other charts.
#' @return Invisibly the path of the file written, or `NULL` when there was
#'   nothing to draw.
#'
#' @examples
#' \donttest{
#' # Needs a result folder with the pivots of the requested class.
#' invisible(NULL)
#' }
#' @export
sem_marker_value_per_area_plot <- function(marker, figure, area, subarea,
                                           phenotype_column, sample_sheet,
                                           areas_selection = NULL,
                                           overwrite = FALSE) {

  ssEnv <- core_get_session_info()

  chart_folder <- io_dir_check_and_create(ssEnv$result_folderChart,
                                          "MARKER_VALUE_PER_AREA")
  # The coordinates are in the name, so two markers or two figures of the same
  # class do not overwrite each other in the same folder.
  filename <- io_file_path_build(
    chart_folder,
    toupper(c("VALUE_PER_AREA", marker, figure, area, subarea,
              core_name_cleaning(phenotype_column))),
    ssEnv$plot_format)

  if (file.exists(filename) && !overwrite) {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " chart already exists, skipped: ", filename)
    return(invisible(filename))
  }

  pivot_data <- io_pivot_to_long_format(marker, figure, area, subarea,
                                        phenotype_column, sample_sheet,
                                        areas_selection)

  if (is.null(pivot_data) || nrow(pivot_data) == 0L) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " no rows to plot for ", marker, "/", figure, "/", area, "/",
              subarea, ": nothing drawn.")
    return(invisible(NULL))
  }

  plot_object <- .sem_area_value_plot_build(pivot_data, marker, figure, area,
                                            subarea, phenotype_column,
                                            ssEnv$color_palette)

  ggplot2::ggsave(filename, plot = plot_object, width = 9, height = 9,
                  units = "in", dpi = as.numeric(ssEnv$plot_resolution_ppi))

  invisible(filename)
}

# The drawing, separated from the door that reads and saves. Two reasons, and the
# second is the one that matters: a function that returns the object can be
# asserted on - the y values it maps, which scale it chose - where a saved PNG
# can only be asserted to exist. And it is the shape the plotting surface is
# moving to, the layout once and a thin door per domain in front of it.
#
# The palette arrives as an argument rather than being read from the session, so
# the drawing needs no session at all.
.sem_area_value_plot_build <- function(pivot_data, marker, figure, area, subarea,
                                       phenotype_column, palette) {

  plot_object <- ggplot2::ggplot(
      pivot_data,
      ggplot2::aes(x = factor(.data$AREA), y = .data$VALUE,
                   fill = .data$phenotype)) +
    ggplot2::geom_point(shape = 21, alpha = 0.7,
                        position = ggplot2::position_jitter(width = 0.2,
                                                            height = 0)) +
    ggplot2::labs(title = paste(marker, figure, "per", area, subarea),
                  x = paste(area, subarea), y = marker, fill = phenotype_column) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

  # The fill scale follows the TYPE of the phenotype rather than assuming one. A
  # phenotype here is usually a group, so a continuous scale was wrong for the
  # common case and not for the rare one.
  if (is.numeric(pivot_data$phenotype))
    return(plot_object +
      ggplot2::scale_fill_gradient(low = "white", high = palette[1]))

  plot_object +
    ggplot2::scale_fill_manual(values = grDevices::colorRampPalette(palette)(
      length(unique(pivot_data$phenotype))))
}
