#' Heatmap of a marker across samples and areas
#'
#' One tile per (sample, area) pair: the value of the marker as the fill, the
#' areas along the x axis and the samples along the y, ordered so that the
#' phenotype groups are contiguous. It is the chart that answers whether a
#' pattern belongs to a group of samples or to a handful of them.
#'
#' It reads through the same path as [sem_marker_value_per_area_plot()] and shows
#' the same numbers; that chart collapses the samples into a cloud per area,
#' this one keeps every sample on its own row. The first answers "which areas",
#' the second "which samples".
#'
#' @section Through 0.99.6 this was `plot_create_heatmap()`, and it was dead six times over:
#' It had a manual page, no caller, and could not have run. Measured on
#' 2026-10-05:
#'
#' * `for (g in seq_len(sample_group_comb))` where `sample_group_comb` is the
#'   matrix `combn()` returns. On a character matrix `seq_len` fails; on a
#'   numeric one it returns `1`, using the first element as the length. Either
#'   way the loop never iterated over the pairs of sample groups it was written
#'   to compare, and `g` was never read inside the body.
#' * it called `read_annotated_bed()`, which does not exist in the package.
#' * the body filtered on the whole `sample_group_comb` matrix rather than on
#'   pair `g`, so it was not a pairwise comparison at all, and the title and the
#'   file name pasted the entire matrix - which means every iteration would have
#'   written the same file.
#' * the reshape cast on `SUBAREA`, a column the line above had just set to a
#'   single constant, so it produced one column; the guard three lines later
#'   requires more than two columns and was therefore always false. **The
#'   heatmap could not be reached even with everything else fixed.**
#' * `width = 2480 / ssEnv$plot_resolution`, and that field holds the string
#'   `"print"`.
#' * the melt dropped `SAMPLEID` into the variables, poisoning the fill column
#'   with sample names.
#'
#' Its data source is gone as well: the per-marker annotated bed files it read
#' were replaced by the parquet pivots, so this is a rewrite on the current
#' storage rather than a repair.
#'
#' @section The fill scale follows the marker, not a fixed ramp:
#' A diverging ramp centred on zero is right for a marker that can be negative,
#' which the deltas are, and misleading for one that cannot, which the counts
#' are not: on non-negative data it spends half its range on values that do not
#' occur and makes a small count look like a mid-point. The scale is chosen from
#' whether the values actually cross zero.
#'
#' @section Why the samples are ordered by phenotype:
#' A heatmap is read by looking for blocks. If the rows are in the order the
#' pivot happens to store them, a pattern that belongs entirely to one group is
#' scattered across the picture and invisible. The predecessor's intent was a
#' colour strip beside the rows; ordering achieves the same reading without a
#' second legend competing with the fill.
#'
#' @param marker character. Marker name, e.g. `"MUTATIONS"`, `"DELTAS"`.
#' @param figure character. `"HYPER"`, `"HYPO"` or `"BOTH"`.
#' @param area,subarea character. The region class, e.g. `"GENE"`, `"WHOLE"`.
#' @param phenotype_column character. Sample sheet column the rows are grouped
#'   by.
#' @param sample_sheet data.frame carrying `Sample_ID` and that column.
#' @param areas_selection character vector of instances to keep, `NULL` for all.
#'   A heatmap of every gene is a stripe, so naming the areas is the normal case.
#' @param aggregation character. The sixth coordinate of the pivot, as in
#'   [sem_marker_value_per_area_plot()].
#' @param overwrite logical. Redraw when the file exists.
#'
#' @return Invisibly the path of the file written, or `NULL` when there was
#'   nothing to draw.
#'
#' @examples
#' # The chart reads a stored pivot, so it needs a result folder from a
#' # completed run. Its arguments are that pivot's coordinates plus the column
#' # the samples are grouped by:
#' sample_sheet <- data.frame(
#'   Sample_ID    = c("CASE_1", "CASE_2", "REF_1"),
#'   Sample_Group = c("Case", "Case", "Reference")
#' )
#' \dontrun{
#' sem_marker_heatmap_plot(
#'   marker = "MUTATIONS", figure = "HYPER", area = "GENE", subarea = "WHOLE",
#'   phenotype_column = "Sample_Group", sample_sheet = sample_sheet,
#'   aggregation = "SUM", areas_selection = c("BRCA1", "TP53")
#' )
#' }
#' @export
sem_marker_heatmap_plot <- function(marker, figure, area, subarea,
                                    phenotype_column, sample_sheet,
                                    areas_selection = NULL,
                                    aggregation = NULL,
                                    overwrite = FALSE) {

  ssEnv <- core_get_session_info()

  chart_folder <- io_dir_check_and_create(ssEnv$result_folderChart,
                                          "MARKER_HEATMAP")
  filename <- io_file_path_build(
    chart_folder,
    toupper(c("HEATMAP", marker, figure, area, subarea, aggregation,
              core_name_cleaning(phenotype_column))),
    ssEnv$plot_format)

  if (file.exists(filename) && !overwrite) {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " chart already exists, skipped: ", filename)
    return(invisible(filename))
  }

  pivot_data <- io_pivot_to_long_format(marker, figure, area, subarea,
                                        phenotype_column, sample_sheet,
                                        areas_selection,
                                        aggregation = aggregation)

  if (is.null(pivot_data) || nrow(pivot_data) == 0L) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " no rows to plot for ", marker, "/", figure, "/", area, "/",
              subarea, ": nothing drawn.")
    return(invisible(NULL))
  }

  plot_object <- .sem_marker_heatmap_build(pivot_data, marker, figure, area,
                                           subarea, phenotype_column,
                                           ssEnv$color_palette)

  ggplot2::ggsave(filename, plot = plot_object, width = 12, height = 9,
                  units = "in", dpi = as.numeric(ssEnv$plot_resolution_ppi))

  invisible(filename)
}

# The drawing, separated from the door, so the scale it chose and the order it
# put the rows in can be asserted on without a file.
.sem_marker_heatmap_build <- function(pivot_data, marker, figure, area, subarea,
                                      phenotype_column, palette,
                                      max_x_labels = 60L,
                                      max_y_labels = 80L) {

  # Rows grouped by phenotype, so a pattern belonging to one group reads as a
  # block instead of being scattered down the picture. Within a group the
  # samples keep the order they arrived in, which is the order of the pivot.
  sample_order <- unique(pivot_data$SAMPLE[
    order(as.character(pivot_data$phenotype))])
  pivot_data$SAMPLE <- factor(pivot_data$SAMPLE, levels = sample_order)
  pivot_data$AREA <- factor(pivot_data$AREA,
                            levels = unique(as.character(pivot_data$AREA)))

  plot_object <- ggplot2::ggplot(
      pivot_data,
      ggplot2::aes(x = .data$AREA, y = .data$SAMPLE, fill = .data$VALUE)) +
    ggplot2::geom_tile() +
    ggplot2::labs(title = paste(marker, figure, "per sample over", area,
                                subarea),
                  x = paste(area, subarea), y = phenotype_column,
                  fill = marker) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text.x = if (nlevels(pivot_data$AREA) > max_x_labels)
        ggplot2::element_blank()
      else
        ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = if (nlevels(pivot_data$SAMPLE) > max_y_labels)
        ggplot2::element_blank()
      else
        ggplot2::element_text())

  finite_values <- pivot_data$VALUE[is.finite(pivot_data$VALUE)]
  crosses_zero <- length(finite_values) > 0L &&
    any(finite_values < 0) && any(finite_values > 0)

  if (crosses_zero)
    return(plot_object +
      ggplot2::scale_fill_gradient2(low = palette[2], mid = "white",
                                    high = palette[1], midpoint = 0))

  plot_object +
    ggplot2::scale_fill_gradient(low = "white", high = palette[1])
}
