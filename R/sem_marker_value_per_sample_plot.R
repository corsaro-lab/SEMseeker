#' Marker value of one probe across samples
#'
#' One point per sample for a single probe: the value of the marker on the y
#' axis, the sample on the x, and the colour saying whether that sample is a
#' hyper-methylation event at this probe, a hypo-methylation event, a reference
#' sample or none of those.
#'
#' It is the transpose of [sem_marker_value_per_probe_plot()], over the same
#' layout, and the pair is what the two halves of the old Manhattan family were
#' trying to be.
#'
#' @section Through 0.99.6 this was `anno_manhattan_plot_marker_per_probe()`:
#' That function had a manual page and could not produce a chart. Measured on
#' 2026-10-04 it called `core_init_env()` from inside a chart and then
#' `io_get_pivot_both()`, where the cache path was a variable that function never
#' assigned, so the normal branch raised `object not found` after completing all
#' of its work. It also scanned every marker to compute five probe statistics,
#' wrote them to a CSV, and drew ten charts per call, which is why asking it for
#' one chart was not possible. The statistics are now
#' [sem_probe_select_by_statistic()] and the chart draws the probe it is given:
#' two calls instead of one, and a chart whose subject is in its own arguments.
#'
#' @section The name says the axis:
#' `per_sample` means one point per sample, which is the x axis, and the probe is
#' what the chart holds fixed. That is the opposite of what the old name's
#' `per_probe` meant, and it is the convention the rest of this family uses: the
#' axis is in the name, the fixed subject is in the signature.
#'
#' @section Reference samples are a class, not a filter:
#' A sample whose `Sample_Group` is `Reference` is drawn in its own colour rather
#' than removed. The thresholds this probe is judged against were computed from
#' those samples, so seeing where they sit is how a reader judges whether an
#' event is a real excursion or an artefact of a narrow reference envelope. The
#' class overrides the direction: a reference sample is marked as a reference
#' even where its own signal sits outside the envelope it helped define.
#'
#' @param marker character. Marker name, e.g. `"MUTATIONS"`, `"DELTAS"`.
#' @param figure character. The figure whose values go on the y axis.
#' @param probe character. The probe held fixed, e.g. `"cg11680158"`.
#' @param sample_sheet data.frame carrying `Sample_ID` and `Sample_Group`.
#' @param samples_selection character vector of sample names, or `NULL` for all
#'   the samples the pivot carries.
#' @param overwrite logical. Redraw when the file exists.
#'
#' @return Invisibly the path of the file written, or `NULL` when there was
#'   nothing to draw.
#'
#' @examples
#' \donttest{
#' # Needs a result folder holding the probe pivots of the marker.
#' invisible(NULL)
#' }
#' @export
sem_marker_value_per_sample_plot <- function(marker, figure, probe,
                                            sample_sheet,
                                            samples_selection = NULL,
                                            overwrite = FALSE) {

  ssEnv <- core_get_session_info()

  chart_folder <- io_dir_check_and_create(ssEnv$result_folderChart,
                                          "MARKER_VALUE_PER_SAMPLE")
  filename <- io_file_path_build(
    chart_folder,
    toupper(c("VALUE_PER_SAMPLE", marker, figure, core_name_cleaning(probe))),
    ssEnv$plot_format)

  if (file.exists(filename) && !overwrite) {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " chart already exists, skipped: ", filename)
    return(invisible(filename))
  }

  values <- .sem_probe_row_get(marker, figure, probe, samples_selection)
  if (is.null(values)) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " no ", marker, "/", figure, " values at probe ", probe,
              ": nothing drawn.")
    return(invisible(NULL))
  }

  hyper <- .sem_probe_row_get(marker, "HYPER", probe, values$SAMPLE)
  hypo  <- .sem_probe_row_get(marker, "HYPO",  probe, values$SAMPLE)

  class_of_sample <- rep("Non-outlier", nrow(values))
  if (!is.null(hyper))
    class_of_sample[values$SAMPLE %in% hyper$SAMPLE[hyper$VALUE != 0]] <- "Hyper"
  if (!is.null(hypo))
    class_of_sample[values$SAMPLE %in% hypo$SAMPLE[hypo$VALUE != 0]] <- "Hypo"

  references <- as.character(
    sample_sheet$Sample_ID[sample_sheet$Sample_Group == "Reference"])
  class_of_sample[values$SAMPLE %in% references] <- "Reference"

  palette <- .sem_outlier_palette(ssEnv)
  drawn_value <- ifelse(class_of_sample %in% c("Hyper", "Hypo"),
                        values$VALUE, 0)

  plot_object <- .plot_lollipop_build(
    data.frame(X = values$SAMPLE, VALUE = drawn_value, CLASS = class_of_sample,
               BASELINE = 0,
               SEGMENT_COLOUR = ifelse(class_of_sample %in% c("Hyper", "Hypo"),
                                       palette[class_of_sample], "white"),
               stringsAsFactors = FALSE),
    x_label = "Sample", y_label = marker, class_colours = palette,
    title = paste(marker, figure, "at", probe))

  ggplot2::ggsave(filename, plot = plot_object, width = 9, height = 4,
                  units = "in", dpi = as.numeric(ssEnv$plot_resolution_ppi))

  invisible(filename)
}

# One probe's row of one pivot, as SAMPLE / VALUE, or NULL when the pivot or the
# probe is not there. The transpose of .sem_probe_column_get(): same read, the
# other axis.
.sem_probe_row_get <- function(marker, figure, probe,
                               samples_selection = NULL) {

  pivot_lazy <- io_read_pivot(marker, figure, area = "PROBE", subarea = "WHOLE")
  if (is.null(pivot_lazy))
    return(NULL)

  pivot <- as.data.frame(pivot_lazy$collect())
  if (nrow(pivot) == 0L)
    return(NULL)

  # Pivot instance names are written upper case, so the comparison is made on
  # one case rather than hoping the caller used it.
  row_index <- which(toupper(as.character(pivot$AREA)) == toupper(probe))
  if (length(row_index) == 0L)
    return(NULL)
  row_index <- row_index[1]

  sample_columns <- setdiff(colnames(pivot), "AREA")
  if (length(sample_columns) == 0L)
    return(NULL)

  values <- data.frame(
    SAMPLE = sample_columns,
    VALUE = as.numeric(unlist(pivot[row_index, sample_columns, drop = TRUE],
                               use.names = FALSE)),
    stringsAsFactors = FALSE)
  values$VALUE[is.na(values$VALUE)] <- 0

  if (!is.null(samples_selection))
    values <- values[values$SAMPLE %in% as.character(samples_selection), ,
                     drop = FALSE]

  if (nrow(values) == 0L) NULL else values
}
