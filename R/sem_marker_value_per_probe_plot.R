#' Marker value of one sample across probes
#'
#' One point per probe for a single sample: the value of the marker on the y
#' axis, the probe on the x, and the colour saying whether that probe is a
#' hyper-methylation event, a hypo-methylation event or neither for this sample.
#'
#' @section Through 0.99.6 this was `sem_manhattan_plot_marker_per_sample()`:
#' That function was exported, documented, and could not run. Measured on
#' 2026-10-04:
#'
#' * it called `core_init_env()` from inside a chart, asking the session for a
#'   figure named `BOTH`, which the figure vocabulary refuses: the call ended in
#'   `Unrecognised argument(s) ... figures = BOTH` before anything was read. A
#'   chart does not re-initialise the session, and this one takes the session as
#'   it finds it, which removes the failure rather than working around it.
#' * past that it read pivots through `io_pivot_file_name()`, the `.csv` variant.
#'   Nothing in the package has written a CSV pivot since the move to parquet, so
#'   `file.exists()` was false for every marker, every marker was dropped, and
#'   the function would have returned having drawn nothing and said nothing.
#' * it dropped markers with `tempKeys <- tempKeys[-k]` while indexing a
#'   different vector by `k`, so after the first removal it removed the wrong
#'   entries.
#' * it read `ssEnv$color_palette` on its third line, before assigning `ssEnv` on
#'   its eighth. The name resolved to a stale object reachable from the
#'   namespace, so the colours came from a frozen snapshot and not from the
#'   session the caller had configured.
#' * it passed `dpi = as.numeric(ssEnv$plot_resolution)`, and that field holds
#'   `"print"`. The resolution was `NA` on every save. The field with a number in
#'   it is `plot_resolution_ppi`.
#'
#' @section The name says the axis:
#' `per_probe` means one point per probe, which is the x axis, and the sample is
#' what the chart holds fixed - the opposite of what the old name's `per_sample`
#' meant. The convention is the one `sem_marker_value_per_area_plot()` already
#' uses, and the argument that is held fixed is named in the signature, so the
#' reading is in the call and not only in the name.
#'
#' @section Which probes, and why not an index range:
#' The old signature took `probes_range = 1000:2000`, positions in the row order
#' of a file. That is not a selection anyone can reproduce or report: it changes
#' with the order the pivot happens to be written in. `probes_selection` names
#' the probes instead. `NULL` draws all of them, which is the normal case for a
#' Manhattan-style chart, and the drawing drops the axis text once the names stop
#' being readable.
#'
#' @section How the direction of each event is decided:
#' The class of a probe comes from the `HYPER` and `HYPO` pivots, which are what
#' define a direction: a non-zero value for this sample in `HYPER` is a
#' hyper-methylation event at that probe. The value drawn comes from `figure`,
#' so the figure chooses what is on the y axis and the direction is read from the
#' data either way. Where neither pivot carries a non-zero value the probe is a
#' non-outlier and its value is drawn at zero, because a marker value with no
#' event behind it is not an event.
#'
#' @param marker character. Marker name, e.g. `"MUTATIONS"`, `"DELTAS"`.
#' @param figure character. The figure whose values go on the y axis: `"HYPER"`,
#'   `"HYPO"` or `"BOTH"`.
#' @param sample_name character. The sample held fixed. Must be a column of the
#'   pivot.
#' @param probes_selection character vector of probe names, or `NULL` for all.
#' @param overwrite logical. Redraw when the file exists.
#'
#' @return Invisibly the path of the file written, or `NULL` when there was
#'   nothing to draw.
#'
#' @examples
#' # One sample held fixed, a point per probe. The chart reads the probe pivots
#' # of the marker, so it needs a result folder from a completed run.
#' probes_of_interest <- c("cg11680158", "cg00000029")
#' \dontrun{
#' sem_marker_value_per_probe_plot(
#'   marker = "DELTAS", figure = "BOTH", sample_name = "CASE_1",
#'   probes_selection = probes_of_interest
#' )
#' }
#' @export
sem_marker_value_per_probe_plot <- function(marker, figure, sample_name,
                                           probes_selection = NULL,
                                           overwrite = FALSE) {

  ssEnv <- core_get_session_info()

  chart_folder <- io_dir_check_and_create(ssEnv$result_folderChart,
                                          "MARKER_VALUE_PER_PROBE")
  filename <- io_file_path_build(
    chart_folder,
    toupper(c("VALUE_PER_PROBE", marker, figure,
              core_name_cleaning(sample_name))),
    ssEnv$plot_format)

  if (file.exists(filename) && !overwrite) {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " chart already exists, skipped: ", filename)
    return(invisible(filename))
  }

  values <- .sem_probe_column_get(marker, figure, sample_name, probes_selection)
  if (is.null(values)) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " no ", marker, "/", figure, " probe values for sample ",
              sample_name, ": nothing drawn.")
    return(invisible(NULL))
  }

  hyper <- .sem_probe_column_get(marker, "HYPER", sample_name, values$AREA)
  hypo  <- .sem_probe_column_get(marker, "HYPO",  sample_name, values$AREA)

  class_of_probe <- rep("Non-outlier", nrow(values))
  if (!is.null(hyper))
    class_of_probe[values$AREA %in% hyper$AREA[hyper$VALUE != 0]] <- "Hyper"
  if (!is.null(hypo))
    class_of_probe[values$AREA %in% hypo$AREA[hypo$VALUE != 0]] <- "Hypo"

  palette <- .sem_outlier_palette(ssEnv)
  # A marker value with no event behind it is not an event: it is drawn at zero
  # rather than being removed, so the width of the axis keeps meaning the extent
  # that was examined.
  drawn_value <- ifelse(class_of_probe == "Non-outlier", 0, values$VALUE)

  plot_object <- .plot_lollipop_build(
    data.frame(X = values$AREA, VALUE = drawn_value, CLASS = class_of_probe,
               BASELINE = 0,
               SEGMENT_COLOUR = ifelse(class_of_probe == "Non-outlier",
                                       "white", palette[class_of_probe]),
               stringsAsFactors = FALSE),
    x_label = "Probe", y_label = marker, class_colours = palette,
    title = paste(marker, figure, "of", sample_name))

  ggplot2::ggsave(filename, plot = plot_object, width = 9, height = 4,
                  units = "in", dpi = as.numeric(ssEnv$plot_resolution_ppi))

  invisible(filename)
}

# One sample's column of one probe pivot, as AREA / VALUE, or NULL when the
# pivot or the sample is not there. Reading the column rather than the whole
# frame is what keeps a per-sample chart from materialising a pivot that can
# hold hundreds of samples across nearly a million probes.
.sem_probe_column_get <- function(marker, figure, sample_name,
                                  probes_selection = NULL) {

  pivot_lazy <- io_read_pivot(marker, figure, area = "PROBE", subarea = "WHOLE")
  if (is.null(pivot_lazy))
    return(NULL)

  pivot <- as.data.frame(pivot_lazy$collect())
  if (nrow(pivot) == 0L || !(sample_name %in% colnames(pivot)))
    return(NULL)

  values <- data.frame(AREA = as.character(pivot$AREA),
                       VALUE = as.numeric(pivot[[sample_name]]),
                       stringsAsFactors = FALSE)
  values$VALUE[is.na(values$VALUE)] <- 0

  if (!is.null(probes_selection))
    values <- values[values$AREA %in% as.character(probes_selection), ,
                     drop = FALSE]

  if (nrow(values) == 0L) NULL else values
}

# The four classes an outlier chart draws, and their colours, taken from the
# session palette so that a study's charts agree with each other. The names are
# the contract between the doors and the drawing.
.sem_outlier_palette <- function(ssEnv) {
  c(Hyper         = ssEnv$color_palette_darker[1],
    Hypo          = ssEnv$color_palette_darker[2],
    "Non-outlier" = "grey",
    Reference     = ssEnv$color_palette[3])
}
