#' Signal of one probe against the reference envelope, across samples
#'
#' One point per sample for a single probe: the raw signal on the y axis, the
#' sample on the x, the reference envelope of that probe drawn as horizontal
#' lines, and each excursion labelled with how far outside it fell.
#'
#' This is the chart that says *why* a sample is an event at a probe, where
#' [sem_marker_value_per_sample_plot()] says *that* it is one. They are two
#' questions and therefore two doors, over the same layout.
#'
#' @section The labels that had never been drawn:
#' Through 0.99.6 the delta annotation of this chart was built with
#' `ifelse(delta_label_value, "")`, a two-argument call. It raises nothing: the
#' character vector coerces to `NA`, `ifelse` returns `NA`, and `geom_text` draws
#' nothing. So the distance of each excursion - the only quantitative content of
#' the chart, and the reason it carries a `label_font_size` parameter - had never
#' appeared on a single image, and no error anywhere said so. A defect that
#' raises is found the first time the code runs; this one needed someone to
#' notice that a label was missing from a picture.
#'
#' @section Where the segment starts:
#' At the threshold the point crossed, not at zero. The length of the segment is
#' then the excursion itself, which is what the chart is about, and it is the
#' same number the label states. A sample inside the envelope gets a segment of
#' zero length rather than a line down to the axis, because it has no excursion
#' to show.
#'
#' @param probe character. The probe held fixed, e.g. `"cg11680158"`.
#' @param sample_sheet data.frame carrying `Sample_ID` and `Sample_Group`.
#' @param samples_selection character vector of sample names, or `NULL` for all.
#' @param show_labels logical. Draw the names of the envelope lines. The deltas
#'   are drawn either way: they are the content, not a decoration.
#' @param label_font_size numeric. Size of the delta labels.
#' @param overwrite logical. Redraw when the file exists.
#'
#' @return Invisibly the path of the file written, or `NULL` when the signal
#'   pivot or the thresholds of that probe are not there.
#'
#' @examples
#' # The raw signal of one probe against the envelope computed from the
#' # reference samples. The envelope is read from the study, never recomputed:
#' # a threshold recomputed on the samples at hand is a different threshold.
#' sample_sheet <- data.frame(
#'   Sample_ID    = c("CASE_1", "CASE_2", "REF_1"),
#'   Sample_Group = c("Case", "Case", "Reference")
#' )
#' \dontrun{
#' sem_signal_threshold_per_sample_plot(
#'   probe = "cg11680158", sample_sheet = sample_sheet
#' )
#' }
#' @export
sem_signal_threshold_per_sample_plot <- function(probe, sample_sheet,
                                                samples_selection = NULL,
                                                show_labels = TRUE,
                                                label_font_size = 3,
                                                overwrite = FALSE) {

  ssEnv <- core_get_session_info()

  chart_folder <- io_dir_check_and_create(ssEnv$result_folderChart,
                                          "SIGNAL_THRESHOLD_PER_SAMPLE")
  filename <- io_file_path_build(
    chart_folder,
    toupper(c("SIGNAL_THRESHOLD", core_name_cleaning(probe))),
    ssEnv$plot_format)

  if (file.exists(filename) && !overwrite) {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " chart already exists, skipped: ", filename)
    return(invisible(filename))
  }

  signal <- .sem_probe_row_get("SIGNAL", io_signal_figure(), probe,
                               samples_selection)
  if (is.null(signal)) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " no signal values at probe ", probe, ": nothing drawn.")
    return(invisible(NULL))
  }

  envelope <- .sem_probe_thresholds_get(probe)
  if (is.null(envelope)) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " no signal thresholds for probe ", probe, ": nothing drawn.")
    return(invisible(NULL))
  }

  plot_object <- .sem_signal_threshold_plot_build(
    signal, envelope, sample_sheet, probe, .sem_outlier_palette(ssEnv),
    show_labels = show_labels, label_font_size = label_font_size)

  ggplot2::ggsave(filename, plot = plot_object, width = 9, height = 5,
                  units = "in", dpi = as.numeric(ssEnv$plot_resolution_ppi))

  invisible(filename)
}

# The reference envelope of one probe, or NULL. The thresholds are written once
# per study by the population analysis, so this is a lookup and not a
# computation: a chart must never recompute a threshold, because a threshold
# recomputed on a different sample set is a different threshold.
.sem_probe_thresholds_get <- function(probe) {

  ssEnv <- core_get_session_info()
  thresholds_file <- io_file_path_build(ssEnv$result_folderData,
                                        "1_signal_thresholds", "parquet")
  if (!file.exists(thresholds_file))
    return(NULL)

  thresholds <- as.data.frame(polars::pl$read_parquet(thresholds_file))
  row_index <- which(toupper(as.character(thresholds$PROBE)) == toupper(probe))
  if (length(row_index) == 0L)
    return(NULL)

  row <- thresholds[row_index[1], , drop = FALSE]
  list(superior = as.numeric(row$signal_superior_thresholds),
       inferior = as.numeric(row$signal_inferior_thresholds),
       median   = as.numeric(row$signal_median_values),
       q1       = as.numeric(row$q1),
       q3       = as.numeric(row$q3))
}

# The drawing of the threshold chart, separated from the door that reads and
# saves, so a test can assert on the layers rather than on the existence of a
# PNG.
.sem_signal_threshold_plot_build <- function(signal, envelope, sample_sheet,
                                            probe, palette,
                                            show_labels = TRUE,
                                            label_font_size = 3) {

  class_of_sample <- ifelse(signal$VALUE > envelope$superior, "Hyper",
                      ifelse(signal$VALUE < envelope$inferior, "Hypo",
                             "Non-outlier"))
  references <- as.character(
    sample_sheet$Sample_ID[sample_sheet$Sample_Group == "Reference"])
  class_of_sample[signal$SAMPLE %in% references] <- "Reference"

  is_excursion <- class_of_sample %in% c("Hyper", "Hypo")

  # The segment starts at the threshold that was crossed, so its length IS the
  # excursion, and a sample inside the envelope gets no segment at all.
  baseline <- ifelse(class_of_sample == "Hyper", envelope$superior,
               ifelse(class_of_sample == "Hypo", envelope$inferior,
                      signal$VALUE))

  plot_data <- data.frame(
    X = signal$SAMPLE, VALUE = signal$VALUE, CLASS = class_of_sample,
    BASELINE = baseline,
    SEGMENT_COLOUR = ifelse(is_excursion, palette[class_of_sample], "white"),
    stringsAsFactors = FALSE)

  reference_lines <- data.frame(
    VALUE = c(envelope$median, envelope$superior, envelope$inferior,
              envelope$q1, envelope$q3),
    LABEL = if (show_labels)
      c("Median", "Upper Limit", "Lower Limit", "Q1", "Q3")
    else
      rep("", 5L),
    COLOUR = c("black", "red", "red", "red", "red"),
    LINETYPE = "dashed",
    stringsAsFactors = FALSE)
  reference_lines <- reference_lines[is.finite(reference_lines$VALUE), ,
                                     drop = FALSE]

  # The label is the distance outside the envelope, prefixed with a delta.
  # The character is written as a \u escape and not literally: R CMD check
  # tolerates non-ASCII in a comment and refuses it in code, and this one is in
  # a string.
  # excursions only: writing one for every sample would annotate "0" across the
  # whole axis.
  distance <- ifelse(class_of_sample == "Hyper",
                     signal$VALUE - envelope$superior,
                     envelope$inferior - signal$VALUE)
  annotations <- data.frame(
    X = signal$SAMPLE[is_excursion],
    VALUE = signal$VALUE[is_excursion],
    LABEL = paste0("\u03B4 ", format(distance[is_excursion],
                                     scientific = TRUE, digits = 2)),
    COLOUR = palette[class_of_sample[is_excursion]],
    stringsAsFactors = FALSE)

  .plot_lollipop_build(
    plot_data, x_label = "Sample", y_label = "Signal",
    class_colours = palette,
    reference_lines = reference_lines,
    annotations = if (nrow(annotations) > 0L) annotations else NULL,
    label_font_size = label_font_size,
    title = paste("Signal against the reference envelope at", probe))
}
