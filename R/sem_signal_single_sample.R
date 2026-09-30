#' sem_signal_single_sample
#'
#' @param values signal values
#' @param sample_detail detais of sample
#' @param probe_features annotation probe
#'
#' @return signal mean
#'
sem_signal_single_sample <- function(values,sample_detail,probe_features)
{
  ssEnv <- core_get_session_info()

  # AI-248: the per-sample folder is named after the (marker, figure) key, and
  # the figure of SIGNAL is the scale - it must agree with io_bed_file_name()
  # and with io_list_bed_files_for_marker_figure(), which rebuild the same name.
  folder_to_save <- io_dir_check_and_create(ssEnv$result_folderData, c(sample_detail$Sample_Group ,paste0("SIGNAL","_", io_signal_figure(), sep = "")))
  # DEBUG (2026-06-09): right before the data.frame() that has
  # been exploding with "arguments imply differing number of rows" since
  # v35. Inspect:
  #   length(values), nrow(probe_features), class(values)
  #   head(names(values)), head(probe_features$PROBE)
  #   is.null(names(values)), all(probe_features$PROBE %in% names(values))
 
  signal_values_annotated <- data.frame(as.data.frame(probe_features), "VALUE" = values, row.names = probe_features$PROBE)[, c("CHR", "START", "END","VALUE")]
  io_dump_sample_as_bed_file(
    data_to_dump = signal_values_annotated,
    fileName = io_file_path_build(baseFolder =  folder_to_save, detailsFilename =  c(sample_detail$Sample_ID,"SIGNAL", io_signal_figure()), extension = "bedgraph", add_gz=TRUE)
  )
}
