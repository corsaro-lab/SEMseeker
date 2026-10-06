#' Core SEMseeker pipeline (internal)
#'
#' Internal entry point used by the public \code{\link{semseeker}} dispatcher.
#' Accepts already-normalised signal data (matrix or data frame with probe-ID
#' rownames) and runs the full SEM analysis. Users should call
#' \code{\link{semseeker}} instead, which handles input normalisation,
#' M-value conversion, and tech/genome_build validation.
#'
#' @param sample_sheet Data frame (or list of data frames) with a \code{Sample_ID} column.
#' @param signal_data Methylation matrix or data frame with probe-ID rownames
#'   (or list of such objects, one per batch).
#' @param result_folder Output directory.
#'
#' @return Invisibly NULL; writes output files to \code{result_folder}.
#' @keywords internal
#' @importFrom doRNG %dorng%
sem_core <- function(sample_sheet,
  signal_data,
  result_folder) {

  # core_init_env() is called by semseeker() before we get here.
  ssEnv <- core_get_session_info()
  core_log_event("BANNER:", format(Sys.time(), "%a %b %d %X %Y"), " SEMseeker will search MARKERS for project \n in ", ssEnv$result_folderData)

  # check if the input is a list of data frames
  if(!is.list(sample_sheet) | is.data.frame(sample_sheet))
    sample_sheet <- list(sample_sheet)
  if(!is.list(signal_data) | is.data.frame(signal_data))
    signal_data <- list(signal_data)


  batch_id <- 1
  ssEnv$batch_count <- length(sample_sheet)
  ssEnv <- core_update_session_info(ssEnv)

  # C-06: write session provenance metadata (genome_build, tech, version, ...)
  total_sample_n <- sum(vapply(sample_sheet, nrow, integer(1)))
  core_session_metadata_write(result_folder, sample_n = total_sample_n)

  for(batch_id in seq_along(sample_sheet))
  {
    start_time <- Sys.time()
    ssEnv$running_batch_id <- batch_id
    ssEnv <- core_update_session_info(ssEnv)
    sample_sheet_local <- io_source_data_get(sample_sheet[[batch_id]])
    # single normalisation point for the sheet side. The signal side is
    # normalised inside sem_analyze_batch(), where the matrix is materialised
    # (keeping it there avoids holding a second reference to a multi-GB object
    # alive in this frame for the whole batch).
    sample_sheet_local <- core_normalize_sample_ids(sample_sheet_local)$sample_sheet
    utils::write.csv2(sample_sheet_local, file = io_file_path_build(ssEnv$result_folderData, paste0(batch_id,"_sample_sheet_original"),"csv",FALSE))
    sem_analyze_batch(io_source_data_get(signal_data[[batch_id]]), sample_sheet_local)
    anno_create_position_pivots(sample_sheet_local,ssEnv$keys_markers_figures)
    core_log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), "Batch Executed in:", difftime(time1 = Sys.time(), time2= start_time,units = "mins") , " minutes.")
  }

  sem_deltaX_get()
  sem_study_summary_total()
  # per-sample statistics sibling (burden + signal descriptors),
  # written next to the sample sheet and joined back on Sample_ID by
  # sem_study_summary_get().
  sem_sample_stats_build()
  anno_annotate_position_pivots()

  # Single point of sidecar materialisation.
  core_ensure_sidecars(ssEnv$result_folderData)

  core_log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), " Saving Sample Sheet with Results! ")

  core_close_env()
}
