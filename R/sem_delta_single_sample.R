#' sem_delta_single_sample
#'
#' @param values data.frame of methylation values with columns CHR, START, END and signal value in column 4
#' @param thresholds data.frame of signal thresholds (from sem_signal_range_values) with columns CHR, START, END,
#'   signal_superior_thresholds, signal_inferior_thresholds, signal_median_values
#' @param sample_detail named list/row with at least Sample_ID and Sample_Group fields
#' @return invisibly NULL; HYPER and HYPO delta results are written as bedgraph.gz files.
#'   Only positions present in both \code{values} and \code{thresholds} are processed
#'   (inner join on CHR, START, END via \code{util_join_values_to_thresholds}).
#'   Returns \code{invisible(NULL)} early if there is no positional overlap.
#'
sem_delta_single_sample <- function(values, thresholds, sample_detail) {

  ssEnv <- core_get_session_info()

  # Polars inner join on (CHR, START, END) - replaces sort-then-positional-zip.
  # util_join_values_to_thresholds() handles type coercion and is shared with
  # sem_mutations_get() and sem_deltar_single_sample().
  joined <- util_join_values_to_thresholds(values, thresholds)

  if (nrow(joined) == 0L) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
      " [sem_delta_single_sample] No overlapping positions - skipping delta",
      " for sample=", sample_detail$Sample_ID)
    return(invisible(NULL))
  }

  high_thresholds <- joined$signal_superior_thresholds
  low_thresholds  <- joined$signal_inferior_thresholds

  if (any(high_thresholds < low_thresholds))
    stop("ERROR: I'm stopping here - some high thresholds have values less than low thresholds!")

  ### get deltas HYPER ###########################################################
  deltas_hyper <- data.frame(
    CHR   = joined$CHR,
    START = joined$START,
    END   = joined$END,
    DELTA = joined$VALUE - high_thresholds,
    stringsAsFactors = FALSE
  )
  deltas_hyper <- anno_sort_by_chr_and_start(deltas_hyper)
  deltas_hyper <- subset(deltas_hyper, deltas_hyper$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- io_dir_check_and_create(ssEnv$result_folderData,
    c(as.character(sample_detail$Sample_Group), "DELTAS_HYPER"))
  io_dump_sample_as_bed_file(
    data_to_dump = deltas_hyper,
    fileName     = io_file_path_build(folder_to_save,
                     c(as.character(sample_detail$Sample_ID), "DELTAS", "HYPER"),
                     "bedgraph", add_gz = TRUE)
  )

  ### get deltas HYPO ############################################################
  deltas_hypo <- data.frame(
    CHR   = joined$CHR,
    START = joined$START,
    END   = joined$END,
    DELTA = low_thresholds - joined$VALUE,
    stringsAsFactors = FALSE
  )
  deltas_hypo <- anno_sort_by_chr_and_start(deltas_hypo)
  deltas_hypo <- subset(deltas_hypo, deltas_hypo$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- io_dir_check_and_create(ssEnv$result_folderData,
    c(as.character(sample_detail$Sample_Group), "DELTAS_HYPO"))
  io_dump_sample_as_bed_file(
    data_to_dump = deltas_hypo,
    fileName     = io_file_path_build(folder_to_save,
                     c(as.character(sample_detail$Sample_ID), "DELTAS", "HYPO"),
                     "bedgraph", add_gz = TRUE)
  )

  deltas_to_check <- c(deltas_hypo$DELTA, deltas_hyper$DELTA)
  if (length(deltas_to_check) > 0)
    if (min(deltas_to_check) < 0)
      stop("ERROR: I'm stopping here - deltas have negative values!")
}
