#' sem_deltar_single_sample
#'
#' @param values data.frame of methylation values with columns CHR, START, END and signal value in column 4
#' @param thresholds data.frame of signal thresholds (from sem_signal_range_values) with columns CHR, START, END,
#'   signal_superior_thresholds, signal_inferior_thresholds, signal_median_values
#' @param sample_detail named list/row with at least Sample_ID and Sample_Group fields
#' @return invisibly NULL; HYPER and HYPO relative-delta results are written as bedgraph.gz files.
#'   Only positions present in both \code{values} and \code{thresholds} are processed
#'   (inner join on CHR, START, END via \code{util_join_values_to_thresholds}).
#'   Returns \code{invisible(NULL)} early if there is no positional overlap.
#'
sem_deltar_single_sample <- function(values, thresholds, sample_detail) {

  ssEnv <- core_get_session_info()

  # Polars inner join on (CHR, START, END) - replaces sort-then-positional-zip.
  # util_join_values_to_thresholds() handles type coercion and is shared with
  # sem_mutations_get() and sem_delta_single_sample().
  joined <- util_join_values_to_thresholds(values, thresholds)

  if (nrow(joined) == 0L) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
      " [sem_deltar_single_sample] No overlapping positions - skipping deltar",
      " for sample=", sample_detail$Sample_ID)
    return(invisible(NULL))
  }

  high_thresholds <- joined$signal_superior_thresholds
  low_thresholds  <- joined$signal_inferior_thresholds

  if (any(is.na(high_thresholds)) || any(is.na(low_thresholds))) {
    # NAs present - logged but not fatal; downstream subset(DELTA > 0) handles them
  }

  dividend <- high_thresholds - low_thresholds

  if (any(is.na(dividend)) || any(dividend < 0))
    stop("ERROR: I'm stopping here - dividend has NA or negative values!")

  # Avoid division by zero
  dividend[dividend == 0] <- 0.000000001

  ### get deltar HYPER ###########################################################
  deltar_hyper <- data.frame(
    CHR   = joined$CHR,
    START = joined$START,
    END   = joined$END,
    DELTA = (joined$VALUE - high_thresholds) / dividend,
    stringsAsFactors = FALSE
  )
  deltar_hyper <- anno_sort_by_chr_and_start(deltar_hyper)
  deltar_hyper <- subset(deltar_hyper, deltar_hyper$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- io_dir_check_and_create(ssEnv$result_folderData,
    c(as.character(sample_detail$Sample_Group), "DELTAR_HYPER"))
  io_dump_sample_as_bed_file(
    data_to_dump = deltar_hyper,
    fileName     = io_file_path_build(folder_to_save,
                     c(as.character(sample_detail$Sample_ID), "DELTAR", "HYPER"),
                     "bedgraph", add_gz = TRUE)
  )

  ### get deltar HYPO ############################################################
  deltar_hypo <- data.frame(
    CHR   = joined$CHR,
    START = joined$START,
    END   = joined$END,
    DELTA = (low_thresholds - joined$VALUE) / dividend,
    stringsAsFactors = FALSE
  )
  deltar_hypo <- anno_sort_by_chr_and_start(deltar_hypo)
  deltar_hypo <- subset(deltar_hypo, deltar_hypo$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- io_dir_check_and_create(ssEnv$result_folderData,
    c(as.character(sample_detail$Sample_Group), "DELTAR_HYPO"))
  io_dump_sample_as_bed_file(
    data_to_dump = deltar_hypo,
    fileName     = io_file_path_build(folder_to_save,
                     c(as.character(sample_detail$Sample_ID), "DELTAR", "HYPO"),
                     "bedgraph", add_gz = TRUE)
  )

  deltar_to_check <- c(deltar_hypo$DELTA, deltar_hyper$DELTA)
  if (length(deltar_to_check) > 0)
    if (min(deltar_to_check) < 0) {
      core_log_event(min(deltar_hypo$DELTA))
      core_log_event(min(deltar_hyper$DELTA))
      stop("ERROR: I'm stopping here - deltar have negative values!")
    }
}
