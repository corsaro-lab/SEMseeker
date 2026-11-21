#--- sem_signal_range_values ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#' calculate the range of signal values to define the outlier
#' @param populationMatrix matrix of methylation for the population under calculation (probes × samples)
#' @param batch_id character string identifying the batch; used to name the cached parquet output file
#' @param probe_features data.frame of probe annotations (from anno_probe_features_get) with columns PROBE, CHR, START, END
#'
#' @return data.frame of per-probe thresholds with columns signal_inferior_thresholds,
#'   signal_superior_thresholds, signal_median_values, iqr, q1, q3, CHR, START, END
#'   and the annotation columns; the probe id is used for the join and not kept,
#'   because results are indexed by coordinates only
#' @importFrom doRNG %dorng%

sem_signal_range_values <- function(populationMatrix, batch_id, probe_features) {

  ssEnv <- core_get_session_info()
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Starting signal thresholds calculation.")
  thresholds_file_name <- io_file_path_build(ssEnv$result_folderData ,c(batch_id, "signal_thresholds"),"parquet")
  if(file.exists(thresholds_file_name))
  {
    result <- io_read_parquet_frame(thresholds_file_name)
    return(result)
  }
  if (sum(is.na(populationMatrix)) > 0)
  {
    msg <- paste0("ERROR:", format(Sys.time(), "%a %b %d %X %Y"), " There are missing values in the population matrix, apply the parameter inpute or remove the missing values.")
    core_log_event(msg)
    stop(msg)
  }

  # Vectorized threshold computation via matrixStats (C-level, ~10x faster
  # than row-wise apply). Scales to 28M positions (whole-genome WGBS).
  # Caller passes a probe-indexed data.frame or matrix (no PROBE column);
  # probe IDs live in rownames. as.matrix() forces numeric matrix needed by
  # matrixStats but is a no-op if input is already a matrix.
  mat <- if (is.matrix(populationMatrix)) populationMatrix else as.matrix(populationMatrix)
  iqr_times <- as.numeric(ssEnv$iqrTimes)

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " Computing thresholds for ", nrow(mat), " positions x ",
            ncol(mat), " samples (vectorized).")

  q1  <- matrixStats::rowQuantiles(mat, probs = 0.25, na.rm = TRUE)
  med <- matrixStats::rowMedians(mat, na.rm = TRUE)
  q3  <- matrixStats::rowQuantiles(mat, probs = 0.75, na.rm = TRUE)
  iqr <- q3 - q1

  # Tukey's fences, and nothing tighter. They used to be clamped to the minimum
  # and maximum the Reference observed: with a small Reference that extreme is
  # tighter than the fence on nearly every probe (all 25 KCNQ1OT1 probes with
  # six controls), so a sample barely below every control was an epimutation.
  # The clamp kept the thresholds inside [0, 1] only as a side effect.
  signal_inferior <- q1 - iqr_times * iqr
  signal_superior <- q3 + iqr_times * iqr

  # A beta value is a proportion: a fence outside [0, 1] is outside what can be
  # measured, and would inflate the high - low width DELTAR divides by. M-values
  # are not bounded and are left as they are.
  if (identical(io_signal_figure(), "BETA")) {
    signal_inferior <- pmax(signal_inferior, 0)
    signal_superior <- pmin(signal_superior, 1)
  }

  result <- data.frame(
    signal_inferior_thresholds = signal_inferior,
    signal_superior_thresholds = signal_superior,
    signal_median_values       = med,
    iqr                        = iqr,
    q1                         = q1,
    q3                         = q3
  )
  rm(mat, q1, med, q3, iqr, signal_inferior, signal_superior)
  gc()

  #
  colnames(result) <- c("signal_inferior_thresholds","signal_superior_thresholds", "signal_median_values","iqr","q1","q3")
  result$PROBE <- rownames(populationMatrix)
  if(nrow(result) != nrow(populationMatrix))
    stop("I'M STOPPING HERE, No thresholds defined for the population.")

  # The thresholds take the probe coordinates and annotation by joining on
  # PROBE: the annotation columns first, then the thresholds, sorted on CHR,
  # START and END with CHR as text. Factors are passed as text, because the
  # data engine would sort a factor by its levels and not lexically.
  as_text <- function(frame) {
    # frame: a data.frame; returns it with every factor column as character
    factor_columns <- vapply(frame, is.factor, logical(1))
    frame[factor_columns] <- lapply(frame[factor_columns], as.character)
    frame
  }
  # PROBE_WHOLE is the canonical (AREA, SUBAREA) name of the probe-level
  # baseline, kept on probe_features for the association lookups; it is a
  # 1:1 alias of PROBE and is not persisted with the thresholds.
  probe_features <- probe_features[, setdiff(colnames(probe_features), "PROBE_WHOLE"), drop = FALSE]
  connection <- io_engine_connection()
  annotation_lazy <- io_register_frame(as_text(probe_features))
  thresholds_lazy <- io_register_frame(as_text(result))
  # The probe id joins the two and is then dropped: what SEMseeker writes is
  # indexed by coordinates, never by the name a manufacturer gave the probe.
  result <- DBI::dbGetQuery(connection, sprintf('SELECT * EXCLUDE (PROBE) FROM (%s) AS annotation JOIN (%s) AS thresholds USING (PROBE) ORDER BY CHR, START, "END"',
                                                dbplyr::sql_render(annotation_lazy), dbplyr::sql_render(thresholds_lazy)))
  io_write_frame_parquet(result, thresholds_file_name)

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Thresholds defined for: ", nrow(result), " probe_features.")
  gc()
  return(result)
}
