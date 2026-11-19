
#' Convert a POSITION-keyed signal pivot into a PROBE-keyed lazy table.
#'
#' The SEM step works on probes, while the SIGNAL pivot is keyed by position.
#' The result is a lazy table of the data engine, not an R data.frame:
#' callers collect only the columns they need (the reference samples, then
#' one population at a time), which keeps a 367k x 4 013 matrix out of R.
#'
#' Rows follow the order of the probe annotation, as the polars version did.
#' The annotation is reduced to the probes of the detected technology and to
#' distinct rows before the join; the pivot has one row per position, so this
#' gives the same rows as removing duplicates after the join.
#'
#' The annotation is registered as a view for the life of the connection,
#' because the lazy table returned refers to it.
#'
#' @param signal_data a tbl_lazy of the data engine (for instance from
#'   io_scan_parquet()) or an R data.frame, with CHR, START, END and one
#'   column per sample
#'
#' @return A tbl_lazy with the PROBE column and one column per sample, in
#'   the order of the probe annotation.
#' @keywords internal
#' @noRd
anno_position_pivot_to_probe <- function(signal_data)
{
  ssEnv     <- core_get_session_info()
  tech_col  <- ssEnv$tech
  probe_features <- anno_probe_annotation_build(ssEnv$tech)[, c("PROBE", "CHR", "START", "END", tech_col)]
  in_tech <- probe_features[[tech_col]]
  probe_features <- unique(probe_features[!is.na(in_tech) & in_tech, c("PROBE", "CHR", "START", "END")])
  probe_features <- data.frame(PROBE = as.character(probe_features$PROBE),
                               CHR   = as.character(probe_features$CHR),
                               START = probe_features$START,
                               END   = probe_features$END,
                               ANNOTATION_ORDER = seq_len(nrow(probe_features)),
                               stringsAsFactors = FALSE)
  annotation <- io_register_frame(probe_features)

  signal_lazy <- if (inherits(signal_data, "tbl_lazy")) signal_data else io_register_frame(signal_data)
  connection <- dbplyr::remote_con(signal_lazy)
  query_text <- sprintf('
    SELECT annotation.PROBE, signal.* EXCLUDE (CHR, START, "END")
    FROM (%s) AS annotation
    JOIN (%s) AS signal
      ON annotation.CHR = CAST(signal.CHR AS VARCHAR) AND annotation.START = signal.START AND annotation."END" = signal."END"
    ORDER BY annotation.ANNOTATION_ORDER',
    dbplyr::sql_render(annotation), dbplyr::sql_render(signal_lazy))
  dplyr::tbl(connection, dbplyr::sql(query_text))
}
