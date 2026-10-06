#' Detect the Illumina methylation array technology from a signal matrix
#'
#' Identifies whether a methylation signal matrix originates from a 27k, 450k,
#' or EPIC 850k array (or WGBS data).  Detection runs in priority order:
#'
#' \enumerate{
#'   \item \strong{Annotation-package overlap} - probe IDs are matched against
#'     each installed \code{IlluminaHumanMethylation*anno} package; the
#'     technology with the most matching probes wins.  Accurate for any subset
#'     size (e.g. 20 k probes out of 866 k).
#'   \item \strong{Row-count heuristics} - last resort for WGBS or unknown
#'     probe naming schemes.
#' }
#'
#' The detected technology and beta/M-value flag are written to the session
#' environment.
#'
#' @param signal_data A numeric matrix or \code{data.frame} with CpG probes as
#'   rows and samples as columns.  Row names must be probe identifiers (e.g.
#'   \code{cg00000029}) unless a \code{PROBE} column is present.
#'
#' @return The updated session environment (\code{ssEnv}), invisibly.  The
#'   detected technology is accessible via \code{core_get_session_info()$tech}.
#'
core_get_meth_tech <- function(signal_data) {

 
  ssEnv    <- core_get_session_info()
  n_probes <- nrow(signal_data)

  # If technology was explicitly declared by the user in core_init_env(), respect it.
  # This is required for LONGREAD (indistinguishable from WGBS by probe-ID pattern)
  # and useful when the user wants to override heuristics.
  .declared_techs <- c("K850", "K450", "K27", "WGBS", "LONGREAD")
  if (!is.null(ssEnv$tech) && ssEnv$tech %in% .declared_techs) {
    core_log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
              "technology pre-declared as '", ssEnv$tech,
              "'; skipping auto-detection.")
    # Still detect beta vs M-values
    exclude_cols <- c("AREA","PROBE","CHR","START","END","K27","K450","K850","k27","k450","k850")
    signal_cols  <- signal_data[seq_len(min(10000L, nrow(signal_data))),!colnames(signal_data) %in% exclude_cols, drop = FALSE]
    max_data   <- max(abs(c(max(signal_cols, na.rm = TRUE),min(signal_cols, na.rm = TRUE))))
    ssEnv$beta <- max_data <= 1
    ssEnv$probes_count <- n_probes
    core_update_session_info(ssEnv)
    return(ssEnv)
  }

  # Informational row-count hints
  if (n_probes == 485512)
    core_log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
              "probe count matches 450k dataset.")
  if (n_probes == 27578)
    core_log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
              "probe count matches 27k dataset.")
  if (n_probes == 866562)
    core_log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
              "probe count matches EPIC 850k dataset.")
  if (n_probes > 866562)
    core_log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
              "probe count exceeds EPIC 850k - treating as WGBS.")

  # Resolve probe identifiers
  probe_ids <- if ("PROBE" %in% colnames(signal_data))
    signal_data$PROBE
  else
    rownames(signal_data)

  tech <- ""
  msg  <- ""

  # ---- Step 1: annotation-package overlap (most accurate) ----
  tech <- .anno_detect_tech_from_anno(probe_ids)
  if (tech != "")
    msg <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
                 "technology identified as", tech,
                 "via annotation-package probe-ID overlap.")

  # ---- Step 2: row-count heuristics (last resort) ----
  if (tech == "") {
    if (n_probes > 866562 ||
        (length(probe_ids) > 0L &&
         all(grepl("_", probe_ids[seq_len(min(1000L, length(probe_ids)))])))) {
      tech <- "WGBS"
      msg  <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
                    "dataset identified as WGBS (row count / probe ID pattern).")
    } else if (n_probes <= 30000L) {
      tech <- "K27"
      msg  <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
                    "dataset identified as 27k array (row-count heuristic).")
    } else if (n_probes <= 510000L) {
      tech <- "K450"
      msg  <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
                    "dataset identified as 450k array (row-count heuristic).")
    } else {
      tech <- "K850"
      msg  <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
                    "dataset identified as EPIC 850k array (row-count heuristic).")
    }
  }

  if (tech == "") {
    msg <- paste("ERROR:", format(Sys.time(), "%a %b %d %X %Y"),
                 "could not determine array technology.")
    core_log_event(msg)
    stop(msg)
  }

  core_log_event(msg)
  ssEnv$tech <- tech

  # ---- Detect beta vs M-values ----
  exclude_cols <- c("AREA", "PROBE", "CHR", "START", "END", "K27", "K450", "K850", "k27", "k450", "k850")
  signal_cols  <- signal_data[
    seq_len(min(10000L, nrow(signal_data))),
    !colnames(signal_data) %in% exclude_cols,
    drop = FALSE
  ]
  max_data   <- max(abs(c(max(signal_cols, na.rm = TRUE),
                          min(signal_cols, na.rm = TRUE))))
  ssEnv$beta <- max_data <= 1

  core_log_event(if (ssEnv$beta)
    paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"), "values are beta (0-1).")
  else
    paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"), "values appear to be M-values."))

  # the figure of SIGNAL is the scale, and the scale is only knowable
  # here - the keys were built before any data was read.
  ssEnv <- util_keys_signal_figure_refresh(ssEnv)

  ssEnv$probes_count <- n_probes
  core_update_session_info(ssEnv)
  core_log_event("JOURNAL: array technology set to:", ssEnv$tech)
  return(ssEnv)
}
