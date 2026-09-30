#' Statistics emitted for a signal scope
#'
#' Single source of truth for which statistics a scope carries:
#' [util_signal_descriptors()] builds its output on this vector, so the computed
#' keys and the declared vocabulary cannot drift apart.
#'
#' It used to serve a second side - the producer of `SAMPLE_STATS_RESULT.csv`
#' composed that file's column names from it via [io_stat_colname()]. AI-255
#' removed the file, and with it that call site; this vector now has one
#' consumer, and [io_stat_colname()] has none.
#'
#' MODELOW / MODEHIGH only make sense on a bounded, bimodal scale (beta in
#' [0,1], the two peaks around 0 and 1). On M-values the distribution is
#' roughly gaussian and the two-mode split carries no meaning, so the columns
#' are omitted rather than filled with a meaningless number.
#'
#' **`N_PROBES` is not here (AI-255).** This vector is the vocabulary of columns
#' an artefact carries, and no artefact carries `N_PROBES` any more: it is a
#' property of the imputation - how many positions of that sample survived the
#' treatment of missing values - and it travels on `SAMPLE_SHEET_RESULT`.
#' [util_signal_descriptors()] still reports it, because a caller reducing a
#' vector legitimately wants to know how many finite values it had, but that is a
#' description of the input, not a column of the output.
#'
#' @param beta logical. TRUE when the signal is on the [0,1] beta scale.
#' @return character vector of statistic names.
#' @keywords internal
#' @noRd
io_signal_stats <- function(beta = TRUE) {
  base <- c("MEDIAN", "MEAN", "VARIANCE", "IQR")
  if (isTRUE(beta)) c(base, "MODELOW", "MODEHIGH") else base
}
