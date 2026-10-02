#' The FIGURE of the SIGNAL marker for this session (internal)
#'
#' AI-248. `SIGNAL` used to carry the figure `MEAN`, which was a **placeholder**
#' that made the key unique - not a statement about the data. It described (badly)
#' a way of aggregating the signal, and now that aggregation has an axis of its
#' own the placeholder has to go.
#'
#' The figure of `SIGNAL` is the **scale** of the value: `BETA` for a bounded
#' proportion in [0,1] - an array beta value, or a methylation fraction from
#' reads, which is the same scale - and `MVALUE` for the logit-transformed one.
#' The scale is already detected once per session by `core_get_meth_tech()`.
#'
#' The technology (`ssEnv$tech`) deliberately does **not** enter here: WGBS,
#' ONT and PacBio produce fractions in [0,1] and are `BETA` like an array is.
#' The package compares array and sequencing on purpose - see
#' `test-cross-format-convergence.R` - and giving the same scale two different
#' labels would misalign exactly those comparisons. Provenance is recorded
#' elsewhere (`ssEnv$tech`, the `TECH` column stamped on results).
#'
#' @param beta logical, defaults to the scale of the current session.
#' @return `"BETA"` or `"MVALUE"`.
#' @keywords internal
#' @noRd
io_signal_figure <- function(beta = NULL) {
  if (is.null(beta)) {
    ssEnv <- core_get_session_info()
    beta <- ssEnv$beta
  }
  # An unset scale means the session has not seen the data yet; the bounded
  # scale is the package default (core_get_meth_tech sets it on first read).
  if (is.null(beta) || length(beta) != 1L || is.na(beta))
    beta <- TRUE

  if (isTRUE(beta)) "BETA" else "MVALUE"
}
