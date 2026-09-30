#' Realign the SIGNAL keys to the scale of the data (internal)
#'
#' AI-248. The figure of `SIGNAL` is the scale of the value (`BETA` / `MVALUE`,
#' see [io_signal_figure()]), but the keys are built by `util_keys_create()`
#' inside `core_init_env()`, before any data has been read - at that point the
#' scale is unknown and the bounded one is assumed.
#'
#' `core_get_meth_tech()` determines the scale on the first batch and calls this
#' to realign the stored keys. Without it an M-value run would carry `SIGNAL`
#' keys labelled `BETA` for its whole life, and would write pivots whose name
#' claims a scale the values do not have.
#'
#' Idempotent: on a beta run, or on a resumed session where the keys are already
#' aligned, it changes nothing.
#'
#' @param ssEnv session environment, already carrying `beta`.
#' @return the session environment, with the SIGNAL keys realigned.
#' @keywords internal
#' @noRd
util_keys_signal_figure_refresh <- function(ssEnv) {

  figure <- io_signal_figure(ssEnv$beta)
  stale  <- c("MEAN", setdiff(c("BETA", "MVALUE"), figure))

  realign <- function(keys) {
    if (is.null(keys) || !all(c("MARKER", "FIGURE") %in% colnames(keys)))
      return(keys)
    target <- keys$MARKER == "SIGNAL" & keys$FIGURE %in% stale
    if (!any(target))
      return(keys)
    keys$FIGURE[target] <- figure
    if ("COMBINED" %in% colnames(keys)) {
      cols <- intersect(c("MARKER", "FIGURE", "AREA", "SUBAREA"), colnames(keys))
      keys$COMBINED[target] <- apply(keys[target, cols, drop = FALSE], 1, function(x)
        core_name_cleaning(gsub(" ", "", paste0(x[x != ""], collapse = "_"))))
    }
    keys
  }

  for (key_set in c("keys_markers_figures", "keys_areas_subareas_markers_figures",
                    "keys_markers_figures_default"))
    ssEnv[[key_set]] <- realign(ssEnv[[key_set]])

  if (!is.null(ssEnv$default$keys_markers_figures_default))
    ssEnv$default$keys_markers_figures_default <-
      realign(ssEnv$default$keys_markers_figures_default)

  ssEnv
}
