#' Compose a signal-descriptor column name - superseded, no live caller
#'
#' `<SCOPE>_<STAT>`, e.g. `SAMPLE_MEDIAN` or `GENE_BODY_IQR`.
#'
#' **Superseded twice over, and kept only until it is removed.**
#' Making the aggregation explicit replaced the two rules of the per-sample
#' table - this one for descriptors and
#' [io_burden_colname()] for burdens - with the single four-segment compositor
#' [io_feature_colname()], where the statistic became the `AGGREGATION` segment:
#' `SAMPLE_MEDIAN` is `SAMPLE_SIGNAL_BETA_MEDIAN` under that rule. The taxonomy unification then
#' removed the artefact these names described, `SAMPLE_STATS_RESULT.csv`.
#' Nothing in `R/` calls this function.
#'
#' Read it as a record of a superseded contract, not as one in force: a caller
#' that composes a name this way builds a string no reader of this package looks
#' up.
#'
#' @param scope scope name: a region class (`<AREA>` or `<AREA>_<SUBAREA>`), or
#'   `"SAMPLE"`.
#' @param stat one of the statistics declared by [io_signal_stats()].
#' @keywords internal
#' @noRd
io_stat_colname <- function(scope, stat) {
  if (missing(scope) || is.null(scope) || !nzchar(scope))
    stop("io_stat_colname(): 'scope' is required (use io_scope_name()).")
  core_name_cleaning(paste(scope, stat, sep = "_"))
}
