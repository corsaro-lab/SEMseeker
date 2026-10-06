#' Compose a per-marker burden column name - superseded, no live caller
#'
#' `<SCOPE>_<MARKER>_<FIGURE>`, e.g. `SAMPLE_MUTATIONS_HYPER`.
#'
#' **Superseded twice over, and kept only until it is removed.**
#' Making the aggregation explicit replaced the two rules of the per-sample
#' table - this one for burdens and
#' [io_stat_colname()] for descriptors - with the single four-segment compositor
#' [io_feature_colname()], which makes the operator explicit instead of implying
#' it from the marker: `SAMPLE_MUTATIONS_HYPER` is `SAMPLE_MUTATIONS_HYPER_SUM`
#' under that rule. The taxonomy unification then removed the artefact these names described,
#' `SAMPLE_STATS_RESULT.csv`. Nothing in `R/` calls this function.
#'
#' The missing segment is the reason it cannot come back as it is: a name
#' without `AGGREGATION` cannot tell the sum of a marker from its mean.
#'
#' @param scope scope name: a region class (`<AREA>` or `<AREA>_<SUBAREA>`), or
#'   `"SAMPLE"`.
#' @param marker marker name.
#' @param figure figure name (`HYPER` / `HYPO`).
#' @keywords internal
#' @noRd
io_burden_colname <- function(scope, marker, figure) {
  if (missing(scope) || is.null(scope) || !nzchar(scope))
    stop("io_burden_colname(): 'scope' is required (use io_scope_name()).")
  core_name_cleaning(paste(scope, marker, figure, sep = "_"))
}
