#' Compose the column name of a computed feature (internal)
#'
#' **The** single compositor of column names for the study summary
#' table, used by the side that builds a column and by the side that looks it
#' up. If the two built the string independently the contract would break at the
#' first divergence, so nothing else in the package is allowed to paste these
#' names together.
#'
#' \preformatted{
#' <SCOPE>_<MARKER>_<FIGURE>_<AGGREGATION>
#'
#' SCOPE        SAMPLE (no restriction) | <AREA> | <AREA>_<SUBAREA>
#' MARKER       SIGNAL | MUTATIONS | LESIONS | DELTA*
#' FIGURE       HYPER / HYPO for the instability markers,
#'              BETA / MVALUE (the scale) for SIGNAL
#' AGGREGATION  SUM | MEAN | MEDIAN | VARIANCE | IQR | MODELOW | MODEHIGH
#' }
#'
#' **This is not the artefact key, and its `SCOPE` is not the `SCOPE`
#' coordinate.** The key of an artefact has six coordinates and lives in
#' [io_pivot_file_name()]: `MARKER`, `FIGURE`, `SCOPE` (`SAMPLE` | `INSTANCE`,
#' the extent one number is valid over), `AREA`, `SUBAREA`, `AGGREGATION`
#'. A *column* of the summary table is a narrower thing: its rows are
#' samples, so the artefact behind it is fixed at `SCOPE = SAMPLE` and there is
#' nothing left to say about the extent. What the reader needs to see in the
#' prefix is the region class the number covers - `GENE_TSS1500` - so `AREA` and
#' `SUBAREA` collapse into the first segment, and `SAMPLE` names the absence of
#' any restriction. One taxonomy, two projections of it; the four segments here
#' are what survives of the six once the row is a sample.
#'
#' The aggregation used to be implicit - one operator per marker, so
#' `MARKER`/`FIGURE` identified the value. With several aggregations over the
#' same scope it no longer does, hence the fourth segment.
#'
#' Axes left `NULL` are dropped, so a quantity that is not the aggregation of a
#' marker can still be named through the same compositor:
#' `io_feature_colname("SAMPLE", aggregation = "N_PROBES")` gives
#' `SAMPLE_N_PROBES`. That mechanism outlived its original occasion - since
#' Since the taxonomy unification **no artefact carries `N_PROBES`** (see [io_signal_stats()]): the
#' count is a property of the imputation and travels on `SAMPLE_SHEET_RESULT`.
#'
#' Live callers in `R/`: one, the composer inside [sem_study_summary_get()].
#'
#' @param scope scope name: a region class (`<AREA>` or `<AREA>_<SUBAREA>`), or
#'   `"SAMPLE"`. Required.
#' @param marker,figure,aggregation the remaining axes; `NULL` to omit.
#' @return the column name, uppercase via [core_name_cleaning()].
#' @keywords internal
#' @noRd
io_feature_colname <- function(scope, marker = NULL, figure = NULL,
                               aggregation = NULL) {
  if (missing(scope) || is.null(scope) || !all(nzchar(scope)))
    stop("io_feature_colname(): 'scope' is required (use io_scope_name()).")

  parts <- list(scope, marker, figure, aggregation)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  parts <- lapply(parts, as.character)

  core_name_cleaning(do.call(paste, c(parts, list(sep = "_"))))
}
