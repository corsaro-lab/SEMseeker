#' Column names of the per-sample features (internal)
#'
#' AI-223. These names are the contract between the `SCOPE = SAMPLE` artefacts
#' and the models fitted on them: [sem_study_summary_get()] composes a column per
#' artefact, `association_analysis()` looks its features up by name.
#'
#' AI-255 removed the file they used to live in (`SAMPLE_STATS_RESULT.csv`); the
#' names outlived it because they name a *feature*, not a file.
#'
#' **The rule stated here is no longer the only place the string is built.**
#' Since AI-255 the region classes of a run come from the registry
#' `ssEnv$keys_areas_subareas`, and `.sem_regions_resolve()` takes the scope
#' string from its `COMBINED` column rather than calling this function, which is
#' left with no caller in `R/` (the tests still exercise it). On a region class
#' the two agree, because `combine_not_empty()` in [util_keys_create()] joins
#' `AREA` and `SUBAREA` by the same rule and cleans the result the same way; the
#' `SAMPLE` case they do *not* share - `COMBINED` gives `""` for an empty pair
#' and `.sem_regions_resolve()` intercepts `"SAMPLE"` before it reaches the
#' registry. So the agreement holds by construction, not by contract, which is
#' the arrangement this helper was written to prevent. Tracked in AI-259.
#'
#' The name of a scope is the region class it covers - `GENE_TSS1500` - or
#' `SAMPLE` when there is no restriction at all. It used to be derived from
#' `depth`, with `depth <= 1` meaning SAMPLE, `2` meaning the area alone and `3`
#' the area with its subarea. That mapping is gone with depth itself: the class
#' is said by the pair, and a number that ordered pairs on a line could only lie
#' about the ones that are not comparable.
#'
#' @param area,subarea the region class; both `NULL` gives `"SAMPLE"`.
#' @return the scope name, uppercase via [core_name_cleaning()].
#' @keywords internal
#' @noRd
io_scope_name <- function(area = NULL, subarea = NULL) {

  if (is.null(area) || !nzchar(as.character(area)))
    return("SAMPLE")

  parts <- as.character(area)
  if (!is.null(subarea) && nzchar(as.character(subarea)))
    parts <- c(parts, as.character(subarea))

  core_name_cleaning(paste(parts, collapse = "_"))
}
