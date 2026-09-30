#' Materialise the per-sample statistics of the run
#'
#' AI-223, rewritten by AI-255. This used to produce `SAMPLE_STATS_RESULT.csv`,
#' a second kind of artefact with its own shape - samples down the rows,
#' features across the columns - and its own producer, mask and aggregation
#' code. That shape is what hid the fact that a burden over the whole sample and
#' a burden per gene are the same marker reduced over different extents.
#'
#' They are now `SCOPE = SAMPLE` artefacts: ordinary pivots, **one row tall**,
#' written by the single derivation path ([io_pivot_build()]). The readable
#' per-sample table is composed at read time by [sem_study_summary_get()], which
#' transposes them and joins them onto the sample sheet.
#'
#' This pass is therefore a **warm-up, not a requirement**: it leaves the
#' statistics of the default region class on disk at the end of a SEM run, so a
#' user who only wants the descriptive table does not pay for it later. Anything
#' it does not materialise is built on first read - which is what removed the
#' old *"produce it and rerun the analysis"*.
#'
#' @return invisibly, the character vector of written artefact paths.
#' @keywords internal
#' @noRd
sem_sample_stats_build <- function() {

  ssEnv <- core_get_session_info()
  keys  <- ssEnv$keys_markers_figures

  if (is.null(keys) || nrow(keys) == 0) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " No marker/figure keys; per-sample statistics not materialised.")
    return(invisible(character(0)))
  }

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " Per-sample statistics: materialising the SAMPLE-scope artefacts ",
            "of the whole sample.")

  written <- character(0)
  for (i in seq_len(nrow(keys))) {
    marker <- as.character(keys$MARKER[i])
    figure <- as.character(keys$FIGURE[i])
    aggs   <- util_aggregations_allowed(marker, figure,
                                        discrete = isTRUE(keys$DISCRETE[i]),
                                        default  = TRUE,
                                        scope    = "SAMPLE",
                                        area     = "PROBE")
    out <- tryCatch(
      io_pivot_build(marker, figure, scope = "SAMPLE", area = "PROBE",
                     subarea = "WHOLE", aggregations = aggs),
      error = function(e) {
        core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                  " Per-sample statistics for ", marker, "_", figure,
                  " could not be built: ", conditionMessage(e))
        character(0)
      })
    written <- c(written, out)
  }

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " Per-sample statistics: ", length(written), " artefact(s) written.")

  invisible(written)
}

#' Usable positions per sample, a property of the imputation (internal)
#'
#' AI-255. `N_PROBES` is **not** an aggregation of a marker and not a property
#' of a scope: it is how many positions of that sample survived the treatment of
#' missing values. It therefore belongs on the sample sheet, with the other
#' descriptive properties of the sample, and not in the taxonomy - which is why
#' [util_aggregate_values()] no longer answers to the name.
#'
#' Nothing is lost for the density: the `MEAN` of a binary marker *is* the
#' density, denominator included.
#'
#' @return a data.frame with `Sample_ID` and `N_PROBES`, or `NULL`.
#' @keywords internal
#' @noRd
.sem_sample_n_probes <- function() {

  pivot <- tryCatch(
    io_read_pivot("SIGNAL", io_signal_figure(), "POSITION", "WHOLE",
                  aggregation = "VALUE", scope = "INSTANCE", build = FALSE),
    error = function(e) NULL)
  if (is.null(pivot))
    return(NULL)

  cols <- setdiff(names(pivot), c("CHR", "START", "END", "PROBE", "AREA"))
  if (length(cols) == 0)
    return(NULL)

  # polars' R bindings take expressions through `...`, not as a list: the list
  # has to be splatted with do.call(). Same idiom as sem_deltaX_get().
  counts <- as.data.frame(
    do.call(pivot$select, lapply(cols, function(cl)
      polars::pl$col(cl)$is_not_null()$sum()$alias(cl)))$collect())

  data.frame(Sample_ID = cols,
             N_PROBES  = as.numeric(counts[1, cols]),
             stringsAsFactors = FALSE)
}
