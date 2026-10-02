#' Read the study sample sheet, optionally joined with the per-sample statistics
#'
#' AI-255. The statistics used to live in a sibling CSV produced during the SEM
#' run, which meant the region classes had to be **foreseen six hours earlier**:
#' asking for one that had not been produced raised *"produce it with
#' semseeker(sample_stats_scopes = ...) and rerun the analysis"*. They are now
#' `SCOPE = SAMPLE` artefacts - pivots one row tall - and **this call is what
#' builds them**. Changing your mind costs one scan of the position pivot; the
#' written artefact is the cache, so the second call is a read.
#'
#' @param sql_sample_selection filter applied via [assoc_filter_sql()].
#' @param regions character vector of region classes to carry, named the way the
#'   columns name them: `"SAMPLE"` (no restriction - every position of the
#'   sample) or any `(AREA, SUBAREA)` pair of the run, e.g. `"GENE_TSS1500"`.
#'   Default `"SAMPLE"`, the historical behaviour.
#' @param with_sample_stats logical. When `TRUE` (default) the statistics are
#'   LEFT JOINed on `Sample_ID`, so consumers see `SAMPLE_<MARKER>_<FIGURE>_<AGG>`
#'   and `GENE_TSS1500_<MARKER>_<FIGURE>_<AGG>` as ordinary columns. Producers
#'   call it with `FALSE` to avoid feeding their own output back in.
#' @keywords internal
#' @noRd
sem_study_summary_get <- function(sql_sample_selection = "", regions = NULL,
                                  with_sample_stats = TRUE)
{
  ssEnv <- core_get_session_info()
  # io_file_path_build() uppercases via core_name_cleaning() - the on-disk file is
  # SAMPLE_SHEET_RESULT.csv, written by sem_study_summary_total(). Going through
  # io_file_path_build() here guarantees the same path is resolved on every OS
  # regardless of file-system case sensitivity.
  summary_file <- io_file_path_build( ssEnv$result_folderData, "sample_sheet_result","csv")
  if(!file.exists(summary_file))
    summary_file <- io_file_path_build( ssEnv$result_folderData, "1_SAMPLE_SHEET_ORIGINAL","csv")

  if(!file.exists(summary_file))
  {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " No sample sheet found in ", ssEnv$result_folderData)
    return()
  }

  study_summary <-   utils::read.csv2(summary_file)
  study_summary <- assoc_filter_sql(sql_sample_selection, study_summary)

  reference <- study_summary[study_summary$Sample_Group == "Reference",]
  no_reference <- study_summary[study_summary$Sample_Group != "Reference",]
  reference <- reference[!(reference$Sample_ID %in% no_reference$Sample_ID), ]

  study_summary <- rbind(no_reference, reference)

  if(nrow(reference)==0)
    study_summary <- no_reference

  if (isTRUE(with_sample_stats))
    study_summary <- .sem_study_summary_join_stats(study_summary, regions)

  return(study_summary)
}

#' Compose the per-sample statistics from the SAMPLE-scope artefacts
#'
#' AI-255. Each artefact is one row wide as the study; transposing it gives one
#' column per `(region, marker, figure, aggregation)`, which is exactly what a
#' model wants. Missing artefacts are built by [io_read_pivot()] on the way in.
#'
#' Silent no-op when nothing can be produced (e.g. a session that only ran the
#' exploratory step): the sample sheet is perfectly usable on its own, and
#' callers that need the burden check for their columns.
#'
#' @keywords internal
#' @noRd
.sem_study_summary_join_stats <- function(study_summary, regions = NULL) {

  if (is.null(study_summary) || !("Sample_ID" %in% colnames(study_summary)))
    return(study_summary)

  ssEnv <- core_get_session_info()
  keys  <- ssEnv$keys_markers_figures
  if (is.null(keys) || nrow(keys) == 0)
    return(study_summary)

  regions <- .sem_regions_resolve(regions)
  if (length(regions) == 0)
    return(study_summary)

  stats <- NULL
  for (region in regions) {
    for (i in seq_len(nrow(keys))) {
      marker <- as.character(keys$MARKER[i])
      figure <- as.character(keys$FIGURE[i])
      aggs <- util_aggregations_allowed(marker, figure,
                                        discrete = isTRUE(keys$DISCRETE[i]),
                                        default  = TRUE,
                                        scope    = "SAMPLE",
                                        area     = region$area)
      for (agg in aggs) {
        col <- .sem_scope_column_get(marker, figure, agg, region)
        if (is.null(col))
          next
        stats <- if (is.null(stats)) col else
          merge(stats, col, by = "Sample_ID", all = TRUE)
      }
    }
  }

  if (is.null(stats) || nrow(stats) == 0)
    return(study_summary)

  # AI-083, moved here by AI-255 with the join. The sample sheet identifiers are
  # the reference; the artefacts name their columns after the pivot columns. If
  # the two sets do not intersect at all, the join would silently produce a table
  # of NAs that looks like "this study has no burden" rather than like the
  # identifier mismatch it is.
  ids_sheet <- unique(as.character(study_summary$Sample_ID))
  ids_stats <- unique(as.character(stats$Sample_ID))
  if (length(intersect(ids_sheet, ids_stats)) == 0)
    stop("no Sample_ID in common between the sample sheet and the computed ",
         "artefacts. Sample sheet ids: ",
         paste(utils::head(ids_sheet, 4), collapse = ", "),
         " | artefact columns: ",
         paste(utils::head(ids_stats, 4), collapse = ", "), ".", call. = FALSE)

  missing_ids <- setdiff(ids_sheet, ids_stats)
  if (length(missing_ids) > 0)
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " ", length(missing_ids), " sample(s) have no computed statistics: ",
              paste(utils::head(missing_ids, 10), collapse = ", "),
              if (length(missing_ids) > 10) ", ..." else "")

  # Never let the statistics shadow a column of the sample sheet.
  duplicated_cols <- setdiff(intersect(colnames(stats), colnames(study_summary)),
                             "Sample_ID")
  if (length(duplicated_cols) > 0)
    stats <- stats[, !(colnames(stats) %in% duplicated_cols), drop = FALSE]

  merge(study_summary, stats, by = "Sample_ID", all.x = TRUE)
}

#' One SAMPLE-scope artefact, transposed into a Sample_ID/value data.frame
#'
#' @keywords internal
#' @noRd
.sem_scope_column_get <- function(marker, figure, aggregation, region) {

  lazy <- tryCatch(
    io_read_pivot(marker, figure, area = region$area, subarea = region$subarea,
                  aggregation = aggregation, scope = "SAMPLE"),
    error = function(e) {
      core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
                " study summary: ", marker, "_", figure, " ", aggregation,
                " on ", region$name, " unavailable: ", conditionMessage(e))
      NULL
    })
  if (is.null(lazy))
    return(NULL)

  row <- as.data.frame(lazy$collect())
  if (nrow(row) == 0)
    return(NULL)

  sample_cols <- setdiff(colnames(row), "AREA")
  if (length(sample_cols) == 0)
    return(NULL)

  out <- data.frame(Sample_ID = sample_cols,
                    value = as.numeric(row[1, sample_cols]),
                    stringsAsFactors = FALSE)
  colnames(out)[2] <- io_feature_colname(region$name, marker, figure, aggregation)
  out
}

#' Resolve requested region classes against the registry of the run
#'
#' AI-255. `"SAMPLE"` means no restriction - every position of the sample, which
#' in the taxonomy is the region class `PROBE_WHOLE`. It stays spelled `SAMPLE`
#' in the column prefixes because that is what a reader of the table expects to
#' see, and because the depth-1 consumer looks its columns up by that name.
#'
#' A region that does not resolve is an error, not a silent drop: a run that
#' quietly ignores the argument looks identical to one that honoured it.
#'
#' @keywords internal
#' @noRd
.sem_regions_resolve <- function(regions = NULL) {

  ssEnv <- core_get_session_info()
  if (is.null(regions) || length(regions) == 0)
    regions <- "SAMPLE"
  regions <- unique(core_name_cleaning(as.character(regions)))

  registry <- ssEnv$keys_areas_subareas
  available <- if (is.null(registry) || !("COMBINED" %in% colnames(registry)))
    character(0) else core_name_cleaning(as.character(registry$COMBINED))

  out <- list()
  for (name in regions) {
    if (identical(name, "SAMPLE")) {
      out[[length(out) + 1L]] <- list(name = "SAMPLE", area = "PROBE",
                                      subarea = "WHOLE")
      next
    }
    idx <- which(available == name)
    if (length(idx) == 0)
      stop("regions: unknown region class '", name,
           "'. A region class is an (AREA, SUBAREA) pair of this run. ",
           "Available: ", paste(c("SAMPLE", available), collapse = ", "), ".",
           call. = FALSE)
    out[[length(out) + 1L]] <- list(name = name,
                                    area    = as.character(registry$AREA[idx[1]]),
                                    subarea = as.character(registry$SUBAREA[idx[1]]))
  }
  out
}
