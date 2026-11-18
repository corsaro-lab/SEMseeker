# ---------------------------------------------------------------------------
# io_read_pivot() - unified dispatcher for marker pivot access
#
# The single read path for marker pivots: every reader of a HYPER / HYPO /
# BOTH pivot goes through here. The writers that build pivots
# (anno_create_position_pivots, the population bulk step) scan the files they
# are building directly: they extend a file, they do not consume a pivot.
#
# Producing the file and reading it are separate steps. io_pivot_ensure()
# makes the parquet pivot available on disk and returns its path:
#   1. parquet pivot already materialised on disk -> its path (fast path)
#   2. derived artefact (aggregated area, SAMPLE scope) -> built on demand
#   3. per-sample .bed / .bedgraph files present  -> merged and written
# and returns NULL when none applies, so callers can decide whether to skip
# or error. io_read_pivot() scans the path lazily with the data engine.
#
# The parquet file is the contract between writers and readers, whatever
# engine each uses: while SEMseeker moves from polars to DuckDB, the callers
# not yet moved read the same file through io_read_pivot_polars().
#
# Rationale: the materialised pivot is treated as an OPTIONAL cache, not as
# the primary storage. Deleting it forces a re-merge from the bed files.
# Analyses that touch only a subset of samples or genomic ranges benefit
# from predicate and projection pushdown: they never read the full pivot.

#' Read a marker pivot as a lazy table
#'
#' Unified accessor for marker pivots: the file is made available by
#' \code{io_pivot_ensure()} (cached parquet, built on demand, or merged from
#' the per-sample bedGraph files) and scanned lazily with the data engine.
#' Nothing is read until the caller collects or sinks the query.
#'
#' @param marker  Character. Marker name, e.g. \code{"MUTATIONS"},
#'   \code{"DELTAR"}, \code{"DELTAS"}.
#' @param figure  Character. Figure name, e.g. \code{"HYPER"}, \code{"HYPO"}.
#' @param area    Character. Area label, default \code{"POSITION"}.
#' @param subarea Character. Subarea label, default \code{"WHOLE"}.
#' @param aggregation Character or NULL. Aggregation of a derived artefact;
#'   NULL takes the default of the key.
#' @param scope   Character. \code{"INSTANCE"} (default) or \code{"SAMPLE"}.
#' @param build   Logical. When TRUE (default) a missing derived artefact is
#'   built from the position pivot; when FALSE only an existing file is read.
#'
#' @return A \code{tbl_lazy} on the data engine connection, or \code{NULL}
#'   when the pivot is neither on disk nor buildable.
#'
#' @keywords internal
#' @noRd
io_read_pivot <- function(marker, figure, area = "POSITION", subarea = "WHOLE",
                          aggregation = NULL, scope = "INSTANCE", build = TRUE) {
  pivot_filename <- io_pivot_ensure(marker, figure, area = area, subarea = subarea,
                                    aggregation = aggregation, scope = scope, build = build)
  if (is.null(pivot_filename))
    return(NULL)
  io_scan_parquet(pivot_filename)
}

#' Read a marker pivot as a lazy polars frame, for the callers not yet moved
#' to the data engine.
#'
#' Transitional: the population and association engines still build polars
#' queries on the pivot. They read the same file io_read_pivot() reads, so the
#' bridge costs no memory, and each caller drops it when it moves to the data
#' engine. Removed with polars.
#'
#' @inheritParams io_read_pivot
#' @return A lazy \code{polars} frame, or \code{NULL} as io_read_pivot().
#' @keywords internal
#' @noRd
io_read_pivot_polars <- function(marker, figure, area = "POSITION", subarea = "WHOLE",
                                 aggregation = NULL, scope = "INSTANCE", build = TRUE) {
  pivot_filename <- io_pivot_ensure(marker, figure, area = area, subarea = subarea,
                                    aggregation = aggregation, scope = scope, build = build)
  if (is.null(pivot_filename))
    return(NULL)
  polars::pl$scan_parquet(pivot_filename)
}

#' Make a marker pivot available on disk and return its path.
#'
#' Separates producing the pivot file from reading it. The parquet file is the
#' contract between the code that writes pivots and the code that reads them,
#' whatever engine each uses, so producing it once here lets every reader
#' choose how to scan it. Three cases, in order:
#'   1. the parquet pivot is already on disk: its path;
#'   2. a derived artefact (any key but the base position pivot) is built from
#'      the position pivot, with every aggregation the key admits in one pass;
#'   3. the base POSITION pivot is merged from the per-sample bedGraph files
#'      and written to its path.
#' The written file is the cache: the next call takes case 1.
#'
#' @inheritParams io_read_pivot
#' @return the absolute path of the parquet pivot, or \code{NULL} when it is
#'   neither on disk nor buildable.
#' @keywords internal
#' @noRd
io_pivot_ensure <- function(marker, figure, area = "POSITION", subarea = "WHOLE",
                            aggregation = NULL, scope = "INSTANCE", build = TRUE) {

  scope <- io_scope_validate(scope)

  # --- Case 1: cached parquet pivot ------------------------------------------
  pivot_filename <- io_pivot_file_name_parquet(marker, figure, area, subarea,
                                               aggregation = aggregation,
                                               scope = scope)
  if (file.exists(pivot_filename)) {
    core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
              " io_pivot_ensure[", marker, "_", figure,
              "] hit cached parquet: ", basename(pivot_filename))
    return(pivot_filename)
  }

  # --- Case 2: derived artefact, built on demand -----------------------------
  # Anything that is not the base position pivot is derived, and it is
  # derived HERE rather than having to be foreseen at SEM time. This is what
  # removes the old "produce it with semseeker(...) and rerun the analysis":
  # changing your mind now costs one scan instead of a whole run.
  is_base <- identical(scope, "INSTANCE") && io_area_is_single_position(area) &&
             identical(toupper(as.character(subarea)), "WHOLE")
  if (isTRUE(build) && !is_base) {
    agg <- .io_aggregation_resolve(aggregation, scope, area)

    # build EVERY aggregation this key admits, not just the one asked
    # for. The scan of the position pivot is the expensive part and it is shared:
    # one group_by emits SUM, MEAN, MEDIAN, VARIANCE and IQR together. Building
    # them one at a time would turn three requests on the same area into three
    # scans - which is the cost the taxonomy is supposed to avoid, and which the
    # producer this replaced already avoided. The extra artefacts are small and
    # are the cache for the next request.
    aggregations <- tryCatch(
      unique(c(agg, util_aggregations_allowed(marker, figure, default = FALSE,
                                              scope = scope, area = area))),
      error = function(e) agg)
    # The two modes need the whole distribution per group and are handled R-side;
    # do not drag them into an opportunistic build.
    aggregations <- setdiff(aggregations, setdiff(c("MODELOW", "MODEHIGH"), agg))

    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " io_pivot_ensure[", marker, "_", figure,
              "] building ", basename(pivot_filename), " on demand, with ",
              length(aggregations), " aggregation(s) in one pass")
    built <- io_pivot_build(marker, figure, scope = scope, area = area,
                            subarea = subarea, aggregations = aggregations)
    if (length(built) > 0 && file.exists(pivot_filename))
      return(pivot_filename)
    return(NULL)
  }

  # --- Case 3: per-sample bed/bedgraph files ---------------------------------
  # NB: this case builds a wide pivot with CHR/START/END as the key columns
  # (a POSITION-shape table). It is only valid when area == "POSITION";
  # PROBE/GENE/ISLAND/... pivots have an AREA key column and are produced by
  # anno_annotate_position_pivots() from the POSITION pivot. Falling back to
  # io_stream_merge_bed for non-POSITION areas would write a POSITION-shape file
  # under a PROBE-shape filename and silently corrupt downstream reads.
  if (!identical(area, "POSITION")) {
    return(NULL)
  }
  bed_files <- io_list_bed_files_for_marker_figure(marker, figure)
  if (length(bed_files) > 0) {
    core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
              " io_pivot_ensure[", marker, "_", figure,
              "] merging ", length(bed_files), " bed/bedgraph files into ",
              basename(pivot_filename))
    dir.create(dirname(pivot_filename), recursive = TRUE, showWarnings = FALSE)
    io_stream_merge_bed(bed_files, marker, figure)$sink_parquet(pivot_filename)
    # Sidecar JSON is materialised by core_ensure_sidecars() at pipeline end.
    return(pivot_filename)
  }

  # --- Nothing to read -------------------------------------------------------
  core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
            " io_pivot_ensure[", marker, "_", figure,
            "] no pivot parquet and no bed/bedgraph files found; returning NULL")
  return(NULL)
}

io_list_bed_files_for_marker_figure <- function(marker, figure) {
  ssEnv <- core_get_session_info()
  data_root <- ssEnv$result_folderData
  if (is.null(data_root) || !dir.exists(data_root)) return(character(0))

  subdir_name <- paste0(marker, "_", figure)
  # Look under each Sample_Group folder (Case, Control, Reference, ...)
  candidates <- list.files(
    data_root,
    pattern    = "\\.(bed|bedgraph)(\\.gz)?$",
    recursive  = TRUE,
    full.names = TRUE
  )
  # Keep only files under .../<sample_group>/<marker>_<figure>/...
  pattern <- paste0("/", subdir_name, "/")
  candidates[grepl(pattern, candidates, fixed = TRUE)]
}
