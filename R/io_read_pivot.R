# ---------------------------------------------------------------------------
# io_read_pivot() - unified dispatcher for marker pivot access
#
# The single read path for marker pivots: every reader of a HYPER / HYPO /
# BOTH pivot goes through here. The writers that build pivots
# (anno_create_position_pivots, the population bulk step) scan the files they
# are building directly: they extend a file, they do not consume a pivot.
#
# Two-branch dispatch:
#   (1) parquet pivot already materialised on disk -> scan_parquet (fast path)
#   (2) per-sample .bed / .bedgraph files present  -> streaming merge on the fly
#       (k-way full join on CHR/START/END, columnar projection per sample,
#        polars Streaming engine sinks the result without materialising it
#        in R memory).
# Falls back to NULL when neither storage is available, so callers can decide
# whether to skip or error.
#
# Rationale: the materialised pivot is treated as an OPTIONAL cache, not as
# the primary storage. Deleting it forces a re-merge from the bed files.
# This also enables analyses that touch only a subset of samples or genomic
# ranges to benefit from polars predicate/projection pushdown - they never
# read the full pivot from disk.

#' Read or stream-build a marker position pivot
#'
#' Unified accessor for the per-marker position pivot. If the pivot is already
#' materialised as parquet (the cache produced by
#' \code{\link{anno_create_position_pivots}}), it is returned as a lazy
#' \code{polars} frame. Otherwise, if the underlying per-sample \code{.bed} /
#' \code{.bedgraph(.gz)} files exist, a streaming row-wise merge over those
#' files is constructed lazily and returned. If neither exists, returns
#' \code{NULL}.
#'
#' This function delegates parallelism and memory management to the polars
#' Rust runtime (Rayon + Streaming engine), so callers see a single
#' \code{LazyFrame} regardless of the backend.
#'
#' @param marker  Character. Marker name, e.g. \code{"MUTATIONS"},
#'   \code{"DELTAR"}, \code{"DELTAS"}.
#' @param figure  Character. Figure name, e.g. \code{"HYPER"}, \code{"HYPO"},
#'   \code{"MEAN"}.
#' @param area    Character. Area label, default \code{"POSITION"}.
#' @param subarea Character. Subarea label, default \code{"WHOLE"}.
#' @param cache   Logical. If \code{TRUE} (default), when CASE 2 is used the
#'   streaming merge result is persisted to the parquet pivot path via
#'   \code{sink_parquet()} (lazy streaming write, no in-memory materialisation).
#'   Subsequent calls will then hit CASE 1 (fast path). Set to \code{FALSE} for
#'   one-shot reads where no future access is expected, to avoid writing the
#'   pivot to disk.
#'
#' @return A lazy \code{polars} frame, or \code{NULL} if no storage is found.
#'
#' @keywords internal
#' @noRd
io_read_pivot <- function(marker, figure, area = "POSITION", subarea = "WHOLE",
                       cache = TRUE, aggregation = NULL, scope = "INSTANCE",
                       build = TRUE) {

  ssEnv <- core_get_session_info()
  scope <- io_scope_validate(scope)

  # --- Branch 1: cached parquet pivot ----------------------------------------
  pivot_filename <- io_pivot_file_name_parquet(marker, figure, area, subarea,
                                               aggregation = aggregation,
                                               scope = scope)
  if (file.exists(pivot_filename)) {
    core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
              " io_read_pivot[", marker, "_", figure,
              "] hit cached parquet: ", basename(pivot_filename))
    return(polars::pl$scan_parquet(pivot_filename))
  }

  # --- Branch 1b: derived artefact, built on demand --------------------------
  # AI-255. Anything that is not the base position pivot is derived, and it is
  # derived HERE rather than having to be foreseen at SEM time. This is what
  # removes the old "produce it with semseeker(...) and rerun the analysis":
  # changing your mind now costs one scan instead of a whole run. The written
  # file is the cache - the second call takes branch 1.
  is_base <- identical(scope, "INSTANCE") && io_area_is_single_position(area) &&
             identical(toupper(as.character(subarea)), "WHOLE")
  if (isTRUE(build) && !is_base) {
    agg <- .io_aggregation_resolve(aggregation, scope, area)

    # AI-255: build EVERY aggregation this key admits, not just the one asked
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
              " io_read_pivot[", marker, "_", figure,
              "] building ", basename(pivot_filename), " on demand, with ",
              length(aggregations), " aggregation(s) in one pass")
    built <- io_pivot_build(marker, figure, scope = scope, area = area,
                            subarea = subarea, aggregations = aggregations)
    if (length(built) > 0 && file.exists(pivot_filename))
      return(polars::pl$scan_parquet(pivot_filename))
    return(NULL)
  }

  # --- Branch 2: per-sample bed/bedgraph files -------------------------------
  # NB: this branch builds a wide pivot with CHR/START/END as the key columns
  # (a POSITION-shape table). It is only valid when area == "POSITION";
  # PROBE/GENE/ISLAND/... pivots have an AREA key column and are produced by
  # anno_annotate_position_pivots() from the POSITION pivot. Trying to fall back to
  # io_stream_merge_bed for non-POSITION areas would write a POSITION-shape file
  # under a PROBE-shape filename and silently corrupt downstream reads.
  if (!identical(area, "POSITION")) {
    return(NULL)
  }
  bed_files <- io_list_bed_files_for_marker_figure(marker, figure)
  if (length(bed_files) > 0) {
    core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
              " io_read_pivot[", marker, "_", figure,
              "] streaming merge over ", length(bed_files), " bed/bedgraph files",
              if (cache) " (will persist via sink_parquet)" else " (one-shot, no cache)")
    lazy <- io_stream_merge_bed(bed_files, marker, figure)

    if (isTRUE(cache)) {
      dir.create(dirname(pivot_filename), recursive = TRUE, showWarnings = FALSE)
      lazy$sink_parquet(pivot_filename)
      # Sidecar JSON is materialised by core_ensure_sidecars() at pipeline end.
      return(polars::pl$scan_parquet(pivot_filename))
    }
    return(lazy)
  }

  # --- Branch 3: nothing to read ---------------------------------------------
  core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
            " io_read_pivot[", marker, "_", figure,
            "] no pivot parquet and no bed/bedgraph files found; returning NULL")
  return(NULL)
}

#' List per-sample bed/bedgraph files for a given marker × figure
#'
#' Discovers all \code{.bed}, \code{.bedgraph}, \code{.bedgraph.gz} files for
#' the requested marker/figure across the \code{Data/<sample_group>/<marker>_<figure>/}
#' directory tree produced by \code{\link{sem_analyze_population}}.
#'
#' @keywords internal
#' @noRd
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
