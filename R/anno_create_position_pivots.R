#' Create per-position pivot tables from per-sample bed/bedgraph files
#'
#' For each (marker, figure) in \code{keys} (restricted to the base markers
#' \code{MUTATIONS}, \code{DELTAR}, \code{DELTAS}), this function ensures the
#' POSITION/WHOLE pivot parquet on disk is up to date with respect to the
#' samples in \code{population}.
#'
#' The new samples are merged from their bed files by
#' \code{\link{io_stream_merge_bed}}. When the pivot already exists with a
#' subset of the samples, the new columns join it in one full outer join on
#' \code{(CHR, START, END)}, so a position found only by the new samples or
#' only by the old ones keeps its coordinates. Everything runs in the data
#' engine and the pivot never enters R: on the largest studies it holds
#' millions of positions for thousands of samples.
#'
#' @param population data.frame with \code{Sample_ID} and \code{Sample_Group}
#'   columns identifying the samples to materialise.
#' @param keys data.frame with \code{MARKER} and \code{FIGURE} columns
#'   describing which (marker, figure) pivots to build.
#'
#' @return Invisible \code{NULL}. Side effect: parquet pivot files and their
#'   JSON sidecars are created or updated under
#'   \code{Data/Pivots/<MARKER>/<MARKER>_<FIGURE>_POSITION_WHOLE_HG19.parquet}.
#'
#' @keywords internal
#' @noRd
anno_create_position_pivots <- function(population, keys) {

  ssEnv <- core_get_session_info()

  # Restrict to base markers - derived markers (DELTAP/DELTARP/DELTAQ/DELTARQ)
  # are built by sem_deltaX_get() post-SEM, not here. LESIONS is treated as a
  # base marker too: its per-sample bed files are written by sem_lesions_get()
  # during sem_analyze_population, and we need a POSITION pivot for annotate_*
  # to derive the PROBE/GENE/... aggregations later.
  selection <- c("MUTATIONS", "DELTAR", "DELTAS", "LESIONS")
  keys <- keys[order(keys$MARKER), ]
  keys <- keys[which(keys$MARKER %in% selection), , drop = FALSE]
  if (nrow(keys) == 0L) return(invisible())

  pop_clean <- population[!is.na(population$Sample_Group), ]

  for (k in seq_len(nrow(keys))) {

    key <- keys[k, ]
    marker <- as.character(key$MARKER)
    figure <- as.character(key$FIGURE)
    if (is.na(marker) || is.na(figure)) next

    area    <- "POSITION"
    subarea <- "WHOLE"
    pivot_filename <- io_pivot_file_name_parquet(marker, figure,
                                                         area, subarea)

    # Figure out which samples are still missing from the existing pivot.
    connection <- io_engine_connection()
    pivot_source <- sprintf("read_parquet(%s)", DBI::dbQuoteString(connection, pivot_filename))
    existing_samples <- character(0)
    if (file.exists(pivot_filename))
      existing_samples <- DBI::dbGetQuery(connection, paste("DESCRIBE SELECT * FROM", pivot_source))$column_name
    missing_pop <- pop_clean[!(pop_clean$Sample_ID %in% existing_samples), ]
    if (nrow(missing_pop) == 0L) next

    # Build the on-disk bed file paths for the missing samples and keep only
    # those that actually exist (some samples have no signal for a given
    # marker/figure → no bedgraph was emitted).
    bed_paths <- vapply(seq_len(nrow(missing_pop)), function(i) {
      io_bed_file_name(missing_pop$Sample_ID[i],
                                missing_pop$Sample_Group[i],
                                marker, figure)
    }, character(1))
    bed_paths <- bed_paths[file.exists(bed_paths)]
    if (length(bed_paths) == 0L && !file.exists(pivot_filename)) next

    dir.create(dirname(pivot_filename), recursive = TRUE, showWarnings = FALSE)
    tmp_filename <- paste0(pivot_filename, ".tmp")

    if (length(bed_paths) == 0L) {

      # No new bed file, only samples that found nothing: the existing pivot
      # gets their zero columns below.
      source_query <- paste("SELECT * FROM", pivot_source)
      present_columns <- existing_samples

    } else {

      core_log_event("INFO: ", Sys.time(),
                " anno_create_position_pivots[", marker, "_", figure,
                "] stream-merging ", length(bed_paths), " bed file(s)")

      merged_filename <- tempfile("merged_pivot_", fileext = ".parquet")
      io_stream_merge_bed(bed_paths, marker, figure, merged_filename)
      merged_source <- sprintf("read_parquet(%s)", DBI::dbQuoteString(connection, merged_filename))
      merged_columns <- DBI::dbGetQuery(connection, paste("DESCRIBE SELECT * FROM", merged_source))$column_name

      if (file.exists(pivot_filename)) {
        # USING merges the three key columns into one, taking each from
        # whichever side has the position; the old samples come first, then
        # the new ones in the order of their files.
        source_query <- sprintf("SELECT * FROM %s AS existing FULL JOIN %s AS merged USING (CHR, START, \"END\")",
                                pivot_source, merged_source)
        present_columns <- c(existing_samples, setdiff(merged_columns, existing_samples))
      } else {
        source_query <- paste("SELECT * FROM", merged_source)
        present_columns <- merged_columns
      }
    }

    # Every sample of the population gets its column. A sample that found
    # nothing writes no bed file, and without this its column would be
    # missing, which a reader cannot tell apart from a lost sample.
    absent <- setdiff(unique(pop_clean$Sample_ID), present_columns)
    zero_columns <- if (length(absent) == 0L) ""
                    else paste0(", ", paste0("CAST(0 AS DOUBLE) AS ", DBI::dbQuoteIdentifier(connection, absent), collapse = ", "))
    DBI::dbExecute(connection, sprintf("COPY (SELECT *%s FROM (%s) ORDER BY CHR, START) TO %s (FORMAT parquet, COMPRESSION zstd)",
                                       zero_columns, source_query, DBI::dbQuoteString(connection, tmp_filename)))

    file.rename(tmp_filename, pivot_filename)
    if (length(bed_paths) > 0L) unlink(merged_filename)
    # Sidecar JSON is now materialised by core_ensure_sidecars() at the end of the
    # pipeline (single point of responsibility).
    gc(verbose = FALSE)
  }

  invisible()
}
