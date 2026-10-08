#' Create per-position pivot tables from per-sample bed/bedgraph files
#'
#' For each (marker, figure) in \code{keys} (restricted to the base markers
#' \code{MUTATIONS}, \code{DELTAR}, \code{DELTAS}), this function ensures the
#' POSITION/WHOLE pivot parquet on disk is up to date with respect to the
#' samples in \code{population}.
#'
#' Implementation note (2026-06-01, unified pivot reader):
#'   The previous R-side \code{for (s in 1:N)} loop with per-sample
#'   \code{readr::read_tsv} + \code{$join(how = "full")} paid the R→Rust FFI
#'   cost N times and ran single-threaded. It has been replaced by a single
#'   call to \code{\link{io_stream_merge_bed}}, which scans all missing bed files
#'   in parallel inside the polars Rust runtime and pivots long→wide in one
#'   shot. When the pivot already exists with a subset of samples, the new
#'   wide block is merged into it via a single full-outer join on
#'   \code{(CHR, START, END)} with coalesce on the key columns.
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
    existing_samples <- character(0)
    if (file.exists(pivot_filename)) {
      existing_samples <- colnames(
        polars::pl$scan_parquet(pivot_filename, n_rows = 1)$collect()
      )
    }
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
      out <- polars::pl$scan_parquet(pivot_filename)

    } else {

      core_log_event("INFO: ", Sys.time(),
                " anno_create_position_pivots[", marker, "_", figure,
                "] stream-merging ", length(bed_paths), " bed file(s)")

      # One-shot lazy pivot built by polars Rust runtime (no R-side loop).
      new_lazy <- io_stream_merge_bed(bed_paths, marker, figure)

      out <- if (file.exists(pivot_filename)) {
        # Merge new sample columns into the existing pivot via one full-outer
        # join on the genomic key. Coalesce CHR/START/END to avoid the
        # _right duplicate columns polars emits on a full join.
        polars::pl$scan_parquet(pivot_filename)$join(
          new_lazy,
          on  = c("CHR", "START", "END"),
          how = "full"
        )$with_columns(
          polars::pl$when(polars::pl$col("CHR")$is_not_null())$
            then(polars::pl$col("CHR"))$
            otherwise(polars::pl$col("CHR_right"))$alias("CHR"),
          polars::pl$when(polars::pl$col("START")$is_not_null())$
            then(polars::pl$col("START"))$
            otherwise(polars::pl$col("START_right"))$alias("START"),
          polars::pl$when(polars::pl$col("END")$is_not_null())$
            then(polars::pl$col("END"))$
            otherwise(polars::pl$col("END_right"))$alias("END")
        )$drop(c("CHR_right", "START_right", "END_right"))
      } else {
        new_lazy
      }
    }

    # Every sample of the population gets its column. A sample that found
    # nothing writes no bed file, and without this its column would be
    # missing, which a reader cannot tell apart from a lost sample.
    absent <- setdiff(unique(pop_clean$Sample_ID), names(out))
    if (length(absent) > 0L)
      out <- do.call(out$with_columns,
                     lapply(absent, function(id) polars::pl$lit(0)$alias(id)))

    out$collect()$
      sort(c("CHR", "START"), descending = FALSE)$
      write_parquet(tmp_filename)

    file.rename(tmp_filename, pivot_filename)
    # Sidecar JSON is now materialised by core_ensure_sidecars() at the end of the
    # pipeline (single point of responsibility).
    gc(verbose = FALSE)
  }

  invisible()
}
