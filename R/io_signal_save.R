#' Save the signal matrix as the SIGNAL pivots, by probe and by position
#'
#' Every later step reads the signal from these two parquet files and not
#' from the matrix in memory: the PROBE pivot keeps the probe ids, the
#' POSITION pivot the genomic coordinates the markers are computed on. The
#' POSITION pivot is sorted by chromosome in natural order (1, 2, ..., 22,
#' X, Y, M) and then by START and END, so each chromosome is one contiguous
#' block: the per-chromosome model fit relies on it. Nothing is written when
#' the POSITION pivot already exists.
#'
#' @param signal_data data.frame of signal values, probes in the rows (probe
#'   ids as row names) and samples in the columns
#' @param sample_sheet data.frame with a Sample_ID column; only its samples
#'   are saved, in its order
#' @param batch_id identifier of the batch; kept for the callers, not used
#' @param probe_features data.frame with PROBE, CHR, START and END, the
#'   coordinates of the probes; NULL takes the attribute "probe_features" of
#'   signal_data and, failing that, the annotation of the technology
#' @return NULL, invisibly; the side effect is the two parquet files
#' @keywords internal
#' @noRd
io_signal_save <- function(signal_data, sample_sheet, batch_id,
                        probe_features = NULL)
{
  ssEnv <- core_get_session_info()
  # Resolve probe_features in priority order: explicit arg > attribute
  # attached by sem_prepare_batch_signal() > legacy anno_probe_features_get() refetch.
  # In the normal pipeline sem_analyze_batch() routes signal_data through
  # sem_prepare_batch_signal() so attr(., "probe_features") is set and the
  # third branch is never taken.
  if (is.null(probe_features))
    probe_features <- attr(signal_data, "probe_features")
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Saving signal data.")
  pivot_file_name_pos <- io_pivot_file_name_parquet("SIGNAL", io_signal_figure(), "POSITION", "WHOLE")
  if (file.exists(pivot_file_name_pos)) {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Signal data already saved.")
    return()
  }

  # normalise both sides and fail with the offending identifiers rather
  # than with R's "undefined columns selected" on the subset below.
  .normalized <- core_normalize_sample_ids(sample_sheet, signal_data,
                                           require_all_ids = TRUE)
  sample_sheet <- .normalized$sample_sheet
  signal_data  <- .normalized$signal_data
  rm(.normalized)

  signal_data <- signal_data[, unique(sample_sheet$Sample_ID), drop = FALSE]

  # ------------------------------------------------------------------
  # WGBS / LONGREAD path - coordinates are encoded in synthetic probe IDs.
  # No Bioconductor annotation join is needed.
  # ------------------------------------------------------------------
  if (ssEnv$tech %in% c("WGBS", "LONGREAD")) {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              "Saving signal data for ", ssEnv$tech,
              " (extracting coordinates from synthetic probe IDs).")

    # Save probe-indexed parquet (synthetic IDs as AREA column)
    signal_probe           <- signal_data
    signal_probe$AREA      <- rownames(signal_data)
    signal_probe           <- signal_probe[, c(ncol(signal_probe),
                                               seq_len(ncol(signal_probe) - 1))]
    pivot_file_name_probe  <- io_pivot_file_name_parquet("SIGNAL", io_signal_figure(), "PROBE", "WHOLE")
    io_write_frame_parquet(signal_probe, pivot_file_name_probe)
    rm(signal_probe)

    # Build position-indexed parquet directly from synthetic probe IDs
    coords       <- io_probe_id_to_coord(rownames(signal_data))
    signal_pos   <- data.frame(
      CHR   = coords$CHR,
      START = coords$START,
      END   = coords$END,
      signal_data,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    # Sort by chromosome (natural order) then position
    chr_order    <- order(
      suppressWarnings(as.integer(signal_pos$CHR)),   # numeric chrs first
      nchar(signal_pos$CHR),                          # X/Y/M after
      signal_pos$CHR,
      signal_pos$START
    )
    signal_pos   <- signal_pos[chr_order, ]
    io_write_frame_parquet(signal_pos, pivot_file_name_pos)
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Saved signal data (", ssEnv$tech, ").")
    gc()
    return()
  }

  # ------------------------------------------------------------------
  # Illumina path - join with Bioconductor annotation to get CHR/START/END
  # ------------------------------------------------------------------
  signal_data$AREA <- rownames(signal_data)
  signal_data      <- signal_data[, c(ncol(signal_data), seq_len(ncol(signal_data) - 1))]

  pivot_file_name_probe <- io_pivot_file_name_parquet("SIGNAL", io_signal_figure(), "PROBE", "WHOLE")
  io_write_frame_parquet(signal_data, pivot_file_name_probe)
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Signal data saved with probe.")

  rm(signal_data)
  gc()
  core_log_event("DEBUG_MEM_SS: ", format(Sys.time(), "%a %b %d %X %Y"), " post-probe-write+gc  mem_MB=", round(sum(gc()[, "(Mb)"]), 1))

  pf <- if (!is.null(probe_features)) probe_features
        else anno_probe_features_get("PROBE")  # legacy fallback (test path)
  # Only the coordinates and the probe id: probe_features also carries the
  # annotation columns (GENE_*, ISLAND_*, ...), which must not reach the
  # POSITION pivot, whose sample columns come right after CHR, START, END.
  pf <- data.frame(CHR   = as.character(pf$CHR),
                   START = pf$START,
                   END   = pf$END,
                   PROBE = as.character(pf$PROBE),
                   stringsAsFactors = FALSE)
  # Natural chromosome order; chromosomes outside 1-22, X, Y, M keep the
  # order in which they first appear, after the others.
  chr_order_key <- function(chromosomes) {
    # chromosomes: CHR values, with or without a leading "chr"; returns the
    # rank of each: its number, 23 for X, 24 for Y, 25 for M or MT, 99 else.
    bare <- sub("^chr", "", chromosomes)
    number <- suppressWarnings(as.integer(bare))
    ifelse(!is.na(number), number,
           ifelse(bare == "X", 23L,
                  ifelse(bare == "Y", 24L,
                         ifelse(bare %in% c("M", "MT"), 25L, 99L))))
  }
  chromosomes <- unique(pf$CHR)
  chromosomes <- chromosomes[order(chr_order_key(chromosomes))]
  pf$CHR_RANK <- match(pf$CHR, chromosomes)

  # io_signal_position_write() sorts in one query and falls back on one
  # chromosome at a time when the sort does not fit in the memory limit.
  connection <- io_engine_connection()
  coordinates_view <- paste0("signal_coordinates_", gsub("[^A-Za-z0-9]", "", basename(tempfile(""))))
  duckdb::duckdb_register(connection, coordinates_view, pf)
  on.exit(duckdb::duckdb_unregister(connection, coordinates_view), add = TRUE)
  io_signal_position_write(connection, coordinates_view, length(chromosomes), pivot_file_name_probe, pivot_file_name_pos)
  rm(pf)

  gc()
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Saved signal data.")
  core_log_event("DEBUG_MEM_SS: ", format(Sys.time(), "%a %b %d %X %Y"), " post-position-write  mem_MB=", round(sum(gc()[, "(Mb)"]), 1))
}

#' Write the POSITION pivot of the signal, falling back on one chromosome at a time
#'
#' One query joins the coordinates to the PROBE pivot, sorts and writes: the
#' fast way. Under a memory limit the sort of the whole pivot can run out of
#' memory (it did under 8 GB on 400 000 probes x 1 000 samples, see the
#' engine migration notes), so an out-of-memory error switches to one sort
#' per chromosome, about three times slower. Any other error stops as
#' before. A partial pivot left by the failed query is removed first, or the
#' next run would find it and skip the save.
#'
#' @param connection the engine connection, from io_engine_connection()
#' @param coordinates_view name of a registered view with CHR, START, END,
#'   PROBE and CHR_RANK, the rank of the chromosome in natural order, from 1
#' @param chromosome_count number of distinct CHR_RANK values in the view
#' @param probe_path path of the SIGNAL PROBE pivot, keyed on AREA
#' @param position_path path of the POSITION pivot to write
#' @return NULL, invisibly; the side effect is the file at position_path
#' @keywords internal
#' @noRd
io_signal_position_write <- function(connection, coordinates_view, chromosome_count, probe_path, position_path) {
  tryCatch(io_signal_position_single(connection, coordinates_view, probe_path, position_path), error = function(e) {
    if (!grepl("Out of Memory", conditionMessage(e), fixed = TRUE)) stop(e)
    unlink(position_path)
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                   " the POSITION pivot did not fit in the engine's memory limit in one sort; writing it one chromosome at a time.")
    io_signal_position_per_chromosome(connection, coordinates_view, chromosome_count, probe_path, position_path)
  })
  invisible(NULL)
}

#' Write the POSITION pivot of the signal with a single join, sort and copy
#'
#' The fast way, used whenever the sort fits in the engine's memory limit.
#'
#' @inheritParams io_signal_position_write
#' @return NULL, invisibly
#' @keywords internal
#' @noRd
io_signal_position_single <- function(connection, coordinates_view, probe_path, position_path) {
  DBI::dbExecute(connection, sprintf("
    COPY (
      SELECT coordinates.CHR, coordinates.START, coordinates.\"END\", signal.* EXCLUDE (AREA)
      FROM %s AS coordinates JOIN read_parquet(%s) AS signal ON coordinates.PROBE = signal.AREA
      ORDER BY coordinates.CHR_RANK, coordinates.START, coordinates.\"END\"
    ) TO %s (FORMAT parquet, COMPRESSION zstd)",
    coordinates_view, DBI::dbQuoteString(connection, probe_path), DBI::dbQuoteString(connection, position_path)))
  invisible(NULL)
}

#' Write the POSITION pivot of the signal one chromosome at a time
#'
#' The fallback when one sort does not fit: each chromosome is joined,
#' sorted on START and END and written to its own file, then the pieces are
#' copied in chromosome order into the pivot, a copy and not a sort. The
#' pieces are read in the order of the list, and preserve_insertion_order
#' (on by default) keeps it in the copy. Each piece rereads the whole PROBE
#' pivot, hence the cost.
#'
#' @inheritParams io_signal_position_write
#' @return NULL, invisibly
#' @keywords internal
#' @noRd
io_signal_position_per_chromosome <- function(connection, coordinates_view, chromosome_count, probe_path, position_path) {
  piece_folder <- tempfile("signal_chr_pieces_")
  dir.create(piece_folder, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(piece_folder, recursive = TRUE), add = TRUE)
  piece_paths <- file.path(piece_folder, sprintf("%03d.parquet", seq_len(chromosome_count)))
  for (chromosome_rank in seq_len(chromosome_count)) {
    DBI::dbExecute(connection, sprintf("
      COPY (
        SELECT coordinates.CHR, coordinates.START, coordinates.\"END\", signal.* EXCLUDE (AREA)
        FROM (SELECT * FROM %s WHERE CHR_RANK = %d) AS coordinates
        JOIN read_parquet(%s) AS signal ON coordinates.PROBE = signal.AREA
        ORDER BY coordinates.START, coordinates.\"END\"
      ) TO %s (FORMAT parquet, COMPRESSION zstd)",
      coordinates_view, chromosome_rank, DBI::dbQuoteString(connection, probe_path),
      DBI::dbQuoteString(connection, piece_paths[chromosome_rank])))
  }
  DBI::dbExecute(connection, sprintf("COPY (SELECT * FROM read_parquet([%s])) TO %s (FORMAT parquet, COMPRESSION zstd)",
                                     paste(DBI::dbQuoteString(connection, piece_paths), collapse = ", "),
                                     DBI::dbQuoteString(connection, position_path)))
  invisible(NULL)
}
