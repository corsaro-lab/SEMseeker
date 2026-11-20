#' Compute the population markers of every sample at once, from the pivots
#'
#' The per-sample engine writes one bedGraph per sample, marker and figure and
#' merges them back into pivots, which on 4 000 samples took about 11 hours.
#' This engine computes the same markers column by column on the SIGNAL
#' POSITION pivot that io_signal_save() wrote, and writes each pivot once:
#'
#'   DELTAS_HYPER = max(SIGNAL - high threshold, 0)
#'   DELTAS_HYPO  = max(low threshold - SIGNAL, 0)
#'   MUTATIONS    = DELTAS > 0, as an integer
#'   DELTAR       = DELTAS / max(high - low, 1e-9)
#'   LESIONS      = the binomial window test of sem_lesions_get_bulk() on the
#'                  MUTATIONS pivot, 200 samples at a time
#'
#' A missing signal stays missing in every marker. The two engines must give
#' the same values: test-e2e_population_path_equivalence.R checks it. A pivot
#' that already exists is not written again, and when all of them exist the
#' function returns before reading anything.
#'
#' @param signal_data kept for the signature of sem_analyze_population();
#'   not used, the signal is read from the SIGNAL POSITION pivot
#' @param sample_sheet kept for the same signature; not used, the samples are
#'   the columns of the SIGNAL pivot
#' @param signal_thresholds data.frame with CHR, START, END,
#'   signal_superior_thresholds (high) and signal_inferior_thresholds (low),
#'   one row per position
#' @param probe_features kept for the same signature; not used
#' @return NULL, invisibly; the side effect is the marker pivots
#'
#' @keywords internal
#' @noRd
sem_analyze_population_bulk <- function(signal_data, sample_sheet,
                                    signal_thresholds, probe_features) {

  ssEnv <- core_get_session_info()
  start_time <- Sys.time()
  core_log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"),
            " [apb] FRAME ENTERED mem_MB=", round(sum(gc()[, "(Mb)"]), 1))
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " [sem_analyze_population_bulk] start (vectorized)")

  # ---- Step 1: prepare paths and inputs ------------------------------------
  pivots_dir <- file.path(ssEnv$result_folderData, "Pivots")
  dir.create(pivots_dir, recursive = TRUE, showWarnings = FALSE)
  for (m in c("SIGNAL", "DELTAS", "DELTAR", "MUTATIONS", "LESIONS")) {
    dir.create(file.path(pivots_dir, m), recursive = TRUE, showWarnings = FALSE)
  }

  # 2026-06-09: EARLY-RETURN if every destination pivot already
  # exists. The per-figure skip checks further down inside this function
  # do guard the actual compute, but the SETUP between [start] and the
  # first per-figure block - io_read_pivot SIGNAL + collect_schema(4014
  # cols) + with_columns(cast Categorical→String) + lazy join - still
  # runs every call. On ewas-scale (367k × 4013) that setup alone
  # peaked ~30 GB R+Polars even in pure-resume mode (v25/v26/v27/v28
  # all crashed there with all downstream pivots already on disk).
  # When the whole function has nothing to produce, skipping it
  # entirely is the only path that scales.
  all_destinations <- character(0)
  for (marker in c("DELTAS", "DELTAR", "MUTATIONS", "LESIONS")) {
    for (figure in c("HYPER", "HYPO")) {
      all_destinations <- c(all_destinations,
        io_pivot_file_name_parquet(marker, figure, "POSITION", "WHOLE"))
    }
  }
  if (length(all_destinations) > 0L && all(file.exists(all_destinations))) {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " [sem_analyze_population_bulk] ALL ", length(all_destinations),
              " destination pivots already exist - skipping bulk pass entirely.")
    return(invisible(NULL))
  }

  signal_pivot_path <- io_pivot_file_name_parquet("SIGNAL", io_signal_figure(), "POSITION", "WHOLE")
  if (!file.exists(signal_pivot_path)) {
    stop("[sem_analyze_population_bulk] SIGNAL POSITION pivot mancante: ",
         signal_pivot_path,
         " - io_signal_save() deve essere chiamato prima.")
  }

  # The thresholds join the signal on the coordinates, as a view of the data
  # engine. CHR as text and START/END as integers, the types of the pivot.
  connection <- io_engine_connection()
  thresholds_view <- paste0("thresholds_", gsub("[^A-Za-z0-9]", "", basename(tempfile(""))))
  duckdb::duckdb_register(connection, thresholds_view, data.frame(
    CHR   = as.character(signal_thresholds$CHR),
    START = as.integer(signal_thresholds$START),
    END   = as.integer(signal_thresholds$END),
    .HIGH = signal_thresholds$signal_superior_thresholds,
    .LOW  = signal_thresholds$signal_inferior_thresholds,
    stringsAsFactors = FALSE,
    check.names      = FALSE
  ))
  on.exit(duckdb::duckdb_unregister(connection, thresholds_view), add = TRUE)
  n_thr_positions <- nrow(signal_thresholds)
  rm(signal_thresholds)
  invisible(gc(verbose = FALSE))

  # Paths and sample ids are free text and enter the SQL quoted. path: a file
  # path, quoted as an SQL string; name: a column name, quoted as an SQL
  # identifier.
  quote_text <- function(path) DBI::dbQuoteString(connection, path)
  quote_name <- function(name) as.character(DBI::dbQuoteIdentifier(connection, name))
  # The column names of a parquet file, read from its footer without reading
  # the data. path: the parquet file.
  pivot_columns <- function(path)
    DBI::dbGetQuery(connection, sprintf("DESCRIBE SELECT * FROM read_parquet(%s)", quote_text(path)))$column_name
  # Writes `select_list` over `source_query` to `path`, keeping the row order
  # of the source. select_list: the SQL expressions after SELECT; source_query:
  # the FROM clause; path: the parquet file.
  write_pivot <- function(select_list, source_query, path)
    DBI::dbExecute(connection, sprintf("COPY (SELECT %s FROM %s) TO %s (FORMAT parquet, COMPRESSION zstd)",
                                       select_list, source_query, quote_text(path)))
  # A NULL stays NULL: DuckDB's greatest() would turn greatest(NULL, 0) into
  # 0, a missing signal into "no deviation". value_sql: the SQL of the value.
  floor_at_zero <- function(value_sql)
    sprintf("CASE WHEN (%1$s) IS NULL THEN NULL ELSE greatest(%1$s, 0) END", value_sql)

  signal_pivot_columns <- pivot_columns(signal_pivot_path)
  sample_cols <- setdiff(signal_pivot_columns, c("CHR", "START", "END", ".HIGH", ".LOW"))
  coordinates_sql <- 'signal.CHR, signal.START, signal."END"'
  signal_join <- sprintf('read_parquet(%s) AS signal JOIN %s AS thresholds
                          ON CAST(signal.CHR AS VARCHAR) = thresholds.CHR AND signal.START = thresholds.START AND signal."END" = thresholds."END"',
                         quote_text(signal_pivot_path), thresholds_view)

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " [sem_analyze_population_bulk] joined SIGNAL with thresholds: ",
            length(sample_cols), " sample columns, ",
            n_thr_positions, " positions")

  # ---- Step 2: DELTAS_HYPER, DELTAS_HYPO bulk -----------------------------
  for (figure in c("HYPER", "HYPO")) {
    deltas_path <- io_pivot_file_name_parquet("DELTAS", figure, "POSITION", "WHOLE")
    if (file.exists(deltas_path)) {
      core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " [bulk] DELTAS_", figure, " pivot already exists, skip")
      next
    }
    t0 <- Sys.time()
    deltas_sql <- vapply(sample_cols, function(sample_name) {
      column <- paste0("signal.", quote_name(sample_name))
      difference <- if (figure == "HYPER") paste(column, "- thresholds.\".HIGH\"") else paste("thresholds.\".LOW\" -", column)
      paste(floor_at_zero(difference), "AS", quote_name(sample_name))
    }, character(1))
    write_pivot(paste(c(coordinates_sql, deltas_sql), collapse = ", "), signal_join, deltas_path)
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " [bulk] DELTAS_", figure, " written in ",
              round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " sec")
  }

  # ---- Step 3: MUTATIONS_{HYPER,HYPO} bulk = (DELTAS > 0) as INTEGER -----
  for (figure in c("HYPER", "HYPO")) {
    mut_path <- io_pivot_file_name_parquet("MUTATIONS", figure, "POSITION", "WHOLE")
    if (file.exists(mut_path)) {
      core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " [bulk] MUTATIONS_", figure, " pivot already exists, skip")
      next
    }
    t0 <- Sys.time()
    deltas_path <- io_pivot_file_name_parquet("DELTAS", figure, "POSITION", "WHOLE")
    mutations_sql <- vapply(sample_cols, function(sample_name)
      sprintf("CAST(%1$s > 0 AS INTEGER) AS %1$s", quote_name(sample_name)), character(1))
    write_pivot(paste(c('CHR, START, "END"', mutations_sql), collapse = ", "),
                sprintf("read_parquet(%s)", quote_text(deltas_path)), mut_path)
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " [bulk] MUTATIONS_", figure, " written in ",
              round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " sec")
  }

  # ---- Step 4: DELTAR_{HYPER,HYPO} bulk = DELTAS / (high - low) ----------
  for (figure in c("HYPER", "HYPO")) {
    deltar_path <- io_pivot_file_name_parquet("DELTAR", figure, "POSITION", "WHOLE")
    if (file.exists(deltar_path)) {
      core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " [bulk] DELTAR_", figure, " pivot already exists, skip")
      next
    }
    t0 <- Sys.time()
    deltas_path <- io_pivot_file_name_parquet("DELTAS", figure, "POSITION", "WHOLE")
    deltas_join <- sprintf('read_parquet(%s) AS signal JOIN %s AS thresholds
                            ON CAST(signal.CHR AS VARCHAR) = thresholds.CHR AND signal.START = thresholds.START AND signal."END" = thresholds."END"',
                           quote_text(deltas_path), thresholds_view)
    # the range of the reference, floored at 1e-9 so a constant probe does
    # not divide by zero, as in sem_deltar_single_sample()
    range_sql <- 'CASE WHEN (thresholds.".HIGH" - thresholds.".LOW") IS NULL THEN NULL ELSE greatest(thresholds.".HIGH" - thresholds.".LOW", 1e-9) END'
    deltar_sql <- vapply(sample_cols, function(sample_name)
      sprintf("signal.%1$s / (%2$s) AS %1$s", quote_name(sample_name), range_sql), character(1))
    write_pivot(paste(c(coordinates_sql, deltar_sql), collapse = ", "), deltas_join, deltar_path)
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " [bulk] DELTAR_", figure, " written in ",
              round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " sec")
  }

  # ---- Step 5: LESIONS_{HYPER,HYPO} per-sample column processing ---------
  # LESIONS calculation via sem_lesions_get_bulk() (multi-sample, bp-window):
  # finestra fisica in bp, soppianta la vecchia logica row-count
  # (sliding_window_size). Il counterpart single-sample e' sem_lesions_get()
  # in R/lesions_get.R (usato dal path legacy sem_analyze_population per-sample).
  # Storia: (kbp arg) → (LESIONS_BP ssEnv) merged 2026-06-10
  # con default 5000 bp (literature-aligned, revisione metodologica in sospeso).
  lesions_bp     <- as.integer(ssEnv$LESIONS_BP)
  bonf_threshold <- as.numeric(ssEnv$bonferroni_threshold)
  CHUNK_SAMPLES  <- 200L  # gruppi di sample per limitare RAM

  for (figure in c("HYPER", "HYPO")) {
    lesions_path <- io_pivot_file_name_parquet("LESIONS", figure, "POSITION", "WHOLE")
    if (file.exists(lesions_path)) {
      core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " [bulk] LESIONS_", figure, " pivot already exists, skip")
      next
    }
    t0 <- Sys.time()
    mut_path <- io_pivot_file_name_parquet("MUTATIONS", figure, "POSITION", "WHOLE")

    sample_chunks <- split(sample_cols,
                          ceiling(seq_along(sample_cols) / CHUNK_SAMPLES))
    tmp_dir <- tempfile("lesions_chunks_")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

    for (ci in seq_along(sample_chunks)) {
      cols_this <- sample_chunks[[ci]]
      # Carica solo le colonne di questo chunk + coordinate (ordinato per CHR/START)
      mut_df <- DBI::dbGetQuery(connection, sprintf('SELECT CHR, START, "END", %s FROM read_parquet(%s) ORDER BY CHR, START',
                                                    paste(vapply(cols_this, quote_name, character(1)), collapse = ", "),
                                                    quote_text(mut_path)))
      mut_df$CHR <- as.character(mut_df$CHR)

      # Delega LESIONS computation a sem_lesions_get_bulk (multi-sample bp window).
      # Il counterpart single-sample (R/lesions_get.R) e' usato dal path
      # legacy sem_analyze_population per-sample loop.
      les_mat <- sem_lesions_get_bulk(mut_df, cols_this,
                                  LESIONS_BP = lesions_bp,
                                  bonf_threshold = bonf_threshold)

      les_df <- data.frame(
        CHR   = mut_df$CHR,
        START = mut_df$START,
        END   = mut_df$END,
        les_mat,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      chunk_file <- file.path(tmp_dir, sprintf("%03d.parquet", ci))
      io_write_frame_parquet(les_df, chunk_file)
      rm(mut_df, les_mat, les_df); invisible(gc(verbose = FALSE))
    }

    # Concat dei chunk per colonna (tutti hanno stesso ordine righe).
    chunk_paths <- list.files(tmp_dir, pattern = "\\.parquet$", full.names = TRUE)
    if (length(chunk_paths) == 1L) {
      file.copy(chunk_paths[1], lesions_path, overwrite = TRUE)
    } else {
      # Every chunk holds the same positions: the samples are joined back on
      # the coordinates, in the order of the chunks.
      chunk_sources <- sprintf("read_parquet(%s) AS chunk%d", vapply(chunk_paths, quote_text, character(1)), seq_along(chunk_paths))
      joined <- paste(c(chunk_sources[1], sprintf('JOIN %s USING (CHR, START, "END")', chunk_sources[-1])), collapse = " ")
      write_pivot("*", joined, lesions_path)
    }

    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " [bulk] LESIONS_", figure, " written in ",
              round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
              " min (",  length(sample_chunks), " chunks of ",
              CHUNK_SAMPLES, " samples, LESIONS_BP=", lesions_bp, ")")
  }

  # ============================================================
  # LEGACY ROW-WINDOW LESIONS CODE (v1, soppiantata da
  # lesions_new il 2026-06-08). Commentato per reference;
  # NON eseguito. Eliminabile in cleanup futuro.
  # ============================================================
  # Algoritmo legacy (row-count window via sliding_window_size):
  #   ENRICHMENT      = rolling_sum(MUT, sliding_window_size, center)
  #   BASEPAIR_COUNT  = rolling_max(START) - rolling_min(START) (variabile per riga)
  #   p0              = sum(MUT) / length(MUT) per-CHR
  #   lesionpValue    = pbinom(ENRICHMENT-1, size=W, prob=p0, lower=FALSE)
  #   LESIONS_BOOL    = lesionpValue < bonferroni / (n_probes * log10(BASEPAIR_COUNT))
  #
  # Differenza chiave vs lesions_new: window in N righe vs ±X kbp fisici.
  # Il row-window non e' biologicamente coerente perche' span fisico varia
  # 400x fra CpG islands e regioni intergeniche.

  end_time <- Sys.time()
  total_min <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " [sem_analyze_population_bulk] completed in ", total_min, " min")
  invisible(NULL)
}
