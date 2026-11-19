# ---------------------------------------------------------------------------
# sem_deltaX_get() - derived markers computed column by column on the pivots
#
# Written 2026-06-01 to replace an R-loop implementation, which was kept alongside it
# for reference and A/B testing and removed on 2026-10-02 once nothing had called it
# for four months. It is in the git history if the comparison is ever wanted again.
#
# Builds the derived pivot tables (DELTAQ, DELTARQ, DELTAP, DELTARP) from the
# corresponding source pivots (DELTAS or DELTAR) in the data engine, without
# bringing the pivots into R (the quantile breaks excepted, see below).
#
# Invariants (confirmed with user 2026-06-01):
#   - breaks GLOBAL on the union of HYPER+HYPO sample values (not per-figure)
#   - zero == missing (excluded from stats, returned NA in the output)
#   - right=TRUE intervals (a, b] - the polars cut default and the R cut
#     default, kept when the engine moved to DuckDB
#   - skip per-sample bedgraph dumps for derived markers (the removed R-loop
#     implementation wrote 4013 single bedgraphs per derived marker x figure - the
#     ~9 min wall-clock bottleneck on ewas_data_hub)
#
# Implementation notes:
#   - For markers ending in 'P' (DELTAP, DELTARP): the bin edges are
#     equal-width on [global_min, global_max]. Only the per-column min/max
#     are pulled out - O(N_samples) RAM, streaming-friendly. No unpivot.
#   - For markers ending in 'Q' (DELTAQ, DELTARQ): quantile bin edges need
#     the full value distribution; the two source pivots are collected and
#     the quantiles computed by R's `quantile`, as before the move to DuckDB,
#     so the breaks do not depend on the engine.
#   - The cut is one CASE per sample column, inside a single COPY.
# ---------------------------------------------------------------------------

#' Compute derived markers (DELTAQ/DELTARQ/DELTAP/DELTARP) direct from
#' source pivots, one SQL expression per sample column.
#'
#' @param markers Optional character vector to restrict which derived markers
#'   to materialise. Default \code{NULL} = all eligible markers from
#'   \code{ssEnv$keys_markers_figures}.
#'
#' @return Invisible \code{NULL}. Side effect: writes
#'   \code{<MARKER>_HYPER_POSITION_WHOLE_HG19.parquet} and
#'   \code{<MARKER>_HYPO_POSITION_WHOLE_HG19.parquet} for each eligible
#'   marker. Per-sample bedgraph files are NOT written.
#'
#' @keywords internal
#' @noRd
sem_deltaX_get <- function(markers = NULL) {

  ssEnv <- core_get_session_info()
  keys  <- ssEnv$keys_markers_figures
  keys  <- keys[, !(colnames(keys) %in% c("FIGURE", "COMBINED")), drop = FALSE]
  keys  <- unique(keys)
  keys  <- keys[which(!is.na(keys$SOURCE) & !is.na(keys$Q) & keys$Q != 1), , drop = FALSE]
  if (!is.null(markers)) keys <- keys[which(keys$MARKER %in% markers), , drop = FALSE]
  if (nrow(keys) == 0L) return(invisible())

  area    <- "POSITION"
  subarea <- "WHOLE"

  connection <- io_engine_connection()
  # Paths and sample ids are free text and enter the SQL quoted. path: a file
  # path; name: a column name.
  quote_text <- function(path) DBI::dbQuoteString(connection, path)
  quote_name <- function(name) as.character(DBI::dbQuoteIdentifier(connection, name))
  # Zero means "no deviation" and is excluded from the breaks and from the
  # output, as `vec[vec == 0] <- NA` did in the R-loop version. column_sql:
  # a quoted column name.
  zero_to_null <- function(column_sql) sprintf("CASE WHEN %1$s = 0 THEN NULL ELSE %1$s END", column_sql)
  # The sample columns of a source pivot. path: the parquet file.
  sample_columns <- function(path)
    setdiff(DBI::dbGetQuery(connection, sprintf("DESCRIBE SELECT * FROM read_parquet(%s)", quote_text(path)))$column_name,
            c("CHR", "START", "END"))

  # Smallest and largest non-zero value over all the sample columns of a
  # pivot, as two scalars: only one row of aggregates enters R.
  # path: the parquet file; cols: its sample columns.
  col_minmax <- function(path, cols) {
    quoted <- vapply(cols, quote_name, character(1))
    aggregates <- c(sprintf("min(%s)", zero_to_null(quoted)), sprintf("max(%s)", zero_to_null(quoted)))
    extremes <- as.numeric(DBI::dbGetQuery(connection, sprintf("SELECT %s FROM read_parquet(%s)",
                                                               paste(aggregates, collapse = ", "), quote_text(path)))[1, ])
    list(
      min = suppressWarnings(min(extremes[seq_along(cols)], na.rm = TRUE)),
      max = suppressWarnings(max(extremes[length(cols) + seq_along(cols)], na.rm = TRUE))
    )
  }

  for (k in seq_len(nrow(keys))) {

    key   <- keys[k, ]
    src   <- as.character(key$SOURCE)
    mar   <- as.character(key$MARKER)
    Q_val <- as.integer(key$Q)

    dest_h <- io_pivot_file_name_parquet(mar, "HYPER", area, subarea)
    dest_o <- io_pivot_file_name_parquet(mar, "HYPO",  area, subarea)
    src_h  <- io_pivot_file_name_parquet(src, "HYPER", area, subarea)
    src_o  <- io_pivot_file_name_parquet(src, "HYPO",  area, subarea)

    if (file.exists(dest_h) && file.exists(dest_o)) {
      core_log_event("INFO: ", Sys.time(),
                " [deltaX_get] skip ", mar, " (both pivots already exist)")
      next
    }
    if (!file.exists(src_h) || !file.exists(src_o)) {
      core_log_event("WARNING: ", Sys.time(),
                " [deltaX_get] source pivot(s) missing for marker=", mar)
      next
    }

    t_start <- Sys.time()
    core_log_event("INFO: ", t_start,
              " [deltaX_get] computing ", mar,
              " (Q=", Q_val, ", suffix=",
              if (endsWith(mar, "P")) "P/equal-width" else "Q/quantile",
              ") direct from ", src, " pivots")

    cols_h <- sample_columns(src_h)
    cols_o <- sample_columns(src_o)

    # ---------------------------------------------------------------------
    # Compute global breaks on HYPER + HYPO (0 → null).
    # ---------------------------------------------------------------------
    if (endsWith(mar, "P")) {
      # Equal-width: O(N_samples) RAM.
      mm_h <- col_minmax(src_h, cols_h)
      mm_o <- col_minmax(src_o, cols_o)
      mn <- min(mm_h$min, mm_o$min, na.rm = TRUE)
      mx <- max(mm_h$max, mm_o$max, na.rm = TRUE)
      if (!is.finite(mn) || !is.finite(mx) || mn == mx) {
        core_log_event("WARNING: ", Sys.time(),
                  " [deltaX_get] degenerate range for ", mar,
                  ": min=", mn, " max=", mx, " - skip")
        next
      }
      breaks_full <- seq(mn, mx, length.out = Q_val + 1L)

    } else if (endsWith(mar, "Q")) {
      # Quantile: need the full distribution. Collect both DataFrames and
      # flatten - costly in RAM (~ N_samples × N_probes × 8 bytes per
      # source figure) but exact.
      probs <- seq(0, 1, length.out = Q_val + 1L)
      read_values <- function(path, cols)
        # path: a source pivot; cols: its sample columns, read whole
        DBI::dbGetQuery(connection, sprintf("SELECT %s FROM read_parquet(%s)",
                                            paste(vapply(cols, quote_name, character(1)), collapse = ", "), quote_text(path)))
      df_h  <- read_values(src_h, cols_h)
      df_o  <- read_values(src_o, cols_o)
      vals  <- c(unlist(df_h, use.names = FALSE),
                 unlist(df_o, use.names = FALSE))
      rm(df_h, df_o); gc(verbose = FALSE)
      vals[vals == 0] <- NA_real_
      breaks_full <- unique(stats::quantile(vals, probs = probs, na.rm = TRUE))
      rm(vals); gc(verbose = FALSE)

    } else {
      core_log_event("WARNING: ", Sys.time(),
                " [deltaX_get] marker ", mar,
                " doesn't end in P or Q - skipping")
      next
    }

    core_log_event("INFO: ", Sys.time(),
              " [deltaX_get] ", mar,
              " breaks=[", paste(signif(breaks_full, 4), collapse = ", "), "]")

    # Inner breakpoints: the cut compares with the points between the bins, not the endpoints.
    inner_breaks <- breaks_full[-c(1L, length(breaks_full))]
    if (length(inner_breaks) == 0L) {
      core_log_event("WARNING: ", Sys.time(),
                " [deltaX_get] no inner breakpoints for ", mar,
                " - skip")
      next
    }

    # ---------------------------------------------------------------------
    # Apply the cut: 0 and NULL give NULL; otherwise the label is 1 plus the
    # number of inner breaks the value exceeds, which makes the intervals
    # (a, b] as in polars' cut and R's cut. The breaks enter the SQL as text
    # cast to DOUBLE with 17 significant digits, so each is the same double
    # R computed.
    # ---------------------------------------------------------------------
    break_sql <- sprintf("CAST('%s' AS DOUBLE)", sprintf("%.17g", inner_breaks))
    # col_name: a sample column; returns its label expression, aliased to it.
    make_cut_expr <- function(col_name) {
      quoted <- quote_name(col_name)
      exceeded <- paste(sprintf("CASE WHEN %s > %s THEN 1 ELSE 0 END", quoted, break_sql), collapse = " + ")
      sprintf("CASE WHEN %1$s IS NULL OR %1$s = 0 THEN NULL ELSE CAST(1 + %2$s AS INTEGER) END AS %1$s", quoted, exceeded)
    }

    # Writes the labels of one figure, sorted on CHR and START. source_path:
    # the source pivot; cols: its sample columns; dest_path: the derived pivot.
    write_one <- function(source_path, cols, dest_path) {
      dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
      DBI::dbExecute(connection, sprintf('COPY (SELECT CHR, START, "END", %s FROM read_parquet(%s) ORDER BY CHR, START) TO %s (FORMAT parquet, COMPRESSION zstd)',
                                         paste(vapply(cols, make_cut_expr, character(1)), collapse = ", "),
                                         quote_text(source_path), quote_text(dest_path)))
    }

    write_one(src_h, cols_h, dest_h)
    write_one(src_o, cols_o, dest_o)

    dt <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    core_log_event("INFO: ", Sys.time(),
              " [deltaX_get] wrote ", mar,
              " HYPER+HYPO in ", round(dt, 1), " sec")
  }

  invisible()
}
