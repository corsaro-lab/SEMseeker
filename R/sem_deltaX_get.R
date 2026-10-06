# ---------------------------------------------------------------------------
# sem_deltaX_get() - polars-native, wide-dataframe derived markers
#
# Written 2026-06-01 to replace an R-loop implementation, which was kept alongside it
# for reference and A/B testing and removed on 2026-10-02 once nothing had called it
# for four months. It is in the git history if the comparison is ever wanted again.
#
# Builds the derived pivot tables (DELTAQ, DELTARQ, DELTAP, DELTARP) from the
# corresponding source pivots (DELTAS or DELTAR) entirely in polars.
#
# Invariants (confirmed with user 2026-06-01):
#   - breaks GLOBAL on the union of HYPER+HYPO sample values (not per-figure)
#   - zero == missing (excluded from stats, returned NA in the output)
#   - right=TRUE intervals (a, b] - matches polars cut default and R cut
#     default
#   - skip per-sample bedgraph dumps for derived markers (the removed R-loop
#     implementation wrote 4013 single bedgraphs per derived marker x figure - the
#     ~9 min wall-clock bottleneck on ewas_data_hub)
#
# Implementation notes:
#   - For markers ending in 'P' (DELTAP, DELTARP): the bin edges are
#     equal-width on [global_min, global_max]. Only the per-column min/max
#     are pulled out - O(N_samples) RAM, streaming-friendly. No unpivot.
#   - For markers ending in 'Q' (DELTAQ, DELTARQ): quantile bin edges need
#     the full value distribution; we collect the two source DataFrames and
#     compute the quantiles in R (R `quantile` is C-vectorised and exact;
#     polars `quantile` on a single big column also works but the unpivot
#     in polars R 1.11 can drop class info on lazy frames).
#   - Cut is applied column-wise inside a single $with_columns(): polars
#     Rust runtime parallelises the expressions via Rayon, so the cost
#     stays wide even though we list the columns from R.
# ---------------------------------------------------------------------------

#' Compute derived markers (DELTAQ/DELTARQ/DELTAP/DELTARP) direct from
#' source pivots via polars wide expressions.
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

  # Replace 0 with null (matches `vec[vec == 0] <- NA` in legacy sem_deltaX_get).
  zero_to_null <- function(col_name) {
    polars::pl$when(polars::pl$col(col_name) == 0)$
      then(polars::pl$lit(NA_real_))$
      otherwise(polars::pl$col(col_name))
  }

  # Per-column min/max ignoring nulls - pulled into R as scalars, then
  # global aggregated. Memory cost: O(N_samples) doubles per source figure.
  # NB: polars R $select() wants expressions as varargs, so we splice the
  # list via do.call().
  col_minmax <- function(lf, cols) {
    mins_expr <- lapply(cols, function(c) zero_to_null(c)$min()$alias(c))
    maxs_expr <- lapply(cols, function(c) zero_to_null(c)$max()$alias(c))
    mins_row <- as.data.frame(do.call(lf$select, mins_expr)$collect())
    maxs_row <- as.data.frame(do.call(lf$select, maxs_expr)$collect())
    list(
      min = suppressWarnings(min(as.numeric(mins_row[1, ]), na.rm = TRUE)),
      max = suppressWarnings(max(as.numeric(maxs_row[1, ]), na.rm = TRUE))
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
                " [deltaX_get_polars] skip ", mar, " (both pivots already exist)")
      next
    }
    if (!file.exists(src_h) || !file.exists(src_o)) {
      core_log_event("WARNING: ", Sys.time(),
                " [deltaX_get_polars] source pivot(s) missing for marker=", mar)
      next
    }

    t_start <- Sys.time()
    core_log_event("INFO: ", t_start,
              " [deltaX_get_polars] computing ", mar,
              " (Q=", Q_val, ", suffix=",
              if (endsWith(mar, "P")) "P/equal-width" else "Q/quantile",
              ") direct from ", src, " pivots")

    # read via unified dispatcher. The file.exists() guard above
    # has already filtered out cases where neither source pivot exists.
    lf_h <- io_read_pivot(src, "HYPER", area, subarea)
    lf_o <- io_read_pivot(src, "HYPO",  area, subarea)
    cols_h <- setdiff(names(lf_h$collect_schema()), c("CHR", "START", "END"))
    cols_o <- setdiff(names(lf_o$collect_schema()), c("CHR", "START", "END"))

    # ---------------------------------------------------------------------
    # Compute global breaks on HYPER + HYPO (0 → null).
    # ---------------------------------------------------------------------
    if (endsWith(mar, "P")) {
      # Equal-width: O(N_samples) RAM.
      mm_h <- col_minmax(lf_h, cols_h)
      mm_o <- col_minmax(lf_o, cols_o)
      mn <- min(mm_h$min, mm_o$min, na.rm = TRUE)
      mx <- max(mm_h$max, mm_o$max, na.rm = TRUE)
      if (!is.finite(mn) || !is.finite(mx) || mn == mx) {
        core_log_event("WARNING: ", Sys.time(),
                  " [deltaX_get_polars] degenerate range for ", mar,
                  ": min=", mn, " max=", mx, " - skip")
        next
      }
      breaks_full <- seq(mn, mx, length.out = Q_val + 1L)

    } else if (endsWith(mar, "Q")) {
      # Quantile: need the full distribution. Collect both DataFrames and
      # flatten - costly in RAM (~ N_samples × N_probes × 8 bytes per
      # source figure) but exact.
      probs <- seq(0, 1, length.out = Q_val + 1L)
      df_h  <- as.data.frame(lf_h$select(cols_h)$collect())
      df_o  <- as.data.frame(lf_o$select(cols_o)$collect())
      vals  <- c(unlist(df_h, use.names = FALSE),
                 unlist(df_o, use.names = FALSE))
      rm(df_h, df_o); gc(verbose = FALSE)
      vals[vals == 0] <- NA_real_
      breaks_full <- unique(stats::quantile(vals, probs = probs, na.rm = TRUE))
      rm(vals); gc(verbose = FALSE)

    } else {
      core_log_event("WARNING: ", Sys.time(),
                " [deltaX_get_polars] marker ", mar,
                " doesn't end in P or Q - skipping")
      next
    }

    core_log_event("INFO: ", Sys.time(),
              " [deltaX_get_polars] ", mar,
              " breaks=[", paste(signif(breaks_full, 4), collapse = ", "), "]")

    # Inner breakpoints (polars cut takes cut points, not endpoints).
    inner_breaks <- breaks_full[-c(1L, length(breaks_full))]
    if (length(inner_breaks) == 0L) {
      core_log_event("WARNING: ", Sys.time(),
                " [deltaX_get_polars] no inner breakpoints for ", mar,
                " - skip")
      next
    }

    # ---------------------------------------------------------------------
    # Apply cut wide: 0→null, then polars cut → int label 1..N.
    # ---------------------------------------------------------------------
    make_cut_expr <- function(col_name) {
      zero_to_null(col_name)$
        cut(breaks = inner_breaks)$
        to_physical()$
        cast(polars::pl$Int32)$
        add(1L)$
        alias(col_name)
    }

    write_one <- function(lf_src, cols, dest_path) {
      exprs <- lapply(cols, make_cut_expr)
      out <- do.call(lf_src$with_columns, exprs)$
        sort(c("CHR", "START"), descending = FALSE)
      dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
      out$sink_parquet(dest_path)
    }

    write_one(lf_h, cols_h, dest_h)
    write_one(lf_o, cols_o, dest_o)

    dt <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    core_log_event("INFO: ", Sys.time(),
              " [deltaX_get_polars] wrote ", mar,
              " HYPER+HYPO in ", round(dt, 1), " sec")
  }

  invisible()
}
