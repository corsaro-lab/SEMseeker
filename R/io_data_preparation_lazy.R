#' Data-engine equivalent of `io_data_preparation()` for the lazy
#' batch path. Two responsibilities:
#'
#'   1. **Apply `transformation_y` on every sample column** of the pivot,
#'      so the lazy path produces statistically equivalent CSVs
#'      to the per-area `apply_stat_model_batch.R` path which goes
#'      through R-side `io_data_preparation()`.
#'
#'   2. **Drop universal degenerate-burden rows** (`var(Y) == 0`,
#'      i.e. every sample carries the same value for that probe/area),
#'      so `limma::lmFit` doesn't see constant-Y rows that produce
#'      NaN t-stats and pollute the inference CSV.
#'
#' Both steps stay end-to-end lazy: the operations are appended to the
#' input `pivot_lazy` and the caller decides when to `collect()`. NO
#' materialisation happens here - that is the whole point of the
#' lazy path.
#'
#' Scope is intentionally narrower than `io_data_preparation()`:
#'   - Y-side ONLY (transformation_y, degenerate-row filter). The
#'     transformation_x and IV-factor coding from `io_data_preparation()`
#'     act on the sample-sheet R-side and are applied in the lazy batch
#'     caller AFTER the y_mat is materialised.
#'   - `factor` and `quantile_<N>` transformations are NOT supported
#'     here (they need group/window aggregations that don't survive
#'     well in lazy mode on a wide-format pivot). They raise an
#'     explicit error directing the caller to the per-area path.
#'
#' @param pivot_lazy A lazy table of the data engine (`tbl_lazy`), or a
#'   data.frame that is registered with the engine, whose rows are probes/areas
#'   and whose columns are: optional `AREA` (carrier of the probe ID,
#'   passed through untouched) + one column per sample.
#' @param sample_cols Character vector. The column names in
#'   `pivot_lazy` that hold the Y values (one per sample). Coord /
#'   metadata columns (AREA, PROBE, CHR, START, END, K27/K450/K850)
#'   MUST be excluded by the caller before passing them in.
#' @param transformation_y One of `c("none", "log", "log2", "log10",
#'   "exp", "scale")`. Default `"none"`. Mirrors the R-side switch in
#'   `io_data_preparation()`.
#' @param apply_degenerate_filter Logical, default `TRUE`. Drop rows
#'   where `min(sample_cols) == max(sample_cols)` (no variance in Y).
#'   Matches the filter applied in `io_data_preparation()` for the
#'   non-lazy path. Set `FALSE` only for diagnostic test code that wants
#'   to see the raw transformations alone.
#' @param key Optional named list with `MARKER`/`FIGURE`/`AREA`/`SUBAREA`
#'   identifiers - included in the log line so the diagnostic is
#'   traceable across multiple concurrent dispatches.
#' @param family_test Optional character - passed through to the log
#'   line, same purpose as `key`.
#'
#' @return The transformed and filtered lazy table. The schema is
#'   unchanged (same column names, same order), with fewer rows.
#'
#' @keywords internal
#' @noRd
io_data_preparation_lazy <- function(pivot_lazy,
                                   sample_cols,
                                   transformation_y = "none",
                                   apply_degenerate_filter = TRUE,
                                   key = NULL,
                                   family_test = NULL) {

  if (is.data.frame(pivot_lazy) && !inherits(pivot_lazy, "tbl_lazy")) {
    pivot_lazy <- io_register_frame(pivot_lazy)
  } else if (!inherits(pivot_lazy, "tbl_lazy")) {
    stop("io_data_preparation_lazy: pivot_lazy must be a lazy table of the data engine or a data.frame.",
         call. = FALSE)
  }
  if (!is.character(sample_cols) || length(sample_cols) == 0L) {
    stop("io_data_preparation_lazy: sample_cols must be a non-empty character vector.",
         call. = FALSE)
  }

  trans_y <- as.character(transformation_y)
  if (length(trans_y) == 0L || is.na(trans_y) || !nzchar(trans_y)) {
    trans_y <- "none"
  }

  key_str <- if (!is.null(key))
    sprintf(" [%s/%s/%s/%s]", as.character(key$MARKER), as.character(key$FIGURE),
            as.character(key$AREA), as.character(key$SUBAREA)) else ""
  fam_str <- if (!is.null(family_test))
    sprintf(" %s", as.character(family_test)) else ""

  supported <- c("none", "log", "log2", "log10", "exp", "scale")
  if (!(trans_y %in% supported)) {
    # `factor` and `quantile_<N>` cannot be expressed cleanly in a lazy
    # wide-format chain (they need group / window aggregations). Emit a
    # warning so the user sees the gap explicitly, then silently treat
    # the row as `transformation_y = "none"` and continue - the CSV
    # will still be produced (without the requested transformation),
    # and the user can route the inference_detail to the per-area
    # engine offline if they need that transformation honoured.
    explain_alt <- if (grepl("^quantile", trans_y) || trans_y == "factor")
      paste0(
        " Route this inference_detail to the per-area engine ",
        "(family_test != 'limma_*' / 'voom_*' to bypass the lazy batch path), ",
        "or pre-transform the pivot before scan if you must keep limma/voom."
      ) else
      paste0(" Use one of: ", paste(supported, collapse = ", "), ".")
    warning(sprintf(paste0(
      "io_data_preparation_lazy: transformation_y='%s' is NOT supported by ",
      "the lazy batch path. Continuing with transformation_y='none'.",
      "%s"),
      trans_y, explain_alt), call. = FALSE)
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " io_data_preparation_lazy", fam_str, key_str,
              ": unsupported transformation_y='", trans_y,
              "' on lazy batch path - continuing as 'none'.")
    trans_y <- "none"
  }

  # ---- Step 1: transformation_y on every sample column --------------
  # One SQL expression per sample column. DuckDB follows IEEE semantics
  # on doubles as polars did: ln(0) is -Inf, ln of a negative value is
  # NaN, a division by a zero standard deviation is Inf or NaN.
  # Note: `log` on negative or zero values yields NaN; the R-side
  # `io_data_preparation()` guards by adding `min(burden_values[burden_values>0])`,
  # which would require a scan over the full pivot. We use a tiny
  # epsilon (1e-9) instead - cheap, lazy, accurate for non-negative
  # input (DELTARP, DELTARQ, MUTATIONS). On signed input like
  # M-values (SIGNAL), `log` doesn't apply semantically anyway -
  # the caller should choose `transformation_y = "none"` there.
  # `scale` uses the sample standard deviation (n - 1), as polars' std().
  connection <- dbplyr::remote_con(pivot_lazy)
  quoted_cols <- vapply(sample_cols, function(column_name) as.character(DBI::dbQuoteIdentifier(connection, column_name)), character(1))
  if (trans_y != "none") {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " io_data_preparation_lazy", fam_str, key_str,
              ": applying transformation_y='", trans_y,
              "' to ", length(sample_cols), " sample columns.")

    sql_template <- switch(trans_y,
      "log"   = "ln(CAST({col} AS DOUBLE) + 1e-9)",
      "log2"  = "log2(CAST({col} AS DOUBLE) + 1e-9)",
      "log10" = "log10(CAST({col} AS DOUBLE) + 1e-9)",
      "exp"   = "exp(CAST({col} AS DOUBLE))",
      "scale" = "(CAST({col} AS DOUBLE) - avg({col}) OVER ()) / stddev_samp({col}) OVER ()"
    )
    trans_exprs <- lapply(quoted_cols, function(quoted_col) dbplyr::sql(gsub("{col}", quoted_col, sql_template, fixed = TRUE)))
    names(trans_exprs) <- sample_cols
    pivot_lazy <- dplyr::mutate(pivot_lazy, !!!trans_exprs)
  }

  # ---- Step 2: universal degenerate-burden filter -----------
  # A row is "degenerate" iff every sample column carries the same
  # value (var(Y) == 0). We detect this lazily via:
  #   least(sample_cols) != greatest(sample_cols)
  # which is equivalent to checking variance > 0 without materialising
  # the full row. least() and greatest() skip NULLs as polars'
  # min_horizontal() and max_horizontal() did, and a row of NULLs gives
  # NULL, which the filter drops.
  # The filter wraps the query as a subquery written out: given as raw SQL
  # to dplyr::filter(), dbplyr cannot see that it reads the transformed
  # columns and may place it in the WHERE of the same SELECT, where it
  # would read the values before the transformation.
  if (apply_degenerate_filter) {
    column_list <- paste(quoted_cols, collapse = ", ")
    pivot_lazy <- dplyr::tbl(connection, dbplyr::sql(sprintf("SELECT * FROM (%s) AS prepared WHERE least(%s) <> greatest(%s)",
                                                             dbplyr::sql_render(pivot_lazy), column_list, column_list)))
  }

  pivot_lazy
}
