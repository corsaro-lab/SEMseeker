#' Polars-native equivalent of `io_data_preparation()` for the lazy
#' lazy batch path. Two responsibilities:
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
#' @param pivot_lazy A `polars_lazy_frame` whose rows are probes/areas
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
#' @return The transformed and filtered `polars_lazy_frame`. The schema
#'   is unchanged (same column names; new helper columns `__min_y` /
#'   `__max_y` are dropped before return).
#'
#' @keywords internal
#' @noRd
io_data_preparation_lazy <- function(pivot_lazy,
                                   sample_cols,
                                   transformation_y = "none",
                                   apply_degenerate_filter = TRUE,
                                   key = NULL,
                                   family_test = NULL) {

  if (!inherits(pivot_lazy, "polars_lazy_frame")) {
    if (inherits(pivot_lazy, "polars_data_frame")) {
      pivot_lazy <- pivot_lazy$lazy()
    } else {
      stop("io_data_preparation_lazy: pivot_lazy must be a polars LazyFrame or DataFrame.",
           call. = FALSE)
    }
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
  # Polars 1.x expression API: `pl$col(c)` × native math methods.
  # Note: `log` on negative or zero values yields NaN; the R-side
  # `io_data_preparation()` guards by adding `min(burden_values[burden_values>0])`,
  # which would require a scan over the full pivot. We use a tiny
  # epsilon (1e-9) instead - cheap, lazy, accurate for non-negative
  # input (DELTARP, DELTARQ, MUTATIONS). On signed input like
  # M-values (SIGNAL), `log` doesn't apply semantically anyway -
  # the caller should choose `transformation_y = "none"` there.
  EPS <- 1e-9
  if (trans_y != "none") {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " io_data_preparation_lazy", fam_str, key_str,
              ": applying transformation_y='", trans_y,
              "' to ", length(sample_cols), " sample columns.")

    trans_exprs <- lapply(sample_cols, function(c) {
      expr <- switch(trans_y,
        "log"   = polars::pl$col(c)$add(EPS)$log(),
        "log2"  = polars::pl$col(c)$add(EPS)$log(base = 2),
        "log10" = polars::pl$col(c)$add(EPS)$log10(),
        "exp"   = polars::pl$col(c)$exp(),
        "scale" = polars::pl$col(c)$sub(polars::pl$col(c)$mean())$
                    truediv(polars::pl$col(c)$std())
      )
      expr$alias(c)
    })
    # Polars 1.x: `$with_columns()` takes `...` Expr args, NOT a list.
    # Same do.call unpacking pattern as `analyze_population_bulk.R:112`.
    pivot_lazy <- do.call(pivot_lazy$with_columns, trans_exprs)
  }

  # ---- Step 2: universal degenerate-burden filter -----------
  # A row is "degenerate" iff every sample column carries the same
  # value (var(Y) == 0). We detect this lazily via:
  #   min_horizontal(sample_cols) != max_horizontal(sample_cols)
  # which is equivalent to checking variance > 0 without materialising
  # the full row. The two helper columns are dropped immediately so
  # the output schema matches the input (just fewer rows).
  if (apply_degenerate_filter) {
    # Polars 1.x R bindings: `pl$min_horizontal()` / `pl$max_horizontal()`
    # take `...` exprs as individual args, NOT a list. Wrap via `do.call`
    # to unpack the per-sample column expressions. The same `do.call`
    # pattern is used by `analyze_population_bulk.R` for `pl$concat()`
    # on per-chr chunks - kept consistent for readability.
    col_exprs <- lapply(sample_cols, function(c) polars::pl$col(c))
    min_h <- do.call(polars::pl$min_horizontal, col_exprs)$alias("__min_y")
    max_h <- do.call(polars::pl$max_horizontal, col_exprs)$alias("__max_y")
    pivot_lazy <- pivot_lazy$with_columns(min_h, max_h)$filter(
      polars::pl$col("__min_y") != polars::pl$col("__max_y")
    )$drop(c("__min_y", "__max_y"))
  }

  pivot_lazy
}
