
#' Convert a POSITION-keyed signal pivot into a PROBE-keyed lazy frame.
#'
#' Returns a Polars LazyFrame rather than an R
#' data.frame. Materialisation downstream is opt-in only - sem_analyze_batch
#' and sem_analyze_population_bulk consume the lazy frame directly. This
#' eliminates the ~12 GB R-side peak that caused silent jetsam kills on
#' ewas-scale matrices (367k × 4013) in resume mode (v18, v21).
#'
#' Contract change: callers that need PROBE as a column see it as a
#' regular polars String column. Callers that previously read `rownames`
#' must switch to selecting / collecting the `PROBE` column.
#'
#' Order: input must already be sorted (sort gate = io_signal_save chunked
#' per-chr by START). This function preserves order via `unique(maintain_order=TRUE)`.
#'
#' @param signal_data polars LazyFrame, polars DataFrame, or R data.frame
#'   carrying CHR/START/END + sample columns. Lazy input keeps the entire
#'   pipeline in zero-copy mode; R data.frame is wrapped before the join.
#'
#' @return A polars LazyFrame with the PROBE column and one column per
#'   sample. CHR/START/END/tech annotation columns are dropped after the
#'   join, since downstream consumers only need the probe-keyed shape.
anno_position_pivot_to_probe <- function(signal_data)
{
 
  ssEnv          <- core_get_session_info()
  tech_col       <- ssEnv$tech
  probe_features <- anno_probe_annotation_build(ssEnv$tech)[, c("PROBE", "CHR", "START", "END", tech_col)]
  probe_features_lf <- polars::as_polars_df(probe_features)$lazy()

  # Accept polars LazyFrame, polars DataFrame, or R data.frame/matrix.
  # Staying lazy here avoids materialising the full SIGNAL pivot (e.g. 367k
  # probes × 4k samples ≈ 12 GB) into R memory before the join.
  signal_lazy <-
    if (inherits(signal_data, "polars_lazy_frame"))      signal_data
    else if (inherits(signal_data, "polars_data_frame")) signal_data$lazy()
    else                                                  polars::as_polars_df(signal_data)$lazy()

  # Single lazy chain: join on coords, dedupe preserving input order
  # (sort gate is io_signal_save per AI-096 §single-sort-gate-at-pivot-save),
  # filter probes matching the detected technology, drop annotation columns.
  # NO collect() - caller takes a LazyFrame. Materialisation, if needed,
  # is the caller's explicit decision.
  result_lf <- probe_features_lf$
    join(signal_lazy, on = c("CHR", "START", "END"), how = "inner")$
    unique(maintain_order = TRUE)$
    filter(polars::pl$col(tech_col)$is_not_null() & polars::pl$col(tech_col))$
    drop(c("CHR", "START", "END", tech_col))

  result_lf
}
