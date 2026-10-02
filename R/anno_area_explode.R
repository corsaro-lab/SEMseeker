#' Split composite area names into one row per area (internal)
#'
#' Extracted by AI-255 from `anno_annotate_position_pivots()` so the on-demand
#' builder ([io_pivot_build()]) and the batch annotation share one implementation
#' - two copies of this logic would drift, and the drift would be silent.
#'
#' AI-050: Bioconductor annotation packages assign some probes to multiple genes
#' (intergenic overlaps, antisense, ...), producing composite `AREA` strings like
#' `"NUDT6;SPATA5"`. Treating the composite as a single gene was a regression
#' that made `assoc_apply_stat_model()` fail to parse (PVALUE=NA) and inflated
#' false positives downstream, because one p-value got smeared across N
#' enrichment hits. Splitting gives one clean mono-gene row per area.
#'
#' AI-061+: `","` and `"/"` are separators too. Comma-separated tokens are
#' already complete HGNC symbols; `"/"` needs the smart split of
#' [.anno_smart_split_area_name()] so `"HLA-A/B/C"` becomes
#' `("HLA-A","HLA-B","HLA-C")` rather than `("HLA-A","B","C")`.
#'
#' **This is correct, and it is also why an exploded frame must never be summed
#' across rows**: a probe on three genes contributes to three rows because those
#' are three different questions. Adding them counts it three times - see
#' [io_pivot_build()], where the `SAMPLE` scope masks positions instead.
#'
#' @param lazy a polars LazyFrame carrying an `AREA` column.
#' @return the LazyFrame with `AREA` exploded and trimmed.
#' @keywords internal
#' @noRd
.anno_area_explode <- function(lazy) {

  slashed_areas <- as.data.frame(
    lazy$select("AREA")$
      filter(polars::pl$col("AREA")$str$contains("/", literal = TRUE))$
      unique()$
      collect()
  )$AREA

  if (length(slashed_areas) > 0L) {
    expansions <- vapply(
      slashed_areas,
      function(s) paste(.anno_smart_split_area_name(s), collapse = ";"),
      character(1)
    )
    mapping_lf <- polars::as_polars_df(data.frame(
      AREA     = slashed_areas,
      AREA_NEW = expansions,
      stringsAsFactors = FALSE
    ))$lazy()
    lazy <- lazy$join(mapping_lf, on = "AREA", how = "left")$with_columns(
      polars::pl$when(polars::pl$col("AREA_NEW")$is_not_null())$
        then(polars::pl$col("AREA_NEW"))$
        otherwise(polars::pl$col("AREA"))$alias("AREA")
    )$drop("AREA_NEW")
  }

  lazy <- lazy$with_columns(
    polars::pl$col("AREA")$str$replace_all(",", ";")$str$split(";")
  )$explode("AREA")

  lazy <- lazy$with_columns(polars::pl$col("AREA")$str$strip_chars())

  # Drop rows where AREA became empty after split/strip (e.g. trailing
  # semicolons from malformed annotations).
  lazy$filter(polars::pl$col("AREA")$str$len_chars() > 0)
}
