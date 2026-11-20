#' Split composite area names into one row per area (internal)
#'
#' Extracted from `anno_annotate_position_pivots()` by the taxonomy unification so the on-demand
#' builder ([io_pivot_build()]) and the batch annotation share one implementation
#' - two copies of this logic would drift, and the drift would be silent.
#'
#' Bioconductor annotation packages assign some probes to multiple genes
#' (intergenic overlaps, antisense, ...), producing composite `AREA` strings like
#' `"NUDT6;SPATA5"`. Treating the composite as a single gene was a regression
#' that made `assoc_apply_stat_model()` fail to parse (PVALUE=NA) and inflated
#' false positives downstream, because one p-value got smeared across N
#' enrichment hits. Splitting gives one clean mono-gene row per area.
#'
#' `","` and `"/"` are separators too. Comma-separated tokens are
#' already complete HGNC symbols; `"/"` needs the smart split of
#' [.anno_smart_split_area_name()] so `"HLA-A/B/C"` becomes
#' `("HLA-A","HLA-B","HLA-C")` rather than `("HLA-A","B","C")`.
#'
#' **This is correct, and it is also why an exploded frame must never be summed
#' across rows**: a probe on three genes contributes to three rows because those
#' are three different questions. Adding them counts it three times - see
#' [io_pivot_build()], where the `SAMPLE` scope masks positions instead.
#'
#' @param connection the data engine connection the query runs on.
#' @param query_text SQL of a query with an `AREA` column.
#' @return SQL of the same query with `AREA` exploded and trimmed, the other
#'   columns unchanged and in their order.
#' @keywords internal
#' @noRd
.anno_area_explode <- function(connection, query_text) {

  slashed_areas <- DBI::dbGetQuery(connection, sprintf(
    "SELECT DISTINCT AREA FROM (%s) AS areas WHERE contains(AREA, '/')", query_text))$AREA

  if (length(slashed_areas) > 0L) {
    expansions <- vapply(
      slashed_areas,
      function(s) paste(.anno_smart_split_area_name(s), collapse = ";"),
      character(1)
    )
    mapping_query <- dbplyr::sql_render(io_register_frame(data.frame(
      AREA     = slashed_areas,
      AREA_NEW = unname(expansions),
      stringsAsFactors = FALSE
    )))
    query_text <- sprintf(
      "SELECT source.* REPLACE (COALESCE(mapping.AREA_NEW, source.AREA) AS AREA)\nFROM (%s) AS source\nLEFT JOIN (%s) AS mapping ON source.AREA = mapping.AREA",
      query_text, mapping_query)
  }

  # Split on ";" (and "," turned into ";"), one row per token, surrounding
  # whitespace removed. Rows whose AREA became empty after split/strip (e.g.
  # trailing semicolons from malformed annotations) are dropped.
  exploded <- sprintf(
    "SELECT * REPLACE (unnest(string_split(replace(AREA, ',', ';'), ';')) AS AREA) FROM (%s) AS composite", query_text)
  sprintf(
    "SELECT * FROM (SELECT * REPLACE (regexp_replace(AREA, '^\\s+|\\s+$', '', 'g') AS AREA) FROM (%s) AS split) AS trimmed WHERE length(AREA) > 0",
    exploded)
}
