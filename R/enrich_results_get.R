# Symmetric counterpart to enrich_result_save(): reads a single
# pathway result CSV from disk and applies optional row filters.
#
# Architectural rule (see backlog AI-020): a `*_get` helper is
# read + filter only. Transforms (p-adjust, category enrichment,
# PHENOTYPE flag, annotation join) belong in a separate pipeline
# step, NOT here. enrich_result_save() currently does mix transform
# with write - that violation is tracked for cleanup in AI-020 and
# kept as-is for now to avoid behaviour drift.
#
# The only label-aware read tweak is `dec = "."` for phenolyzer
# output, which writes US-locale CSV while other backends use the
# European default.
#
# @param file_name           path to a pathway result CSV
# @param label               backend tag from
#                            ssEnv$key_enrichment_format[, "label"]
#                            ("WebGestalt", "pathfindR", "phenolyzer",
#                            "ctdR"). Drives the decimal separator.
# @param significance_column optional column name to filter rows on
# @param alpha               numeric: rows kept when
#                            significance_column < alpha
# @param top                 integer: after sorting by significance_column
#                            ascending, keep only the top N rows
# @param required_columns    character vector: rows where any of these
#                            columns is missing from the CSV are
#                            rejected (returns empty data.frame).
# @return data.frame with the (optionally filtered) CSV contents, or
#         an empty data.frame when the file does not exist, the
#         required columns are missing, or every row was filtered out.
# @keywords internal
enrich_results_get <- function(file_name,
                                label = "default",
                                significance_column = NULL,
                                alpha = NULL,
                                top = NULL,
                                required_columns = NULL) {
  if (!file.exists(file_name))
    return(data.frame())

  df <- if (label == "phenolyzer")
    utils::read.csv2(file_name, dec = ".")
  else
    utils::read.csv2(file_name)

  if (!is.null(required_columns) &&
      !all(required_columns %in% colnames(df)))
    return(data.frame())

  if (!is.null(significance_column) && !is.null(alpha) &&
      significance_column %in% colnames(df))
    df <- df[df[, significance_column] < alpha, , drop = FALSE]

  if (!is.null(top) && !is.null(significance_column) &&
      significance_column %in% colnames(df) &&
      nrow(df) > 0) {
    df <- df[order(df[, significance_column]), , drop = FALSE]
    df <- utils::head(df, as.integer(top))
  }

  df
}
