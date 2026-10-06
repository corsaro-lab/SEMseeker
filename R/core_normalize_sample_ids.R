#' Normalise sample identifiers on both sides of the pipeline (internal)
#'
#' Single point of truth for sample-identifier normalisation. Every entry point
#' that receives a user-supplied sample sheet together with a signal matrix must
#' route both through this helper BEFORE any column subsetting.
#'
#' Background: the pipeline cleans `colnames(signal_data)` with
#' [core_name_cleaning()] (uppercase, non-alphanumeric -> "_") but historically
#' compared/subset those columns with the RAW `Sample_ID` values coming from the
#' sample sheet. With identifiers containing "-", "." or spaces the match is
#' empty and R fails with the opaque `undefined columns selected`, or - worse -
#' silently produces an empty result (`%in%` filters). The bug stayed invisible
#' because the reference datasets use `GSM\\d+` identifiers, on which
#' [core_name_cleaning()] is a no-op.
#'
#' The transformation is idempotent: applying it twice yields the same result,
#' so calling this helper again downstream is always safe.
#'
#' @param sample_sheet Data frame with a sample-identifier column.
#' @param signal_data Optional matrix / data frame whose columns are samples
#'   (probe identifiers live in the rownames). Structural columns are left
#'   untouched.
#' @param sample_id_column Name of the identifier column (default
#'   `"Sample_ID"`).
#' @param require_all_ids If `TRUE`, every sample-sheet identifier must have a
#'   matching column in `signal_data`, otherwise the call fails with an explicit
#'   message listing the offenders. Use it right before any positional/name
#'   subsetting of `signal_data`.
#' @param structural_columns Column names that are NOT samples and must be
#'   preserved as-is.
#'
#' @return A list with the normalised `sample_sheet` and `signal_data`, plus
#'   `unmatched_ids` (sample-sheet identifiers with no signal column) and
#'   `unmatched_columns` (signal columns absent from the sample sheet).
#'
#' @keywords internal
#' @noRd
core_normalize_sample_ids <- function(sample_sheet,
                                      signal_data = NULL,
                                      sample_id_column = "Sample_ID",
                                      require_all_ids = FALSE,
                                      structural_columns = c("PROBE", "AREA", "CHR", "START", "END")) {

  if (!is.null(sample_sheet)) {
    if (!is.data.frame(sample_sheet))
      sample_sheet <- as.data.frame(sample_sheet, stringsAsFactors = FALSE)
    if (!(sample_id_column %in% colnames(sample_sheet)))
      stop("core_normalize_sample_ids(): the sample sheet has no '",
           sample_id_column, "' column. Available columns: ",
           paste(colnames(sample_sheet), collapse = ", "))

    raw_ids <- as.character(sample_sheet[, sample_id_column])
    clean_ids <- core_name_cleaning(raw_ids)
    .core_stop_on_id_collision(raw_ids, clean_ids, "sample sheet")
    sample_sheet[, sample_id_column] <- clean_ids
  } else {
    clean_ids <- character(0)
  }

  if (!is.null(signal_data) && !is.null(colnames(signal_data))) {
    raw_cols <- colnames(signal_data)
    is_structural <- raw_cols %in% structural_columns
    clean_cols <- raw_cols
    clean_cols[!is_structural] <- core_name_cleaning(raw_cols[!is_structural])
    .core_stop_on_id_collision(raw_cols[!is_structural],
                               clean_cols[!is_structural],
                               "signal data columns")
    colnames(signal_data) <- clean_cols
  }

  sample_columns <- if (is.null(signal_data)) character(0) else
    setdiff(colnames(signal_data), structural_columns)

  unmatched_ids <- setdiff(unique(clean_ids), sample_columns)
  unmatched_columns <- setdiff(sample_columns, unique(clean_ids))

  if (isTRUE(require_all_ids) && !is.null(signal_data) && length(unmatched_ids) > 0)
    stop("core_normalize_sample_ids(): ", length(unmatched_ids),
         " sample sheet identifier(s) have no matching column in the signal data ",
         "even after normalisation: ",
         paste(utils::head(unmatched_ids, 10), collapse = ", "),
         if (length(unmatched_ids) > 10) ", ..." else "",
         ". First signal data columns: ",
         paste(utils::head(sample_columns, 4), collapse = ", "), ".")

  list(sample_sheet = sample_sheet,
       signal_data = signal_data,
       unmatched_ids = unmatched_ids,
       unmatched_columns = unmatched_columns)
}

#' Fail when normalisation collapses distinct identifiers onto the same name
#'
#' Forcing [core_name_cleaning()] on user-supplied identifiers can merge two
#' distinct samples (e.g. `"S-1"` and `"S.1"` both become `"S_1"`). That would
#' silently mix samples, so it is always an error.
#'
#' @keywords internal
#' @noRd
.core_stop_on_id_collision <- function(raw, clean, origin) {
  if (length(raw) == 0) return(invisible(NULL))
  keep <- !is.na(raw) & !is.na(clean)
  raw <- raw[keep]; clean <- clean[keep]
  if (length(raw) == 0) return(invisible(NULL))

  n_distinct <- tapply(raw, clean, function(x) length(unique(x)))
  colliding <- names(n_distinct)[n_distinct > 1]
  if (length(colliding) == 0) return(invisible(NULL))

  details <- vapply(colliding, function(cl) {
    paste0(cl, " <- {", paste(sort(unique(raw[clean == cl])), collapse = ", "), "}")
  }, character(1))

  stop("core_normalize_sample_ids(): name normalisation collapses distinct ",
       origin, " identifiers onto the same name: ",
       paste(utils::head(details, 5), collapse = "; "),
       if (length(details) > 5) ", ..." else "",
       ". Rename them upstream so they stay distinct after ",
       "uppercase + non-alphanumeric replacement.")
}
