# util_join_values_to_thresholds
#
# Private helper - inner-join a per-sample values data.frame to a thresholds
# data.frame on (CHR, START, END) using Polars.
#
# This is the canonical join used by sem_mutations_get(), sem_delta_single_sample(),
# and sem_deltar_single_sample() - all three need the same intersection.
#
# WHY POLARS:
#   Nanopore bedmethyl files can have 28M+ rows.  Base-R merge()/match()
#   is too slow at that scale.  Polars handles it with ease via lazy evaluation.
#
# WHY INNER JOIN (not positional zip):
#   The old sort-then-zip assumed values and thresholds had identical position
#   sets. Cross-run analysis (e.g. Nanopore sample vs Illumina reference batch
#   passed via populationControlRangeBetaValues) breaks that assumption silently.
#   An inner join is correct regardless of overlap.
#
# @param values     data.frame: CHR, START, END, VALUE (col 4)
# @param thresholds data.frame: CHR, START, END + threshold columns
#                   (signal_inferior_thresholds, signal_superior_thresholds,
#                    signal_median_values, iqr, q1, q3 - any subset is fine)
# @return data.frame - inner join result; columns CHR, START, END, VALUE plus
#         whichever threshold columns were present in the input.
#         Returns 0-row data.frame if there is no positional overlap.
util_join_values_to_thresholds <- function(values, thresholds) {

  # E-13: Normalise CHR - strip "chr" prefix so bed-file values (chr1) match
  # threshold values (1). io_dump_sample_as_bed_file() prepends "chr" when writing
  # bed files, but signal_thresholds retains bare chromosome numbers from the
  # probe annotation. Without this normalisation the inner join returns 0 rows.
  strip_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)

  vals <- data.frame(
    CHR   = strip_chr(values$CHR),
    START = as.integer(values$START),
    END   = as.integer(values$END),
    VALUE = values[, 4L],
    stringsAsFactors = FALSE
  )

  # Select only positional key + known threshold columns - avoids carrying
  # large probe-annotation columns (GENE_*, ISLAND_*, ...) through the join.
  keep_cols <- intersect(
    c("CHR", "START", "END",
      "signal_inferior_thresholds", "signal_superior_thresholds",
      "signal_median_values", "iqr", "q1", "q3"),
    colnames(thresholds)
  )
  thr <- thresholds[, keep_cols, drop = FALSE]
  thr$CHR   <- strip_chr(thr$CHR)
  thr$START <- as.integer(thr$START)
  thr$END   <- as.integer(thr$END)

  vals_lf <- polars::as_polars_df(vals)$lazy()
  thr_lf  <- polars::as_polars_df(thr)$lazy()

  as.data.frame(
    vals_lf$join(thr_lf, on = c("CHR", "START", "END"), how = "inner")$collect()
  )
}
