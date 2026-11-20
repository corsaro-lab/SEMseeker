# util_join_values_to_thresholds
#
# Private helper - inner-join a per-sample values data.frame to a thresholds
# data.frame on (CHR, START, END) in the data engine.
#
# This is the canonical join used by sem_mutations_get(), sem_delta_single_sample(),
# and sem_deltar_single_sample() - all three need the same intersection.
#
# WHY THE DATA ENGINE:
#   Nanopore bedmethyl files can have 28M+ rows.  Base-R merge()/match()
#   is too slow at that scale.  The join runs in DuckDB, and the rows keep
#   the order of `values`, as the polars join did.
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

  vals$VALUES_ORDER <- seq_len(nrow(vals))
  connection <- io_engine_connection()
  vals_view <- paste0("values_", gsub("[^A-Za-z0-9]", "", basename(tempfile(""))))
  thr_view  <- paste0(vals_view, "_thresholds")
  duckdb::duckdb_register(connection, vals_view, vals)
  duckdb::duckdb_register(connection, thr_view, thr)
  on.exit({
    duckdb::duckdb_unregister(connection, vals_view)
    duckdb::duckdb_unregister(connection, thr_view)
  }, add = TRUE)

  joined <- DBI::dbGetQuery(connection, sprintf('SELECT * FROM %s JOIN %s USING (CHR, START, "END") ORDER BY VALUES_ORDER',
                                                vals_view, thr_view))
  joined$VALUES_ORDER <- NULL
  joined
}
