# ---------------------------------------------------------------------------
# io_stream_merge_bed() - lazy long→wide merge of per-sample bed/bedgraph files
#
# AI-027 (2026-06-01). Returns a polars LazyFrame that, when collected,
# produces a wide pivot table:
#
#   CHR, START, END, <sample_id_1>, <sample_id_2>, ..., <sample_id_N>
#
# Two implementations were tried on the ewas_data_hub DELTAS_HYPO marker
# (~4 k bed files, ~280 M long rows after concat):
#
#   v1 (lapply + concat)        : ~209 s
#   v2 (single scan_csv +       : ~267 s
#        include_file_paths +
#        regex extract)
#
# v1 wins because the regex-extract on 280 M rows in v2 is expensive enough
# to offset the gain of a single scan_csv. The per-file scan in v1 still
# runs in parallel via the polars Rust runtime; the R-side lapply only
# constructs ~4 k lazy plan nodes (cheap), not actual I/O.
# ---------------------------------------------------------------------------

#' Stream-merge a list of per-sample bed/bedgraph files into a wide lazy pivot
#'
#' Each input file is scanned lazily; the basename minus the
#' \code{_<marker>_<figure>.<ext>} suffix is used as the column alias for
#' that sample's value column. All per-sample lazy frames are concatenated
#' vertically and then pivoted long→wide on \code{(CHR, START, END)}.
#'
#' @param bed_files Character vector of file paths to \code{.bed},
#'   \code{.bedgraph}, or \code{.bedgraph.gz} files. Each file must follow the
#'   per-sample naming convention produced by \code{\link{sem_analyze_population}}:
#'   \code{<sample_id>_<marker>_<figure>.bedgraph[.gz]}.
#' @param marker,figure Used to strip the trailing \code{_<marker>_<figure>}
#'   from the basename when deriving the \code{sample_id} alias.
#'
#' @return A lazy \code{polars} frame ready for \code{$collect()},
#'   \code{$sink_parquet()} or further chained operations.
#'
#' @keywords internal
#' @noRd
io_stream_merge_bed <- function(bed_files, marker, figure) {

  if (length(bed_files) == 0L)
    stop("io_stream_merge_bed(): bed_files is empty")

  suffix_pattern <- paste0("_", marker, "_", figure, "\\.(bed|bedgraph)(\\.gz)?$")

  # Build one lazy per-sample frame. The R-side lapply only constructs
  # plan nodes; the actual CSV parsing happens in parallel inside the
  # polars Rust runtime on $collect().
  inputs <- lapply(bed_files, function(f) {
    sample_id <- sub(suffix_pattern, "", basename(f))
    polars::pl$scan_csv(
      f,
      separator  = "\t",
      has_header = FALSE,
      skip_rows  = 0
    )$select(
      polars::pl$col("column_1")$alias("CHR"),
      polars::pl$col("column_2")$alias("START"),
      polars::pl$col("column_3")$alias("END"),
      polars::pl$col("column_4")$alias("value")
    )$with_columns(
      polars::pl$col("CHR")$str$replace("^(?i)chr", "")$alias("CHR"),
      polars::pl$lit(sample_id)$alias("sample_id")
    )
  })

  long <- do.call(polars::pl$concat, c(inputs, list(how = "vertical")))

  # Pivot long → wide. Polars Rust runs this as a parallel hash aggregation
  # (Rayon). No Reduce(join) over thousands of frames involved.
  wide <- long$collect()$pivot(
    values = "value",
    index  = c("CHR", "START", "END"),
    on     = "sample_id"
  )

  wide$sort(c("CHR", "START"), descending = FALSE)$lazy()
}
