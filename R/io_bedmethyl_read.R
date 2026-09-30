#' Read bedmethyl files (modkit/nanopolish) into a SEMseeker coordinate data frame
#'
#' Parses one or more bedmethyl files (produced by \code{modkit pileup} or
#' \code{nanopolish call-methylation}) and returns a wide data frame in
#' SEMseeker coordinate format: \code{CHR, START, END, sample1, sample2, ...}.
#' Values are methylation fractions in [0, 1]. Positions with coverage below
#' \code{min_coverage} are set to \code{NA}.
#'
#' Expected bedmethyl schema (modkit, tab-separated, no header):
#' \enumerate{
#'   \item \code{chrom}
#'   \item \code{start_position}
#'   \item \code{end_position}
#'   \item \code{modified_base_code} (e.g. \code{m} for 5mC)
#'   \item \code{score}
#'   \item \code{strand}
#'   \item \code{start_code}
#'   \item \code{end_code}
#'   \item \code{color}
#'   \item \code{N_valid_cov}
#'   \item \code{fraction_modified} (percent, 0-100)
#'   \item ... remaining modkit columns (ignored)
#' }
#'
#' @param file_paths Character vector of paths to bedmethyl files.
#' @param sample_ids Optional character vector of sample IDs, one per file.
#'   If \code{NULL} (default), IDs are derived from the file basenames
#'   (extension stripped).
#' @param min_coverage Integer. Positions with \code{N_valid_cov < min_coverage}
#'   are dropped. Default 5 (common threshold for Nanopore methylation calls).
#'
#' @return A data frame with columns \code{CHR, START, END} followed by one
#'   numeric column per sample (values in \code{[0, 1]}, \code{NA} if dropped
#'   by coverage filter or absent in that sample).
#' @keywords internal
io_bedmethyl_read <- function(file_paths, sample_ids = NULL, min_coverage = 5L) {
  if (!length(file_paths))
    stop("io_bedmethyl_read(): no file_paths provided.")

  missing_files <- file_paths[!file.exists(file_paths)]
  if (length(missing_files))
    stop("io_bedmethyl_read(): file(s) not found: ",
         paste(missing_files, collapse = ", "))

  if (is.null(sample_ids)) {
    sample_ids <- tools::file_path_sans_ext(basename(file_paths))
  } else if (length(sample_ids) != length(file_paths)) {
    stop("io_bedmethyl_read(): sample_ids length must match file_paths length.")
  }

  # bedmethyl uses positional columns - name only the ones we keep
  bedmethyl_cols <- c("CHR", "START", "END", "mod", "score", "strand",
                      "start_code", "end_code", "color",
                      "N_valid_cov", "fraction_modified")

  per_file <- vector("list", length(file_paths))
  for (i in seq_along(file_paths)) {
    lf <- polars::pl$scan_csv(
      file_paths[i],
      separator = "\t",
      has_header = FALSE,
      n_rows = NULL
    )
    # select & rename first 11 columns (modkit layout)
    df <- lf$select(
      polars::pl$nth(0L)$alias("CHR"),
      polars::pl$nth(1L)$alias("START"),
      polars::pl$nth(2L)$alias("END"),
      polars::pl$nth(9L)$alias("N_valid_cov"),
      polars::pl$nth(10L)$alias("fraction_modified")
    )$filter(polars::pl$col("N_valid_cov") >= min_coverage
    )$with_columns(
      (polars::pl$col("fraction_modified") / 100)$alias(sample_ids[i])
    )$select(
      polars::pl$col("CHR"),
      polars::pl$col("START"),
      polars::pl$col("END"),
      polars::pl$col(sample_ids[i])
    )$collect()
    per_file[[i]] <- as.data.frame(df)
  }

  # Outer-join all samples on (CHR, START, END)
  result <- per_file[[1L]]
  if (length(per_file) > 1L) {
    for (i in seq.int(2L, length(per_file))) {
      result <- merge(result, per_file[[i]],
                      by = c("CHR", "START", "END"), all = TRUE)
    }
  }

  # Stable ordering (chromosome-aware sort delegated to downstream pipeline)
  result[order(result$CHR, result$START), , drop = FALSE]
}
