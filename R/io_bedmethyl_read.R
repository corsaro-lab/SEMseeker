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
#'   \item \code{modified_base_code} (\code{m} for 5mC, \code{h} for 5hmC)
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
#' Only the 5mC rows (code \code{m}) are read. \code{modkit pileup} writes
#' one row per modification at the same position, 5mC and 5hmC, and SEMseeker
#' measures methylation: keeping both would give a position two values in the
#' same sample. A file with no 5mC row is refused, because its sample would
#' enter the analysis as a column with no value. A position that still
#' appears twice in one sample after the filter is refused too, because the
#' wide table has one cell per position and sample and keeping either value
#' would be a silent choice.
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

  if (anyDuplicated(sample_ids))
    stop("io_bedmethyl_read(): sample ids must be unique; repeated: ",
         paste(unique(sample_ids[duplicated(sample_ids)]), collapse = ", "), call. = FALSE)

  connection <- io_engine_connection()
  file_samples <- data.frame(filename  = normalizePath(file_paths, winslash = "/", mustWork = TRUE),
                             sample_id = sample_ids,
                             stringsAsFactors = FALSE)
  mapping_view <- paste0("bedmethyl_files_", gsub("[^A-Za-z0-9]", "", basename(tempfile(""))))
  duckdb::duckdb_register(connection, mapping_view, file_samples)
  long_view <- paste0(mapping_view, "_long")
  on.exit({
    DBI::dbExecute(connection, paste("DROP VIEW IF EXISTS", long_view))
    duckdb::duckdb_unregister(connection, mapping_view)
  }, add = TRUE)

  # bedmethyl uses positional columns: the first eleven are named, the
  # remaining modkit columns keep DuckDB's names and are not read. union_by_name
  # lets files with a different number of trailing columns be read together.
  file_list <- paste(DBI::dbQuoteString(connection, file_samples$filename), collapse = ", ")
  DBI::dbExecute(connection, sprintf("
    CREATE TEMP VIEW %s AS
    SELECT bed.CHR, bed.START, bed.\"END\", bed.modification, bed.N_valid_cov,
           bed.fraction_modified / 100 AS value, files.sample_id
    FROM read_csv([%s], delim = '\t', header = false, filename = true, union_by_name = true,
                  names = ['CHR', 'START', 'END', 'modification', 'score', 'strand', 'start_code',
                           'end_code', 'color', 'N_valid_cov', 'fraction_modified'],
                  types = {'CHR': 'VARCHAR', 'START': 'BIGINT', 'END': 'BIGINT', 'modification': 'VARCHAR',
                           'N_valid_cov': 'BIGINT', 'fraction_modified': 'DOUBLE'}) AS bed
    JOIN %s AS files ON bed.filename = files.filename",
    long_view, file_list, mapping_view))

  without_methylation <- DBI::dbGetQuery(connection, sprintf("
    SELECT files.sample_id FROM %s AS files
    WHERE NOT EXISTS (SELECT 1 FROM %s AS bed WHERE bed.sample_id = files.sample_id AND bed.modification = 'm')",
    mapping_view, long_view))$sample_id
  if (length(without_methylation) > 0L)
    stop("io_bedmethyl_read(): no 5mC row (modification code 'm') for sample(s) ",
         paste(without_methylation, collapse = ", "), ".", call. = FALSE)

  duplicates <- DBI::dbGetQuery(connection, sprintf("
    SELECT sample_id, CHR, START, \"END\" FROM %s WHERE modification = 'm'
    GROUP BY sample_id, CHR, START, \"END\" HAVING count(*) > 1 LIMIT 5", long_view))
  if (nrow(duplicates) > 0L)
    stop("io_bedmethyl_read(): 5mC position(s) repeated within one sample, for instance ",
         paste(sprintf("%s at %s:%s-%s", duplicates$sample_id, duplicates$CHR, duplicates$START, duplicates$END), collapse = "; "),
         ". Each file must hold one 5mC value per position.", call. = FALSE)

  # The coverage filter drops the row before the pivot, so a position below
  # the threshold in one sample is NA there and kept for the others.
  sample_list <- paste(DBI::dbQuoteString(connection, sample_ids), collapse = ", ")
  result <- DBI::dbGetQuery(connection, sprintf("
    PIVOT (SELECT CHR, START, \"END\", value, sample_id FROM %s
           WHERE modification = 'm' AND N_valid_cov >= %d)
    ON sample_id IN (%s) USING first(value) GROUP BY CHR, START, \"END\"",
    long_view, as.integer(min_coverage), sample_list))

  # Stable ordering (chromosome-aware sort delegated to downstream pipeline)
  result[order(result$CHR, result$START), , drop = FALSE]
}
