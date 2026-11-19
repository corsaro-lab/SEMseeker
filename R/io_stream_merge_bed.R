# ---------------------------------------------------------------------------
# io_stream_merge_bed() - long-to-wide merge of per-sample bed/bedgraph files
#
# Writes a wide pivot table to parquet:
#
#   CHR, START, END, <sample_id_1>, <sample_id_2>, ..., <sample_id_N>
#
# It is the heaviest operation of the per-sample engine: on the largest study
# run so far about 4 000 files and 280 million long rows. The data engine
# reads every file in one read_csv, pivots and copies the result to parquet
# without bringing it into R. Measured on 1 000 synthetic files of 70 000 rows
# (70 million long rows) before the move from polars: 28.6 s against polars'
# 79.3 s with no memory limit, and 72 s with an 8 GB limit at 13.5 GB
# resident, where polars took 30.4 GB and offers no limit.
#
# The sample columns follow the order of the files, as the polars merge did,
# and the rows are sorted on CHR and START: the pivot is the same file
# whichever engine wrote it.
# ---------------------------------------------------------------------------

#' Merge per-sample bed/bedgraph files into a wide parquet pivot
#'
#' Each file contributes one sample column, named after the file's basename
#' minus the \code{_<marker>_<figure>.<ext>} suffix; the rows are the union of
#' the positions \code{(CHR, START, END)} of all files, with \code{NULL} where
#' a sample has no value. A leading \code{chr} is removed from CHR.
#'
#' A position that appears twice in the same sample is refused: the wide
#' table has one cell per position and sample, and keeping either value would
#' be a silent choice.
#'
#' @param bed_files Character vector of file paths to \code{.bed},
#'   \code{.bedgraph}, or \code{.bedgraph.gz} files, tab separated, without a
#'   header, with the columns chromosome, start, end and value. Each file must
#'   follow the per-sample naming convention produced by
#'   \code{\link{sem_analyze_population}}:
#'   \code{<sample_id>_<marker>_<figure>.bedgraph[.gz]}.
#' @param marker,figure Used to strip the trailing \code{_<marker>_<figure>}
#'   from the basename when deriving the sample column name.
#' @param parquet_path The parquet file to write; an existing file is
#'   replaced.
#'
#' @return \code{parquet_path}, invisibly.
#'
#' @keywords internal
#' @noRd
io_stream_merge_bed <- function(bed_files, marker, figure, parquet_path) {

  if (length(bed_files) == 0L)
    stop("io_stream_merge_bed(): bed_files is empty")

  connection <- io_engine_connection()
  suffix_pattern <- paste0("_", marker, "_", figure, "\\.(bed|bedgraph)(\\.gz)?$")
  # The file-to-sample mapping is computed in R, with the same expression the
  # polars merge used, and joined on the file name DuckDB reports.
  file_samples <- data.frame(filename  = normalizePath(bed_files, winslash = "/", mustWork = TRUE),
                             sample_id = sub(suffix_pattern, "", basename(bed_files)),
                             stringsAsFactors = FALSE)
  mapping_view <- paste0("merge_files_", gsub("[^A-Za-z0-9]", "", basename(tempfile(""))))
  duckdb::duckdb_register(connection, mapping_view, file_samples)
  long_view <- paste0(mapping_view, "_long")
  on.exit({
    DBI::dbExecute(connection, paste("DROP VIEW IF EXISTS", long_view))
    duckdb::duckdb_unregister(connection, mapping_view)
  }, add = TRUE)

  file_list <- paste(DBI::dbQuoteString(connection, file_samples$filename), collapse = ", ")
  DBI::dbExecute(connection, sprintf("
    CREATE TEMP VIEW %s AS
    SELECT regexp_replace(bed.CHR, '^(?i)chr', '') AS CHR, bed.START, bed.\"END\", bed.value, files.sample_id
    FROM read_csv([%s], delim = '\t', header = false, filename = true, auto_detect = false,
                  columns = {'CHR': 'VARCHAR', 'START': 'BIGINT', 'END': 'BIGINT', 'value': 'DOUBLE'}) AS bed
    JOIN %s AS files ON bed.filename = files.filename",
    long_view, file_list, mapping_view))

  duplicates <- DBI::dbGetQuery(connection, sprintf("
    SELECT sample_id, CHR, START, \"END\", count(*) AS copies FROM %s
    GROUP BY sample_id, CHR, START, \"END\" HAVING count(*) > 1 LIMIT 5", long_view))
  if (nrow(duplicates) > 0L)
    stop("io_stream_merge_bed(): position(s) repeated within one sample, for instance ",
         paste(sprintf("%s at %s:%s-%s", duplicates$sample_id, duplicates$CHR, duplicates$START, duplicates$END), collapse = "; "),
         ". Each bed file must hold one value per position.", call. = FALSE)

  # IN (...) fixes the sample columns to the order of the files; without it
  # PIVOT would order them alphabetically.
  sample_list <- paste(DBI::dbQuoteString(connection, unique(file_samples$sample_id)), collapse = ", ")
  DBI::dbExecute(connection, sprintf("
    COPY (
      SELECT * FROM (PIVOT %s ON sample_id IN (%s) USING first(value) GROUP BY CHR, START, \"END\")
      ORDER BY CHR, START
    ) TO %s (FORMAT parquet, COMPRESSION zstd)",
    long_view, sample_list, DBI::dbQuoteString(connection, parquet_path)))

  invisible(parquet_path)
}
