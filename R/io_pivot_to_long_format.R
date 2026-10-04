#' A pivot in long format: one row per instance and sample
#'
#' Reads the pivot of one class and returns it melted, so that a plotting
#' function gets one row per (instance, sample) pair with the phenotype of that
#' sample already attached. Wide pivots carry one column per sample, which is
#' the shape the pipeline writes and the wrong shape for every grammar-of-graphics
#' layer.
#'
#' @section What this used to do, and why none of it worked:
#' Through 0.99.5 this function could not return anything. It called
#' `io_read_pivot()`, which returns a **lazy polars frame**, and then used the
#' result as if it were a data.frame: `subset(area_pivot, area_pivot$SAMPLEID ==
#' ...)` on a lazy frame raises `SAMPLEID is not a member of this polars object`,
#' and `area_pivot[-1, ]` raises `Cannot subset rows of a LazyFrame`. The second
#' line of the body was already an error, so no caller of this function had ever
#' drawn anything from a real pivot.
#'
#' Four more defects were behind that one, and they matter because fixing only
#' the first would have produced an empty answer instead of an error:
#'
#' * it filtered on a column called `SAMPLEID`. Pivots have `AREA`. The name
#'   came from the CSV pivots, which carried a first row of sample groups under
#'   that header; the parquet pivots the pipeline writes have no such row.
#' * it dropped the first row unconditionally, for the same historical reason,
#'   so the first instance of every real pivot would have been lost.
#' * `areas_selection = NULL`, the documented way to ask for every area, made
#'   the filter `%in% NULL`, which is `FALSE` everywhere. The default asked for
#'   nothing.
#' * it accumulated with `if (exists("res"))`, and `exists()` searches the
#'   calling frames: a variable named `res` in any caller would have been
#'   `rbind`-ed to the first sample's rows.
#'
#' It also called `sem_study_summary_get()`, which reads a file, and assigned the
#' result to a variable it never used.
#'
#' @section The sixth coordinate:
#' A pivot is identified by six coordinates, and `aggregation` is one of them: a
#' block holding more than one position was reduced somehow, and the name of the
#' file has to say how. This function used to pass four of the six to
#' `io_read_pivot()`, so for any area class whose blocks hold more than one
#' position - `GENE`, `DMR`, every real area - the file name could not even be
#' built and the read ended in `aggregation is required for scope 'INSTANCE'`.
#' It happened to work for the per-position classes only because there the
#' aggregation resolves to `VALUE` on its own.
#'
#' @section The phenotype of a sample:
#' The phenotype is looked up by `Sample_ID` and the **first** matching row wins.
#' A sample sheet may legitimately carry the same `Sample_ID` in more than one
#' row - a reference reused across groups is the normal case - and the type of
#' the column is preserved, so a numeric phenotype stays numeric and the caller
#' can choose a continuous scale for it.
#'
#' A sample column with no row in the sample sheet keeps its values and gets
#' `NA` as its phenotype, with the count logged. Dropping it would remove data
#' from a chart without the chart saying so; an `NA` level is visible.
#'
#' @param marker character. Marker name, e.g. `"MUTATIONS"`, `"DELTAS"`.
#' @param figure character. `"HYPER"`, `"HYPO"` or `"BOTH"`.
#' @param area character. Area class, e.g. `"GENE"`, `"PROBE"`, `"DMR"`.
#' @param subarea character. Subarea, e.g. `"WHOLE"`, `"TSS1500"`.
#' @param phenotype_column character. Column of the sample sheet to pair to each
#'   sample. Must exist in `sample_sheet`.
#' @param sample_sheet data.frame. Must carry `Sample_ID`.
#' @param areas_selection character vector. Instances to keep; `NULL`, the
#'   default, keeps all of them.
#' @param aggregation character. How a block holding more than one position was
#'   reduced: `"SUM"`, `"MEAN"`, `"MEDIAN"`, `"VARIANCE"`, `"IQR"`, `"MODELOW"`,
#'   `"MODEHIGH"`, `"VALUE"`. Required for any area whose blocks hold more than
#'   one position, which is every area class except the per-position ones; there
#'   it resolves on its own and may be left `NULL`.
#'
#' @return A data.frame with `AREA`, `VALUE`, `SAMPLE` and `phenotype`, or
#'   `NULL` when no pivot of that class exists or nothing survives the
#'   selection. `NULL` and a zero-row answer mean the same thing to a caller:
#'   there is nothing to draw.
#'
#' @keywords internal
#' @noRd
io_pivot_to_long_format <- function(marker, figure, area, subarea,
                                    phenotype_column, sample_sheet,
                                    areas_selection = NULL,
                                    aggregation = NULL)
{
  # A request naming a column the sheet does not have is malformed, and the
  # answer is a refusal rather than a column of NA: the chart would be drawn,
  # saved and indistinguishable from one whose phenotype was genuinely missing.
  if (!(phenotype_column %in% colnames(sample_sheet)))
    stop("phenotype_column '", phenotype_column,
         "' is not a column of the sample sheet. Available: ",
         paste(colnames(sample_sheet), collapse = ", "))

  area_pivot_lazy <- io_read_pivot(marker, figure, area, subarea,
                                   aggregation = aggregation)
  if (is.null(area_pivot_lazy))
    return(NULL)

  area_pivot <- as.data.frame(area_pivot_lazy$collect())
  if (nrow(area_pivot) == 0L || ncol(area_pivot) < 2L)
    return(NULL)

  # Column 1 holds the instance. It is `AREA` in every pivot the pipeline
  # writes; reading the name rather than asserting it keeps this working for a
  # pivot whose instance column is called something else.
  instance_column <- colnames(area_pivot)[1]
  sample_columns  <- colnames(area_pivot)[-1]

  if (!is.null(areas_selection)) {
    keep <- area_pivot[[instance_column]] %in% areas_selection
    area_pivot <- area_pivot[keep, , drop = FALSE]
    if (nrow(area_pivot) == 0L)
      return(NULL)
  }

  # match() gives the FIRST matching row, which is what the reuse of a reference
  # across groups requires, and indexing the column rather than subsetting the
  # data.frame preserves its type.
  sheet_row <- match(sample_columns, as.character(sample_sheet$Sample_ID))
  phenotype_values <- sample_sheet[[phenotype_column]][sheet_row]

  unknown <- sum(is.na(sheet_row))
  if (unknown > 0L)
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " ", unknown, " of ", length(sample_columns),
              " sample columns of the ", marker, "/", figure, " pivot have no",
              " row in the sample sheet: their phenotype is NA.")

  n_instances <- nrow(area_pivot)
  data.frame(
    AREA      = rep(as.character(area_pivot[[instance_column]]),
                    times = length(sample_columns)),
    VALUE     = as.numeric(unlist(area_pivot[, sample_columns, drop = FALSE],
                                   use.names = FALSE)),
    SAMPLE    = rep(sample_columns, each = n_instances),
    phenotype = rep(phenotype_values, each = n_instances),
    stringsAsFactors = FALSE)
}
