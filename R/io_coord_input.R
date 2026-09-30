# coord_input.R - Transparent conversion of coordinate-based input
# (WGBS and long-read) to SEMseeker's internal probe-ID-indexed format.
#
# Design:
#   Illumina input: rownames = probe IDs (e.g. "cg00000029")
#   WGBS/LONGREAD:  columns CHR, START, [END,] sample1, sample2, ...
#
#   io_normalize_signal_input() is called transparently inside sem_analyze_batch()
#   before core_get_meth_tech(). If the input has CHR/START columns, it converts
#   the data frame to a probe-ID-indexed matrix using a synthetic probe ID:
#
#     synthetic probe ID = "{CHR_no_prefix}_<START>"
#     e.g. "chr1" / 10000  →  "1_10000"
#          "chrX" / 543200 →  "X_543200"
#
#   CHR normalisation: "chr1" → "1", "chrX" → "X" (strips leading "chr").
#   This matches the convention already used everywhere in SEMseeker where CHR
#   is stored without the "chr" prefix.
#
#   Area support for WGBS/LONGREAD in this release (B-01/B-02):
#     Supported:  CHR_WHOLE, PROBE_WHOLE (coordinate-based, no annotation)
#     Not yet:    GENE_*, ISLAND_*, DMR_* → requires anno_area_granges_build() (C-04)

# ---------------------------------------------------------------------------
# Column-name aliases recognised for coordinate columns
# ---------------------------------------------------------------------------

.COORD_CHR_NAMES   <- c("CHR", "chr", "chrom", "CHROM",
                         "chromosome", "Chromosome", "seqnames")
.COORD_START_NAMES <- c("START", "start", "chromStart",
                         "pos", "POS", "position", "Position")
.COORD_END_NAMES   <- c("END", "end", "chromEnd")

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Find coordinate column names in df.
# Returns list(chr=, start=, end=) or NULL if not a coordinate frame.
.io_find_coord_cols <- function(df) {
  chr_col   <- intersect(.COORD_CHR_NAMES,   colnames(df))
  start_col <- intersect(.COORD_START_NAMES, colnames(df))
  end_col   <- intersect(.COORD_END_NAMES,   colnames(df))
  if (length(chr_col) == 0 || length(start_col) == 0) return(NULL)
  list(
    chr   = chr_col[1],
    start = start_col[1],
    end   = if (length(end_col) > 0) end_col[1] else NA_character_
  )
}

# Build synthetic probe ID.
.io_make_probe_id <- function(chr_vec, start_vec) {
  paste0(anno_normalize_chr(chr_vec, "internal"), "_", as.character(start_vec))
}

# ---------------------------------------------------------------------------
# Exported-internal functions
# ---------------------------------------------------------------------------

#' Detect whether a data frame is in coordinate format (CHR + START columns).
#'
#' @param df A data frame.
#' @return Logical scalar.
io_is_coord_format <- function(df) {
  is.data.frame(df) && !is.null(.io_find_coord_cols(df))
}

#' Convert a coordinate-based methylation data frame to SEMseeker internal format.
#'
#' The input is a wide data frame where each row is a CpG position:
#'   CHR | START | [END] | sample1 | sample2 | ...
#'
#' The output is a data frame with synthetic probe IDs as rownames and sample
#' columns only (CHR/START/END are consumed and encoded in the rowname).
#'
#' @param df Data frame with CHR and START columns (END optional).
#' @return Data frame with rownames = synthetic probe IDs, columns = samples.
io_coord_to_semseeker <- function(df) {
  cols <- .io_find_coord_cols(df)
  if (is.null(cols))
    stop("io_coord_to_semseeker: no CHR/START columns found in input.")

  probe_ids <- .io_make_probe_id(df[[cols$chr]], df[[cols$start]])

  drop_cols <- c(cols$chr, cols$start)
  if (!is.na(cols$end)) drop_cols <- c(drop_cols, cols$end)
  drop_cols <- intersect(drop_cols, colnames(df))

  result <- df[, !colnames(df) %in% drop_cols, drop = FALSE]
  rownames(result) <- probe_ids
  result
}

#' Transparently normalise signal input for any technology.
#'
#' If the input data frame has CHR/START columns (WGBS or long-read format),
#' it is converted to the probe-ID-indexed format expected by SEMseeker.
#' Illumina matrices (rownames = probe IDs) are returned unchanged.
#'
#' Called inside \code{sem_analyze_batch()} before \code{core_get_meth_tech()}.
#'
#' @param signal_data A data frame (probe-indexed or coordinate-based).
#' @return A data frame with rownames = probe IDs (real or synthetic).
io_normalize_signal_input <- function(signal_data) {
  if (io_is_coord_format(signal_data)) {
    # Logging is done by the caller (sem_analyze_batch) which has an active session.
    signal_data <- io_coord_to_semseeker(signal_data)
  }
  signal_data
}

# ---------------------------------------------------------------------------
# Parsing synthetic probe IDs back to genomic coordinates
# ---------------------------------------------------------------------------

#' Parse synthetic probe IDs back to a CHR / START / END data frame.
#'
#' Synthetic probe ID format: "<CHR>_<START>" where CHR has no "chr" prefix.
#' E.g. "1_10000" → CHR = "1", START = 10000L, END = 10001L.
#'
#' @param probe_ids Character vector of synthetic probe IDs.
#' @return data.frame with columns CHR (character), START (integer), END (integer).
io_probe_id_to_coord <- function(probe_ids) {
  # Split on the LAST underscore (START is always a pure integer suffix).
  chr_part   <- sub("_[0-9]+$", "", probe_ids)
  start_part <- as.integer(sub("^.*_([0-9]+)$", "\\1", probe_ids))
  data.frame(
    CHR   = chr_part,
    START = start_part,
    END   = start_part + 1L,
    stringsAsFactors = FALSE
  )
}

#' Build a minimal probe_features data frame from synthetic probe IDs.
#'
#' Used inside \code{sem_analyze_batch()} for WGBS/LONGREAD data in place of the
#' Bioconductor-annotation-based \code{anno_probe_features_get("PROBE")} call.
#'
#' @param probe_ids Character vector of synthetic probe IDs.
#' @return data.frame with columns PROBE, CHR, START, END.
io_coord_probe_features <- function(probe_ids) {
  coords <- io_probe_id_to_coord(probe_ids)
  data.frame(
    PROBE = probe_ids,
    CHR   = coords$CHR,
    START = coords$START,
    END   = coords$END,
    stringsAsFactors = FALSE
  )
}
