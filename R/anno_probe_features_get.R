#' Retrieve probe feature annotations for a given genomic area
#'
#' Returns a \code{data.frame} of CpG probe coordinates and feature annotations
#' for the requested area/subarea combination.
#'
#' For \strong{Illumina} data (K850/K450/K27), annotations are built from
#' Bioconductor array annotation packages (see \code{\link{anno_probe_annotation_build}})
#' and cached in the session environment.
#'
#' For \strong{WGBS} and \strong{LONGREAD} data, coordinates are read directly
#' from the saved POSITION pivot parquet; semantic areas (GENE_*, ISLAND_*,
#' CHR_CYTOBAND, DMR_*) are resolved via \code{\link{anno_area_granges_build}} and
#' \code{GenomicRanges::findOverlaps()}.
#'
#' @section PROBE_WHOLE vs POSITION_WHOLE - technology semantics:
#' The area \code{PROBE_WHOLE} has different meanings depending on technology:
#' \describe{
#'   \item{Illumina}{Each row identifies a specific \emph{array probe} by its
#'     manufacturer ID (e.g. \code{cg00000029}).  The statistical test is
#'     performed at the individual-probe level.  Probe identity is meaningful
#'     here: two studies using the same array share the exact same set of probe
#'     IDs and can be directly compared.}
#'   \item{WGBS / LONGREAD}{There are no probe IDs.  \code{PROBE_WHOLE} is
#'     treated as \strong{\code{POSITION_WHOLE}}: each row identifies a CpG by
#'     its genomic coordinate (\code{CHR\_START}, e.g. \code{"1\_10000"}).
#'     The statistical test is performed at the individual-position level.
#'     Two WGBS datasets can be compared only if they share the same reference
#'     genome (\code{ssEnv\$genome_build}) - mismatches are detected by the
#'     session provenance guard (C-06).}
#' }
#' In both cases the downstream analysis pipeline is identical; the distinction
#' is purely in what the \code{PROBE} column \emph{means} to the researcher.
#'
#' Probes on sex chromosomes are removed when \code{ssEnv$sex_chromosome_remove}
#' is \code{TRUE}.
#'
#' @param area_subarea Character scalar: area and subarea joined by an
#'   underscore (e.g. \code{"GENE_BODY"}, \code{"CHR_WHOLE"},
#'   \code{"ISLAND_N_SHORE"}).  If no underscore is present, \code{"_WHOLE"}
#'   is appended automatically.
#'
#' @return A \code{data.frame} with columns \code{PROBE}, \code{CHR},
#'   \code{START}, \code{END}, and the requested feature column.
#'
anno_probe_features_get <- function(area_subarea) {

  ssEnv <- core_get_session_info()

  # Contract: sem_prepare_batch_signal() (fresh path) or core_get_meth_tech() (resume
  # path) should set ssEnv$tech before any anno_probe_features_get() call site.
  # Fallback: read the SIGNAL PROBE pivot and re-derive - kept for tests and
  # legacy callers that bypass sem_prepare_batch_signal(). Loud WARNING so the
  # drift is observable.
  if (is.null(ssEnv$tech) || ssEnv$tech == "") {
    core_log_event("WARNING: anno_probe_features_get() called before ssEnv$tech is set. ",
              "Falling back to lazy detection from SIGNAL PROBE pivot. ",
              "sem_prepare_batch_signal() should run first in the normal pipeline.")
    signal_pivot_lazy <- io_read_pivot("SIGNAL", io_signal_figure(), "PROBE", "WHOLE")
    if (is.null(signal_pivot_lazy))
      stop("SIGNAL PROBE pivot not available - cannot detect technology.")
    signal_data_r <- as.data.frame(signal_pivot_lazy$collect())
    if ("AREA" %in% colnames(signal_data_r)) {
      rownames(signal_data_r) <- signal_data_r$AREA
      signal_data_r$AREA <- NULL
    }
    ssEnv <- core_get_meth_tech(signal_data_r)
    if (is.null(ssEnv$tech) || ssEnv$tech == "")
      stop("anno_probe_features_get: could not determine array technology.")
  }

  if (!grepl("_", area_subarea))
    area_subarea <- paste0(area_subarea, "_WHOLE")

  # -----------------------------------------------------------------------
  # WGBS / LONGREAD path
  # Coordinate-only areas (CHR_WHOLE, PROBE_WHOLE) are handled inline.
  # All semantic areas (GENE_*, ISLAND_*, CHR_CYTOBAND, DMR_*) are built
  # via anno_area_granges_build() and CpG assignments resolved with findOverlaps().
  # -----------------------------------------------------------------------
  if (ssEnv$tech %in% c("WGBS", "LONGREAD")) {

    # read via unified dispatcher. NULL means neither cached
    # nor per-sample bed files exist for SIGNAL - same failure
    # mode as the previous file.exists() check.
    sig_pivot_lazy <- io_read_pivot("SIGNAL", io_signal_figure(), "POSITION", "WHOLE")
    if (is.null(sig_pivot_lazy))
      stop("POSITION pivot not found. Ensure io_signal_save() completed before ",
           "calling anno_probe_features_get() for area '", area_subarea, "'.")

    pos <- as.data.frame(sig_pivot_lazy$collect())[, c("CHR", "START", "END")]
    probe_ids <- paste0(pos$CHR, "_", pos$START)

    # --- Coordinate-only areas (no annotation needed) ---
    if (area_subarea %in% c("CHR_WHOLE", "PROBE_WHOLE",
                             "POSITION_WHOLE", "PROBE")) {
      probe_features <- data.frame(
        PROBE = probe_ids, CHR = pos$CHR, START = pos$START, END = pos$END,
        stringsAsFactors = FALSE)
      # Canonical (AREA, SUBAREA) columns: PROBE_WHOLE and CHR_WHOLE are
      # required by association/enrichment code that does dynamic
      # `probe_features[[area_subarea]]` lookup (same shape as GENE_BODY,
      # ISLAND_N_SHORE, CHR_CYTOBAND, DMR_*).
      if (grepl("CHR", area_subarea)) probe_features$CHR_WHOLE   <- paste0("chr", probe_features$CHR)
      if (grepl("PROBE", area_subarea)) probe_features$PROBE_WHOLE <- probe_features$PROBE
      if (isTRUE(ssEnv$sex_chromosome_remove))
        probe_features <- probe_features[
          !(probe_features$CHR %in% c("X", "Y")), ]
      return(probe_features)
    }

    # --- Semantic areas via anno_area_granges_build() + findOverlaps() ---
    for (pkg in c("GenomicRanges", "IRanges", "S4Vectors")) {
      if (!requireNamespace(pkg, quietly = TRUE))
        stop("Package '", pkg, "' is required for WGBS/LONGREAD area analysis.\n",
             "Install: BiocManager::install(c('GenomicRanges','IRanges','S4Vectors'))")
    }

    # Build GRanges for the CpG positions (seqnames need "chr" prefix for TxDb)
    cpg_gr <- GenomicRanges::GRanges(
      seqnames = paste0("chr", pos$CHR),
      ranges   = IRanges::IRanges(pos$START, pos$START)
    )

    # Build area GRanges (cached after first call per session)
    area_gr <- anno_area_granges_build(area_subarea,
                                  genome_build = ssEnv$genome_build)

    # Assign each CpG to its overlapping area (ignore strand for methylation)
    hits <- GenomicRanges::findOverlaps(cpg_gr, area_gr,
                                        ignore.strand = TRUE)

    q_hits <- S4Vectors::queryHits(hits)
    s_hits <- S4Vectors::subjectHits(hits)

    probe_features <- data.frame(
      PROBE = probe_ids[q_hits],
      CHR   = pos$CHR[q_hits],
      START = pos$START[q_hits],
      END   = pos$END[q_hits],
      stringsAsFactors = FALSE
    )
    probe_features[[area_subarea]] <-
      as.character(GenomicRanges::mcols(area_gr)$label[s_hits])

    if (nrow(probe_features) == 0)
      core_log_event("WARNING: no CpG positions overlap area '", area_subarea,
                "' for genome_build = '", ssEnv$genome_build, "'.")

    if (isTRUE(ssEnv$sex_chromosome_remove))
      probe_features <- probe_features[
        !(probe_features$CHR %in% c("X", "Y")), ]

    return(probe_features)
  }

  # -----------------------------------------------------------------------
  # Illumina path - Bioconductor annotation packages
  # -----------------------------------------------------------------------
  pkg <- .ANNO_PKGS[[ssEnv$tech]]
  if (is.null(pkg) || !requireNamespace(pkg, quietly = TRUE)) {
    stop("Annotation package '", pkg, "' is not installed. ",
         "Install it with: BiocManager::install('", pkg, "')")
  }
  probe_features <- anno_probe_annotation_build(ssEnv$tech)

  # Keep only probes matching the current technology
  probe_features <- probe_features[
    !is.na(probe_features[[ssEnv$tech]]) & probe_features[[ssEnv$tech]], ]
  probe_features$END <- probe_features$START

  if ((grepl("CHR", area_subarea) || grepl("PROBE", area_subarea)) &&
      !grepl("CHR_CYTOBAND", area_subarea)) {
    probe_features <- dplyr::distinct(
      probe_features[, c(ssEnv$tech, "PROBE", "CHR", "START", "END")])
  } else {
    cols_needed <- c(ssEnv$tech, "PROBE", "CHR", "START", "END", area_subarea)
    cols_needed <- intersect(cols_needed, colnames(probe_features))
    probe_features <- dplyr::distinct(probe_features[, cols_needed, drop = FALSE])
    # Dead code removed (2026-06-10): a stray
    #   signal_data <- signal_data[rownames(signal_data) %in% signal_data$PROBE, ]
    # was sitting here referencing a variable that does not exist in
    # anno_probe_features_get()'s scope. It never ran when ssEnv$tech was set
    # by the legacy lazy-detection path (the function returned early via
    # the PROBE / CHR branch), but with the AREA-based call sites added
    # by anno_annotate_position_pivots() - area_subarea = "GENE_BODY",
    # "ISLAND_N_SHORE", etc. - the else-branch is now reached and the
    # broken reference halts the run.
  }

  # Canonical (AREA, SUBAREA) columns required by downstream
  # association/enrichment code that iterates over area_subarea names
  # and does dynamic `probe_features[[area_subarea]]` lookup (same shape
  # as GENE_BODY, ISLAND_N_SHORE, CHR_CYTOBAND, DMR_*).
  if (grepl("CHR", area_subarea) && !grepl("CHR_CYTOBAND", area_subarea))
    probe_features$CHR_WHOLE <- paste0("chr", probe_features$CHR)

  if (grepl("PROBE", area_subarea))
    probe_features$PROBE_WHOLE <- probe_features$PROBE

  # Drop the technology flag column - not needed downstream
  probe_features <- probe_features[
    , -which(colnames(probe_features) %in% ssEnv$tech), drop = FALSE]

  # Remove sex chromosomes if requested
  if (isTRUE(ssEnv$sex_chromosome_remove))
    probe_features <- probe_features[
      !(probe_features$CHR %in% c("X", "Y")), ]

  return(probe_features)
}
