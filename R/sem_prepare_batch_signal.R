#' Prepare a SIGNAL matrix for one sem_analyze_batch() invocation
#'
#' Single source of truth for normalising a raw input matrix into a shape
#' that is consistent with the probe-level annotation used downstream. The
#' invariant guaranteed by this function is:
#'
#'   `nrow(signal_data) == nrow(attr(signal_data, "probe_features"))`
#'
#' and `rownames(signal_data)` is exactly `attr(., "probe_features")$PROBE`
#' in the same order.
#'
#' This replaces ~30 lines of scattered annotation/filter/align logic that
#' previously lived in `sem_analyze_batch()` fresh-path and was a source of
#' silent drift between the SIGNAL matrix and the probe_features used to
#' compute thresholds, write the POSITION pivot, and run downstream
#' analyses. The classic failure (visible v35-v43) was
#'
#'   `Error in data.frame(probe_features, VALUE = values,
#'                        row.names = probe_features$PROBE):
#'    arguments imply differing number of rows: 366948, 374705`
#'
#' which the alignment step in this function eliminates by construction.
#'
#' @param signal_data A data.frame whose rownames are probe identifiers
#'   (Illumina probe IDs like `"cg00050873"` for K27/K450/K850, or
#'   coordinate-encoded `"<CHR>_<START>"` strings for WGBS/LONGREAD).
#'   Sample columns follow. Must be PROBE-keyed, i.e. already passed
#'   through `io_normalize_signal_input()`.
#' @param tech Character scalar. One of `"K27"`, `"K450"`, `"K850"`,
#'   `"WGBS"`, `"LONGREAD"`. If `NULL` (default) the function calls
#'   `core_get_meth_tech()` to detect it.
#' @param sex_chromosome_remove Logical. If `TRUE` (default), drop probes
#'   on `CHR == "X"` or `CHR == "Y"` from both `probe_features` and
#'   `signal_data`. Applied uniformly across all techs.
#'
#' @return The input `signal_data` filtered to the autosomal manifest
#'   intersection, with two attributes attached:
#'   \itemize{
#'     \item `probe_features` - data.frame with columns
#'       `PROBE, CHR, START, END` (+ any extra columns produced by the
#'       tech-specific annotation builder), one row per surviving probe,
#'       in the same order as `rownames(signal_data)`.
#'     \item `tech` - the resolved tech string.
#'   }
#'
#' @section Invariants:
#' Post-call:
#' \itemize{
#'   \item `nrow(signal_data) == nrow(attr(signal_data, "probe_features"))`
#'   \item `identical(rownames(signal_data),
#'                    attr(signal_data, "probe_features")$PROBE)`
#'   \item `!any(attr(signal_data, "probe_features")$CHR %in% c("X","Y"))`
#'         when `sex_chromosome_remove = TRUE`
#'   \item `!anyDuplicated(attr(signal_data, "probe_features")$PROBE)`
#' }
#'
#' @keywords internal
sem_prepare_batch_signal <- function(signal_data,
                                 tech = NULL,
                                 sex_chromosome_remove = TRUE) {

  # ---- 1. Resolve tech --------------------------------------------------
  if (is.null(tech) || !nzchar(tech)) {
    ssEnv <- core_get_meth_tech(signal_data)
    tech  <- ssEnv$tech
  }
  if (is.null(tech) || !nzchar(tech))
    stop("sem_prepare_batch_signal: could not resolve methylation technology.")

  # ---- 2. Build probe_features (tech-specific) --------------------------
  # Illumina path uses the Bioconductor manifest (with CpG island, gene,
  # cytoband, DMR annotation). WGBS / LONGREAD path uses the coordinate-
  # encoded PROBE IDs directly. Both routes produce a data.frame keyed by
  # PROBE with CHR / START / END columns at minimum.
  if (tech %in% c("WGBS", "LONGREAD")) {
    probe_features <- io_coord_probe_features(rownames(signal_data))
  } else {
    probe_features <- anno_probe_annotation_build(tech)
    # Drop the tech-specific TRUE/FALSE flag column once it's served its
    # purpose; downstream consumers only need the coordinate / annotation
    # columns.
    if (tech %in% colnames(probe_features))
      probe_features <- probe_features[, setdiff(colnames(probe_features), tech),
                                        drop = FALSE]
  }

  # ---- 3. Collapse duplicate PROBE rows --------------------------------
  # anno_probe_annotation_build() does a left-join with dmr_annotation, and
  # dmr_annotation has ~744 duplicate PROBE entries (one CpG can belong
  # to multiple DMRs). Without collapsing, the merge inflates
  # probe_features and the downstream `nrow(signal_data) ==
  # nrow(probe_features)` invariant fails by coincidence.
  if (anyDuplicated(probe_features$PROBE) > 0L) {
    probe_features <- probe_features[!duplicated(probe_features$PROBE),
                                      , drop = FALSE]
  }

  # ---- 4. Intersect with the input -------------------------------------
  # Keep only probes present in BOTH the manifest and the input signal.
  probe_features <- probe_features[probe_features$PROBE %in% rownames(signal_data),
                                    , drop = FALSE]

  # ---- 5. Sex-chromosome removal (uniform for all techs) ---------------
  if (isTRUE(sex_chromosome_remove)) {
    probe_features <- probe_features[!(probe_features$CHR %in% c("X", "Y")),
                                      , drop = FALSE]
  }

  if (nrow(probe_features) == 0L)
    stop("sem_prepare_batch_signal: probe_features became empty after ",
         "intersection / sex-chr filter. Check input PROBE format and tech.")

  # ---- 6. Align signal_data to probe_features --------------------------
  # Single subset that BOTH (a) keeps only probes present in
  # probe_features and (b) reorders signal_data to match the
  # probe_features row order. After this line the two data structures
  # share a strict bijection of rows.
  signal_data <- signal_data[probe_features$PROBE, , drop = FALSE]

  # ---- 7. Sanity checks (cheap, catch silent contract breaks) ----------
  stopifnot(
    "sem_prepare_batch_signal: row count mismatch after alignment" =
      nrow(signal_data) == nrow(probe_features),
    "sem_prepare_batch_signal: row order mismatch after alignment" =
      identical(rownames(signal_data), as.character(probe_features$PROBE))
  )

  # ---- 8. Convenience AREA_SUBAREA columns -----------------------------
  # PROBE_WHOLE and CHR_WHOLE are the canonical (AREA, SUBAREA) names
  # for the PROBE-level and CHR-level baselines (same pattern as
  # GENE_BODY, ISLAND_N_SHORE, CHR_CYTOBAND, etc.). Downstream code
  # iterates over (AREA, SUBAREA) combinations and does dynamic lookup
  # via `probe_features[[area_subarea]]`, so these columns must exist
  # even when their value is a 1:1 alias of PROBE / "chr"+CHR.
  if (!"PROBE_WHOLE" %in% colnames(probe_features))
    probe_features$PROBE_WHOLE <- probe_features$PROBE
  if (!"CHR_WHOLE" %in% colnames(probe_features))
    probe_features$CHR_WHOLE <- paste0("chr", probe_features$CHR)

  # ---- 9. Attach annotation as attributes ------------------------------
  attr(signal_data, "probe_features") <- probe_features
  attr(signal_data, "tech")            <- tech
  signal_data
}
