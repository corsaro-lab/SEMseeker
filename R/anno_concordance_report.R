#' Cross-pipeline annotation concordance report (Illumina manifest vs WGBS)
#'
#' Validates \code{\link{anno_area_granges_build}} (the WGBS / long-read annotation
#' path) against the Illumina manifest annotation used by
#' \code{\link{anno_probe_annotation_build}}. The two pipelines are expected to
#' assign CpGs to the same semantic areas (\code{GENE_*}, \code{ISLAND_*},
#' \code{CHR_CYTOBAND}, \code{DMR_*}). This function takes a subset of Illumina
#' probes of known annotation, runs the WGBS pipeline on the same coordinates,
#' and reports the concordance rate per area.
#'
#' Concordance categories (informational, returned in the \code{category} column):
#' \itemize{
#'   \item \code{"bundled"} - both pipelines use the same bundled data
#'     (\code{cytoband_hg19.rda}, \code{dmr_annotation.rda}). Expected rate: 1.0.
#'   \item \code{"txdb"} - WGBS path uses \code{TxDb} UCSC \code{knownGene}
#'     while the Illumina manifest uses RefSeq. Some probes map to different
#'     gene symbols between the two sources.
#'   \item \code{"annotationhub"} - WGBS path uses CpG island BED from
#'     \code{AnnotationHub} (UCSC). Expected to be near-identical to the
#'     manifest but boundary probes may shift across shore/shelf strata.
#' }
#'
#' Two concordance rates are returned:
#' \itemize{
#'   \item \code{concordance_rate_strict}: label sets (split on \code{;})
#'     must match exactly (\code{setequal}).
#'   \item \code{concordance_rate_intersection}: non-empty set intersection
#'     (at least one shared gene / island between the two pipelines).
#' }
#' NA-only probes (missing in both pipelines) are counted in \code{n_both_na}
#' and excluded from the rate denominators to avoid inflating the score.
#'
#' @param tech Character: one of \code{"K27"}, \code{"K450"}, \code{"K850"}.
#'   Default \code{"K850"}.
#' @param n_probes Integer (or \code{NULL} for all probes). Default 1000L.
#' @param genome_build Character: \code{"hg19"} (default), \code{"hg38"}, or
#'   \code{"mm10"}.
#' @param areas Character vector of \code{area_subarea} names to benchmark
#'   (e.g. \code{"GENE_BODY"}, \code{"ISLAND_N_SHORE"}). \code{NULL} (default)
#'   = all supported areas.
#' @param seed Integer, random seed for probe subsampling. Default 42L.
#' @param csv_out Optional path where the report is written via
#'   \code{write.csv2()}. \code{NULL} (default) = no file.
#'
#' @return A \code{data.frame} with one row per area and columns:
#'   \code{area}, \code{category}, \code{n_probes}, \code{n_both_na},
#'   \code{n_both_labeled}, \code{n_only_illumina}, \code{n_only_wgbs},
#'   \code{n_label_match_strict}, \code{n_label_match_intersection},
#'   \code{concordance_rate_strict}, \code{concordance_rate_intersection}.
#'
#' @seealso \code{\link{anno_probe_annotation_build}},
#'   \code{\link{anno_area_granges_build}}.
#'
#' @keywords internal
anno_concordance_report <- function(tech         = "K850",
                                          n_probes     = 1000L,
                                          genome_build = "hg19",
                                          areas        = NULL,
                                          seed         = 42L,
                                          csv_out      = NULL) {

  for (pkg in c("GenomicRanges", "IRanges", "S4Vectors"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("anno_concordance_report(): package '", pkg,
           "' is required.")

  if (is.null(areas)) {
    areas <- c(
      "GENE_BODY", "GENE_TSS200", "GENE_TSS1500", "GENE_1STEXON",
      "GENE_5UTR", "GENE_3UTR",
      "ISLAND_WHOLE", "ISLAND_ISLAND",
      "ISLAND_N_SHORE", "ISLAND_S_SHORE",
      "ISLAND_N_SHELF", "ISLAND_S_SHELF", "ISLAND_OPENSEA",
      "CHR_CYTOBAND", "DMR_WHOLE", "DMR_DMR"
    )
  }

  # --- Ground truth: Illumina manifest annotation ------------------------
  anno <- anno_probe_annotation_build(tech)

  # Bioconductor packages must not call set.seed() in package code (it
  # mutates the user's RNG state); withr::with_seed scopes the seed to
  # the sampling expression only.
  n_total <- nrow(anno)
  idx <- if (is.null(n_probes) || n_probes >= n_total) seq_len(n_total)
         else withr::with_seed(seed, sort(sample.int(n_total, n_probes)))
  subset_df <- anno[idx, , drop = FALSE]

  # The CpG GRanges is rebuilt per area because different area_gr builders use
  # different seqname styles ("chr1" vs "1") depending on their data source
  # (AnnotationHub uses "chr1"; bundled cytoband/dmr tables use "1").
  # .anno_wgbs_labels_for_area() aligns seqnames to match the area's style.

  # --- Per-area comparison -----------------------------------------------
  rows <- lapply(areas, function(area) {
    .anno_concordance_for_area(area, subset_df, genome_build)
  })
  report <- do.call(rbind, rows)

  if (!is.null(csv_out))
    utils::write.csv2(report, csv_out, row.names = FALSE)

  report
}


# =========================================================================
# Internal helpers
# =========================================================================

#' @keywords internal
.anno_concordance_for_area <- function(area_subarea, subset_df, genome_build) {

  category <- .anno_area_category(area_subarea)

  # --- Illumina labels (from manifest column) ---------------------------
  illumina_raw <- if (area_subarea %in% colnames(subset_df))
    subset_df[[area_subarea]]
  else
    rep(NA_character_, nrow(subset_df))

  illumina_sets <- lapply(illumina_raw, .anno_normalize_label_set)

  # --- WGBS labels (via anno_area_granges_build + findOverlaps) --------------
  wgbs_sets <- .anno_wgbs_labels_for_area(area_subarea, subset_df, genome_build)

  # --- Per-CpG comparison -----------------------------------------------
  n       <- length(illumina_sets)
  il_empty <- vapply(illumina_sets, function(s) length(s) == 0L, logical(1))
  wg_empty <- vapply(wgbs_sets,     function(s) length(s) == 0L, logical(1))

  n_both_na       <- sum(il_empty & wg_empty)
  n_only_illumina <- sum(!il_empty & wg_empty)
  n_only_wgbs     <- sum(il_empty & !wg_empty)
  n_both_labeled  <- sum(!il_empty & !wg_empty)

  match_strict <- vapply(seq_len(n), function(i) {
    !il_empty[i] && !wg_empty[i] && setequal(illumina_sets[[i]], wgbs_sets[[i]])
  }, logical(1))
  match_inters <- vapply(seq_len(n), function(i) {
    !il_empty[i] && !wg_empty[i] &&
      length(intersect(illumina_sets[[i]], wgbs_sets[[i]])) > 0L
  }, logical(1))

  denom <- n_both_labeled
  rate_strict <- if (denom > 0L) sum(match_strict) / denom else NA_real_
  rate_inters <- if (denom > 0L) sum(match_inters) / denom else NA_real_

  data.frame(
    area                             = area_subarea,
    category                         = category,
    n_probes                         = n,
    n_both_na                        = n_both_na,
    n_both_labeled                   = n_both_labeled,
    n_only_illumina                  = n_only_illumina,
    n_only_wgbs                      = n_only_wgbs,
    n_label_match_strict             = sum(match_strict),
    n_label_match_intersection       = sum(match_inters),
    concordance_rate_strict          = rate_strict,
    concordance_rate_intersection    = rate_inters,
    stringsAsFactors                 = FALSE
  )
}

#' @keywords internal
.anno_wgbs_labels_for_area <- function(area_subarea, subset_df, genome_build) {

  area_gr <- tryCatch(
    anno_area_granges_build(area_subarea, genome_build = genome_build),
    error = function(e) NULL
  )
  if (is.null(area_gr))
    return(rep(list(character(0)), nrow(subset_df)))

  # Align seqname style ("chr1" vs "1") to the area_gr's own style.
  # bundled cytoband / dmr tables use "1"; AnnotationHub + TxDb use "chr1".
  area_seqs <- as.character(GenomicRanges::seqnames(area_gr))
  use_chr_prefix <- length(area_seqs) > 0L && any(grepl("^chr", area_seqs))
  probe_chr <- if (use_chr_prefix) paste0("chr", subset_df$CHR) else subset_df$CHR

  cpg_gr <- GenomicRanges::GRanges(
    seqnames = probe_chr,
    ranges   = IRanges::IRanges(start = subset_df$START, end = subset_df$START)
  )

  hits <- suppressWarnings(
    GenomicRanges::findOverlaps(cpg_gr, area_gr, ignore.strand = TRUE))
  labels <- as.character(GenomicRanges::mcols(area_gr)$label)

  # Build a per-CpG list of label sets
  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)
  out <- replicate(length(cpg_gr), character(0), simplify = FALSE)
  if (length(qh)) {
    by_q <- split(labels[sh], qh)
    for (q in names(by_q)) {
      qi <- as.integer(q)
      out[[qi]] <- .anno_normalize_label_set(by_q[[q]])
    }
  }
  out
}

#' @keywords internal
.anno_normalize_label_set <- function(x) {
  if (length(x) == 0L) return(character(0))
  # Accept vector of strings, possibly `;`-delimited
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(character(0))
  parts <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  # Strip "chr" prefix from island coordinates for consistency
  parts <- sub("^chr", "", parts, ignore.case = TRUE)
  parts <- toupper(parts)    # case-insensitive gene symbol comparison
  sort(unique(parts))
}

#' @keywords internal
.anno_area_category <- function(area_subarea) {
  parts <- strsplit(area_subarea, "_", fixed = TRUE)[[1]]
  area  <- parts[1]
  switch(area,
    GENE   = "txdb",
    ISLAND = "annotationhub",
    CHR    = "bundled",        # CHR_CYTOBAND uses bundled cytoband table
    DMR    = "bundled",        # DMR_* uses bundled dmr_annotation table
    "unknown"
  )
}
