# island_opensea.R - shared CpG-island context semantics
#
# Single source of truth for the ISLAND subarea definitions used by BOTH
# annotation backends, so they cannot drift apart:
#   * anno_probe_annotation_build()  - Illumina array packages (Relation_to_Island)
#   * anno_area_granges_build()       - coordinate / AnnotationHub (WGBS, long-read)
#
# Illumina's Relation_to_Island has six categories: Island, N_Shore, S_Shore,
# N_Shelf, S_Shelf, OpenSea. Historically SEMseeker kept only the four
# shore/shelf contexts plus a WHOLE that actually meant "Island core" - the
# explicit Island core and the OpenSea compartment were lost. This module
# restores both and aligns the semantics:
#
#   WHOLE    = whole CpG-island neighbourhood = Island + shores + shelves
#              (mirrors GENE_WHOLE = whole gene). Labelled by island coordinate.
#   ISLAND   = island core only (Relation_to_Island == "Island").
#   N_SHORE / S_SHORE / N_SHELF / S_SHELF = ±2 kb / ±2-4 kb flanks (unchanged).
#   OPENSEA  = CpGs outside every neighbourhood. OpenSea regions are the
#              genomic gaps BETWEEN neighbourhoods, each labelled by its own
#              coordinate "chr:start-end" - exactly as islands are labelled by
#              their coordinate. A gap never crosses a chromosome boundary, so
#              no OpenSea region spans chromosomes.
#
# Neighbourhood = island core extended by .ISLAND_FLANK on each side
# (2 kb shore + 2 kb shelf = 4 kb).

# Shore (2 kb) + shelf (2 kb) flank that defines the island neighbourhood.
.ISLAND_FLANK <- 4000L

#' OpenSea regions as the genomic gaps between island neighbourhoods.
#'
#' Computed per chromosome with IRanges so it does not depend on GRanges
#' seqlengths / GenomeInfoDb. Each gap stays within one chromosome, so no
#' OpenSea region ever spans a chromosome boundary.
#'
#' @param island_gr GRanges of island cores.
#' @param flank Integer half-width of the neighbourhood (default 4 kb).
#' @param chrom_ends Optional named (by seqname) integer vector/list giving the
#'   universe end per chromosome (e.g. seqlength, or the furthest CpG). Trailing
#'   gaps past the last neighbourhood only appear for chromosomes listed here;
#'   otherwise the universe ends at the last neighbourhood.
#' @return GRanges of OpenSea regions with \code{mcols()$label} = "chr:start-end".
#' @keywords internal
#' @noRd
.anno_opensea_gaps <- function(island_gr, flank = .ISLAND_FLANK, chrom_ends = NULL) {
  chrs   <- as.character(GenomicRanges::seqnames(island_gr))
  starts <- GenomicRanges::start(island_gr)
  ends   <- GenomicRanges::end(island_gr)

  out_chr <- character(0); out_start <- integer(0); out_end <- integer(0)
  for (chr in unique(chrs)) {
    sel <- chrs == chr
    nb  <- IRanges::reduce(IRanges::IRanges(
      pmax(1L, starts[sel] - flank), ends[sel] + flank))
    uend <- if (!is.null(chrom_ends) && !is.null(chrom_ends[[chr]]) &&
                !is.na(chrom_ends[[chr]]))
              as.integer(chrom_ends[[chr]]) else max(IRanges::end(nb))
    g <- IRanges::gaps(nb, start = 1L, end = uend)
    if (length(g) == 0L) next
    out_chr   <- c(out_chr, rep(chr, length(g)))
    out_start <- c(out_start, IRanges::start(g))
    out_end   <- c(out_end, IRanges::end(g))
  }

  gr <- GenomicRanges::GRanges(out_chr, IRanges::IRanges(out_start, out_end))
  GenomicRanges::mcols(gr)$label <- paste0(out_chr, ":", out_start, "-", out_end)
  gr
}

#' Assign each probe/CpG the coordinate label of the OpenSea gap containing it.
#'
#' Used by the Illumina backend, where annotation is probe-centric. Probes that
#' fall inside a neighbourhood (i.e. are not OpenSea) overlap no gap and get NA.
#'
#' @param probe_chr Character vector of seqnames (must match \code{island_gr}'s
#'   seqlevel style, e.g. both "chr1" or both "1").
#' @param probe_pos Integer vector of CpG positions.
#' @param island_gr GRanges of island cores.
#' @return Character vector of OpenSea region labels (NA where not OpenSea).
#' @keywords internal
#' @noRd
.anno_assign_opensea_labels <- function(probe_chr, probe_pos, island_gr) {
  # Extend the per-chromosome universe to the furthest CpG so trailing OpenSea
  # probes land in a gap rather than falling off the last neighbourhood.
  max_per_chr <- tapply(probe_pos, probe_chr, max, na.rm = TRUE)
  chrom_ends  <- as.list(max_per_chr + .ISLAND_FLANK + 1L)

  gaps_gr  <- .anno_opensea_gaps(island_gr, chrom_ends = chrom_ends)
  probe_gr <- GenomicRanges::GRanges(
    probe_chr, IRanges::IRanges(probe_pos, probe_pos))

  hit <- GenomicRanges::findOverlaps(probe_gr, gaps_gr, select = "first")
  out <- rep(NA_character_, length(probe_gr))
  ok  <- !is.na(hit)
  out[ok] <- GenomicRanges::mcols(gaps_gr)$label[hit[ok]]
  out
}

#' Build the seven ISLAND_* annotation columns from Illumina island context.
#'
#' Pure (no Bioconductor annotation package needed): operates on the four
#' per-probe vectors the Illumina manifest provides. Shared by
#' \code{anno_probe_annotation_build()} and unit tests so the ISLAND/OPENSEA recoding
#' is exercised without loading an array annotation package.
#'
#' @param island_rel Character vector: \code{Relation_to_Island} per probe.
#' @param island_name Character vector: \code{Islands_Name} per probe
#'   ("chr:start-end", shared by an island and its shores/shelves).
#' @param chr Character vector: chromosome WITHOUT the "chr" prefix (as in
#'   \code{anno_df$CHR}).
#' @param start Integer vector: CpG position.
#' @return Named list of seven character vectors: \code{ISLAND_WHOLE},
#'   \code{ISLAND_ISLAND}, \code{ISLAND_N_SHORE}, \code{ISLAND_S_SHORE},
#'   \code{ISLAND_N_SHELF}, \code{ISLAND_S_SHELF}, \code{ISLAND_OPENSEA}.
#' @keywords internal
#' @noRd
.anno_island_columns <- function(island_rel, island_name, chr, start) {
  island_rel  <- as.character(island_rel)
  island_name <- as.character(island_name)
  ctx <- c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf")

  out <- list(
    # WHOLE = whole neighbourhood (core + shores + shelves), like GENE_WHOLE.
    ISLAND_WHOLE   = ifelse(island_rel %in% ctx,    island_name, NA_character_),
    ISLAND_ISLAND  = ifelse(island_rel == "Island",  island_name, NA_character_),
    ISLAND_N_SHORE = ifelse(island_rel == "N_Shore", island_name, NA_character_),
    ISLAND_S_SHORE = ifelse(island_rel == "S_Shore", island_name, NA_character_),
    ISLAND_N_SHELF = ifelse(island_rel == "N_Shelf", island_name, NA_character_),
    ISLAND_S_SHELF = ifelse(island_rel == "S_Shelf", island_name, NA_character_),
    ISLAND_OPENSEA = rep(NA_character_, length(island_rel))
  )

  # OPENSEA: label each open-sea CpG by the inter-neighbourhood gap holding it.
  is_opensea <- island_rel == "OpenSea" & !is.na(island_rel)
  if (any(is_opensea)) {
    island_gr <- .anno_islands_gr_from_names(island_name)
    if (length(island_gr) > 0L) {
      out$ISLAND_OPENSEA[is_opensea] <- .anno_assign_opensea_labels(
        probe_chr = paste0("chr", chr[is_opensea]),
        probe_pos = start[is_opensea],
        island_gr = island_gr
      )
    }
  }
  out
}

#' Parse Illumina \code{Islands_Name} coordinate strings into a GRanges.
#'
#' \code{Islands_Name} is the UCSC CpG-island coordinate (e.g.
#' "chr1:2004858-2005346") shared by an island and its shores/shelves.
#'
#' @param names_vec Character vector of "chrN:start-end" strings (NA/"" ignored).
#' @return GRanges of unique island cores (seqlevels carry the "chr" prefix).
#' @keywords internal
#' @noRd
.anno_islands_gr_from_names <- function(names_vec) {
  u <- unique(names_vec[!is.na(names_vec) & nzchar(names_vec)])
  m <- regmatches(u, regexec("^(chr[0-9XYM]+):([0-9]+)-([0-9]+)$", u))
  ok <- vapply(m, length, integer(1)) == 4L
  m <- m[ok]
  if (length(m) == 0L)
    return(GenomicRanges::GRanges())
  chr   <- vapply(m, `[[`, character(1), 2L)
  start <- as.integer(vapply(m, `[[`, character(1), 3L))
  end   <- as.integer(vapply(m, `[[`, character(1), 4L))
  GenomicRanges::GRanges(chr, IRanges::IRanges(start, end))
}
