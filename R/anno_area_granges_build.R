# area_granges_build.R - Semantic area GRanges for WGBS / long-read data
#
# Reconstructs the same genomic boundaries that Illumina uses for its area
# definitions (TSS200, TSS1500, gene body, CpG islands, shores, shelves ...)
# using TxDb annotation packages and AnnotationHub, so that WGBS / long-read
# analyses use exactly the same region semantics as Illumina array analyses.
#
# Reference genome is taken from ssEnv$genome_build (set in core_init_env()) or
# passed explicitly; default "hg19" matches the Illumina annotation packages.
#
# NOTE - PROBE_WHOLE vs POSITION_WHOLE (technology semantics)
# -----------------------------------------------------------
# For Illumina data, "PROBE_WHOLE" means individual array probes identified
# by manufacturer IDs (e.g. cg00000029).  Statistical tests run at the
# probe level.  Probe IDs are comparable across studies using the same array.
#
# For WGBS / long-read data there are no probe IDs.  The unit of analysis is
# a genomic POSITION (CHR, START), encoded as a synthetic ID "CHR_START"
# (e.g. "1_10000").  PROBE_WHOLE is therefore treated as POSITION_WHOLE:
# each row identifies a CpG by coordinate, not by manufacturer identity.
# Cross-study comparisons require the same reference genome (ssEnv$genome_build).
#
# This function builds GRanges for all areas EXCEPT PROBE/POSITION_WHOLE,
# which are handled inline in anno_probe_features_get() without annotation.
#
# Supported area/subarea values:
#   GENE:    TSS200, TSS1500, 1STEXON, 5UTR, 3UTR, BODY, EXONBND, WHOLE
#   ISLAND:  WHOLE, ISLAND, N_SHORE, S_SHORE, N_SHELF, S_SHELF, OPENSEA
#   CHR:     WHOLE, CYTOBAND
#   DMR:     WHOLE, DMR
#   PROBE:   WHOLE  (coordinate-only, handled by anno_probe_features_get())
#
# All returned GRanges carry mcols()$label - the subarea identifier used
# downstream to group CpGs (gene symbol, island coordinate, cytoband name ...).

# ---------------------------------------------------------------------------
# Package-level in-memory cache (survives the R session, cleared on restart)
# ---------------------------------------------------------------------------
.area_gr_cache <- new.env(parent = emptyenv())

# ---------------------------------------------------------------------------
# TxDb helpers
# ---------------------------------------------------------------------------

.TXDB_PKGS <- list(
  hg19 = "TxDb.Hsapiens.UCSC.hg19.knownGene",
  hg38 = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  mm10 = "TxDb.Mmusculus.UCSC.mm10.knownGene"
)

.anno_get_txdb <- function(genome_build) {
  pkg <- .TXDB_PKGS[[genome_build]]
  if (is.null(pkg))
    stop("No TxDb package defined for genome_build = '", genome_build, "'. ",
         "Supported: ", paste(names(.TXDB_PKGS), collapse = ", "))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("TxDb package '", pkg, "' is not installed.\n",
         "Install with: BiocManager::install('", pkg, "')")
  get(pkg, envir = asNamespace(pkg))
}

# Share of unmappable keys above which the vector is treated as the wrong kind
# of identifier rather than as a few missing annotations.
#
# Measured on hg19: gene Entrez ids of the TxDb fail to map for 0.02 per cent
# of the set, 5 out of 28,622, which are ids that genuinely carry no symbol.
# The same call with transcript ids fails for 94 per cent. Anything between the
# two separates the cases; 20 per cent sits a thousandfold above the legitimate
# rate and far below the broken one.
.SYMBOL_UNMAPPED_MAX <- 0.20

# Map Entrez gene IDs to gene symbols.
#
# This used to fall back to the id, twice: to the id string when the symbol
# database was absent, and to the id of every key that did not map. Both were
# silent, and together they hid a caller that was passing transcript ids. The
# label it produced was a real gene symbol of an unrelated gene, because
# transcript ids and Entrez ids are both small integers and most of them
# collide. Nothing downstream could tell that label from a correct one.
#
# So the failures are loud now. The residual fallback stays for the handful of
# Entrez ids that carry no symbol: at 0.02 per cent it is the thing it was
# meant for.
.anno_entrez_to_symbol <- function(entrez_ids) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE))
    stop("a gene region has to be labelled with a gene symbol, and that needs ",
         "org.Hs.eg.db and AnnotationDbi. Without them the label would stay the ",
         "Entrez id the transcript database keys on, while the array path ",
         "labels with symbols, so the two would name the same gene two ways and ",
         "no result of one could be joined to the other. ",
         "BiocManager::install(c(\"org.Hs.eg.db\", \"AnnotationDbi\")).",
         call. = FALSE)

  entrez_ids <- as.character(entrez_ids)
  # AnnotationDbi does not return NAs when NONE of the keys are valid: it raises
  # .testForValidKeys("None of the keys entered are valid keys for 'ENTREZID'").
  # That is the same diagnosis as the threshold below, arrived at by a different
  # route, so it is caught here and reported in the same words. Leaving it to
  # propagate would be loud enough but would name AnnotationDbi's internals
  # instead of the mistake, which is that the caller passed the wrong kind of
  # identifier. The partial case - a vector where a few keys happen to be valid,
  # which is exactly what transcript ids look like - does not raise at all and is
  # what the threshold is for.
  syms <- tryCatch(
    suppressMessages(
      AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys      = entrez_ids,
        column    = "SYMBOL",
        keytype   = "ENTREZID",
        multiVals = "first"
      )
    ),
    error = function(e) {
      if (grepl("valid keys", conditionMessage(e), fixed = TRUE))
        return(stats::setNames(rep(NA_character_, length(entrez_ids)), entrez_ids))
      stop(e)
    })

  unmapped <- mean(is.na(syms))
  if (length(syms) > 0 && unmapped > .SYMBOL_UNMAPPED_MAX)
    stop(sprintf(
      paste0("%.1f%% of the identifiers passed for symbol lookup are not Entrez ",
             "gene ids of this database (%d of %d). That is not a gap in the ",
             "annotation, it is the wrong kind of key: either the caller passed ",
             "identifiers of something else, such as transcript ids, or the ",
             "genome build is not the one org.Hs.eg.db describes. Labelling ",
             "anyway is the worse outcome, because most wrong keys collide with ",
             "a valid Entrez id and come back as the symbol of an unrelated ",
             "gene. First few: %s."),
      100 * unmapped, sum(is.na(syms)), length(syms),
      paste(utils::head(entrez_ids[is.na(syms)], 5), collapse = ", ")),
      call. = FALSE)

  syms[is.na(syms)] <- entrez_ids[is.na(syms)]
  unname(syms)
}

#' Rename a per-transcript GRanges by the Entrez id of its parent gene
#'
#' Everything TxDb groups \code{by = "tx"} comes back named by transcript id:
#' \code{exonsBy(by = "tx")}, \code{fiveUTRsByTranscript()},
#' \code{threeUTRsByTranscript()}. The label block downstream reads
#' \code{names()} and maps it with \code{.anno_entrez_to_symbol()}, which
#' expects an Entrez gene id.
#'
#' Handing it a transcript id does not fail. Both are small integers, so most
#' transcript ids collide with a valid Entrez id and the mapper returns the
#' symbol of an unrelated gene: transcript 43 comes back as whatever gene
#' Entrez 43 is. The label is then a real, well-formed, wrong gene symbol, and
#' nothing downstream can tell it from a right one.
#'
#' Ranges whose transcript has no parent gene are dropped: a range that cannot
#' name its gene has nothing to contribute to a per-gene question.
#'
#' @param txdb a TxDb object.
#' @param gr a GRanges named by transcript id.
#' @return \code{gr}, named by Entrez gene id, without the unmappable ranges.
#' @keywords internal
#' @noRd
.anno_name_by_parent_gene <- function(txdb, gr) {
  by_gene <- GenomicFeatures::transcriptsBy(txdb, by = "gene")
  tx_ids  <- as.character(GenomicRanges::mcols(unlist(by_gene))$tx_id)
  gene_of <- stats::setNames(
    rep(names(by_gene), S4Vectors::elementNROWS(by_gene)), tx_ids)

  genes <- unname(gene_of[as.character(names(gr))])
  gr    <- gr[!is.na(genes)]
  names(gr) <- genes[!is.na(genes)]
  gr
}

# ---------------------------------------------------------------------------
# CpG island helper (AnnotationHub + disk cache)
# ---------------------------------------------------------------------------

.anno_get_cpg_islands <- function(genome_build) {
  if (!requireNamespace("AnnotationHub", quietly = TRUE))
    stop("AnnotationHub is required for island areas.\n",
         "Install with: BiocManager::install('AnnotationHub')")

  cache_dir  <- tools::R_user_dir("SEMseeker", "cache")
  cache_file <- file.path(cache_dir,
                          paste0("cpg_islands_", genome_build, ".rds"))

  if (file.exists(cache_file)) {
    core_log_event("INFO: loading cached CpG islands for ", genome_build,
              " from ", cache_file)
    return(readRDS(cache_file))
  }

  core_log_event("INFO: downloading CpG islands for ", genome_build,
            " from AnnotationHub (cached for future use).")
  ah <- AnnotationHub::AnnotationHub()
  # Search for UCSC cpgIslandExt track for the requested assembly
  q  <- AnnotationHub::query(ah, c("cpgIslandExt", genome_build))
  if (length(q) == 0)
    q <- AnnotationHub::query(ah, c("CpG Islands", genome_build, "UCSC"))
  if (length(q) == 0)
    stop("Could not find a CpG island resource for genome_build = '",
         genome_build, "' in AnnotationHub.")

  islands <- ah[[names(q)[1]]]

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(islands, cache_file)
  core_log_event("INFO: CpG islands cached at ", cache_file)
  islands
}

# ---------------------------------------------------------------------------
# Individual area builders - each returns a GRanges with mcols()$label
# ---------------------------------------------------------------------------

.anno_build_gene_area <- function(subarea, txdb) {
  # TSS = single-base GRanges at each gene's transcription start (strand-aware)
  all_genes <- GenomicFeatures::genes(txdb, single.strand.genes.only = FALSE)
  # genes() can return a GRangesList for multi-strand genes; keep only GRanges
  if (methods::is(all_genes, "GRangesList"))
    all_genes <- unlist(all_genes)

  # After unlist(), gene_ids live in names(all_genes); direct GRanges from
  # genes() also has them in mcols$gene_id. Prefer names when available.
  gene_ids <- if (!is.null(names(all_genes)) &&
                  length(names(all_genes)) == length(all_genes))
    names(all_genes)
  else
    GenomicRanges::mcols(all_genes)$gene_id
  symbols <- .anno_entrez_to_symbol(gene_ids)

  tss <- GenomicRanges::resize(all_genes, 1L, fix = "start")

  gr <- switch(subarea,
    TSS200 = {
      GenomicRanges::flank(tss, 200L)
    },
    TSS1500 = {
      # Ring: 201-1500 bp upstream of TSS (exclude the TSS200 window).
      # Build strand-aware by flanking 1500 bp then narrowing the TSS-proximal
      # 200 bp off. narrow() is NOT strand-aware, so split by strand:
      #   +strand: TSS is at higher coords → TSS-proximal end = end(range)
      #   -strand: TSS is at lower coords  → TSS-proximal end = start(range)
      # Preserves 1-to-1 correspondence with tss, so per-gene symbols remain
      # aligned with the output ranges.
      up_1500 <- GenomicRanges::flank(tss, 1500L)
      is_plus <- as.character(GenomicRanges::strand(up_1500)) %in% c("+", "*")
      out <- up_1500
      if (any(is_plus))
        out[is_plus]  <- GenomicRanges::narrow(up_1500[is_plus],  end   = -201L)
      if (any(!is_plus))
        out[!is_plus] <- GenomicRanges::narrow(up_1500[!is_plus], start =  201L)
      out
    },
    `1STEXON` = {
      # exon_rank is a property of a transcript, not of a gene: exonsBy(by =
      # "gene") does not carry it, so ranking exons there cannot work. Group by
      # transcript, take rank 1 (the first exon in transcription order, which on
      # the minus strand is the rightmost), then name by the parent gene.
      ex <- unlist(GenomicFeatures::exonsBy(txdb, by = "tx"))
      ex <- ex[GenomicRanges::mcols(ex)$exon_rank == 1L]
      .anno_name_by_parent_gene(txdb, ex)
    },
    `5UTR` = {
      .anno_name_by_parent_gene(
        txdb, unlist(GenomicFeatures::fiveUTRsByTranscript(txdb)))
    },
    `3UTR` = {
      .anno_name_by_parent_gene(
        txdb, unlist(GenomicFeatures::threeUTRsByTranscript(txdb)))
    },
    BODY = {
      all_genes
    },
    EXONBND = {
      exon_gr  <- unlist(GenomicFeatures::exonsBy(txdb, by = "gene"))
      starts   <- GenomicRanges::resize(exon_gr, 1L, fix = "start")
      ends     <- GenomicRanges::resize(exon_gr, 1L, fix = "end")
      bnds     <- c(starts, ends)
      bnds + 50L  # expand symmetrically by 50 bp on each side
    },
    WHOLE = {
      # Full gene span = BODY + TSS1500 upstream
      tss1500 <- GenomicRanges::setdiff(
        GenomicRanges::flank(tss, 1500L),
        GenomicRanges::flank(tss, 200L)
      )
      GenomicRanges::reduce(c(all_genes, tss1500))
    },
    stop("Unknown GENE subarea: '", subarea, "'. ",
         "Supported: TSS200, TSS1500, 1STEXON, 5UTR, 3UTR, BODY, EXONBND, WHOLE")
  )

  # For subareas that retain the per-gene structure, attach gene symbols
  if (subarea %in% c("TSS200", "TSS1500", "BODY")) {
    GenomicRanges::mcols(gr)$label <- symbols
  } else if (subarea == "WHOLE") {
    # After reduce(), per-gene identity is lost - use region coordinates as label
    GenomicRanges::mcols(gr)$label <-
      paste0(GenomicRanges::seqnames(gr), ":",
             GenomicRanges::start(gr), "-",
             GenomicRanges::end(gr))
  } else {
    # 1STEXON, 5UTR, 3UTR, EXONBND: label by parent gene via names
    gene_ids <- names(gr)
    if (!is.null(gene_ids) && length(gene_ids) == length(gr)) {
      syms <- .anno_entrez_to_symbol(gene_ids)
      GenomicRanges::mcols(gr)$label <- syms
    } else {
      GenomicRanges::mcols(gr)$label <- paste0("GENE_", subarea, "_",
                                                seq_along(gr))
    }
  }
  gr
}

.anno_build_island_area <- function(subarea, genome_build, islands = NULL) {
  # `islands` is injectable (a GRanges of island cores) so the WHOLE/ISLAND/
  # OPENSEA/shore/shelf logic can be unit-tested without AnnotationHub.
  if (is.null(islands)) islands <- .anno_get_cpg_islands(genome_build)

  # OPENSEA: inter-neighbourhood gaps, self-labelled by coordinate. Returned
  # early because its ranges (and count) do not correspond to island cores.
  # See island_opensea.R for the shared semantics with the Illumina backend.
  if (subarea == "OPENSEA")
    return(.anno_opensea_gaps(islands,
                         chrom_ends = as.list(GenomeInfoDb::seqlengths(islands))))

  # Unique island label: "seqname:start-end"
  island_labels <- paste0(
    GenomicRanges::seqnames(islands), ":",
    GenomicRanges::start(islands), "-",
    GenomicRanges::end(islands)
  )

  gr <- switch(subarea,
    # WHOLE = whole island neighbourhood (core + shores + shelves), mirroring
    # GENE_WHOLE = whole gene. ISLAND = the core alone (former WHOLE).
    WHOLE   = islands + .ISLAND_FLANK,
    ISLAND  = islands,
    N_SHORE = GenomicRanges::flank(islands, 2000L, start = TRUE),
    S_SHORE = GenomicRanges::flank(islands, 2000L, start = FALSE),
    N_SHELF = GenomicRanges::flank(
                GenomicRanges::flank(islands, 2000L, start = TRUE),
                2000L, start = TRUE),
    S_SHELF = GenomicRanges::flank(
                GenomicRanges::flank(islands, 2000L, start = FALSE),
                2000L, start = FALSE),
    stop("Unknown ISLAND subarea: '", subarea, "'. ",
         "Supported: WHOLE, ISLAND, N_SHORE, S_SHORE, N_SHELF, S_SHELF, OPENSEA")
  )
  GenomicRanges::mcols(gr)$label <- island_labels
  gr
}

.anno_build_chr_area <- function(subarea, genome_build, txdb) {
  if (subarea == "WHOLE") {
    si      <- GenomeInfoDb::seqinfo(txdb)
    lens    <- GenomeInfoDb::seqlengths(si)
    valid   <- !is.na(lens)
    gr <- GenomicRanges::GRanges(
      seqnames = names(lens)[valid],
      ranges   = IRanges::IRanges(1L, lens[valid])
    )
    GenomicRanges::mcols(gr)$label <- as.character(
      GenomicRanges::seqnames(gr))
    return(gr)
  }

  if (subarea == "CYTOBAND") {
    # Use bundled cytoband_hg19 data (available in data/)
    cb_obj <- tryCatch(
      get("cytoband_hg19", envir = asNamespace("SEMseeker")),
      error = function(e) NULL
    )
    if (is.null(cb_obj) || !is.data.frame(cb_obj))
      stop("cytoband_hg19 data object not found in SEMseeker package.")
    # cytoband_hg19$CHR is stored as factor with an empty level (""); convert
    # to character and drop rows with missing/empty seqnames - GenomicRanges
    # refuses GRanges construction if any seqlevel is NA or "".
    cb_chr <- as.character(cb_obj$CHR)
    valid  <- !is.na(cb_chr) & nzchar(cb_chr)
    cb_obj <- cb_obj[valid, , drop = FALSE]
    cb_chr <- cb_chr[valid]
    gr <- GenomicRanges::GRanges(
      seqnames = cb_chr,
      ranges   = IRanges::IRanges(cb_obj$START, cb_obj$END)
    )
    GenomicRanges::mcols(gr)$label <- as.character(cb_obj$CYTOBAND)
    return(gr)
  }

  stop("Unknown CHR subarea: '", subarea, "'. Supported: WHOLE, CYTOBAND")
}

.anno_build_dmr_area <- function(subarea) {
  # dmr_annotation maps Illumina probe IDs to DMR names (PROBE, DMR_WHOLE, DMR_DMR).
  # To get genomic coordinates we join with the K850 "Locations" table.
  dmr_obj <- tryCatch(
    get("dmr_annotation", envir = asNamespace("SEMseeker")),
    error = function(e) NULL
  )
  if (is.null(dmr_obj) || !is.data.frame(dmr_obj))
    stop("dmr_annotation data object not found in SEMseeker package.")

  label_col <- if (subarea %in% c("WHOLE", "DMR_WHOLE")) "DMR_WHOLE"
               else if (subarea == "DMR") "DMR_DMR"
               else colnames(dmr_obj)[ncol(dmr_obj)]

  # Load probe coordinates from K850 Locations table
  k850_pkg <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
  if (!requireNamespace(k850_pkg, quietly = TRUE))
    stop("DMR area coordinates require the K850 annotation package.\n",
         "Install: BiocManager::install('", k850_pkg, "')")

  locs <- .anno_pkg_load_table(k850_pkg, "Locations")
  locs$PROBE <- rownames(locs)

  # Join DMR labels with probe coordinates
  merged <- merge(
    dmr_obj[, c("PROBE", label_col)],
    locs[, c("PROBE", "chr", "pos")],
    by = "PROBE", all.x = FALSE
  )
  if (nrow(merged) == 0)
    stop("No probe coordinates found for DMR annotation. ",
         "Check that K850 annotation package is installed.")

  # Aggregate: DMR genomic extent = min-to-max of probe positions per DMR name
  merged$chr_clean <- sub("^chr", "", merged$chr)
  dmr_ranges <- do.call(rbind, lapply(
    split(merged, merged[[label_col]]),
    function(d) {
      data.frame(
        label = d[[label_col]][1],
        chr   = d$chr_clean[1],
        start = min(d$pos, na.rm = TRUE),
        end   = max(d$pos, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  ))

  gr <- GenomicRanges::GRanges(
    seqnames = paste0("chr", dmr_ranges$chr),
    ranges   = IRanges::IRanges(dmr_ranges$start, dmr_ranges$end)
  )
  GenomicRanges::mcols(gr)$label <- dmr_ranges$label
  gr
}

# ---------------------------------------------------------------------------
# Main exported function
# ---------------------------------------------------------------------------

#' Build a GRanges object for a given genomic area/subarea
#'
#' Constructs the same region boundaries used by Illumina array annotation
#' (TSS200, TSS1500, gene body, CpG islands, shores, shelves ...) from TxDb
#' packages and AnnotationHub, so that WGBS and long-read analyses share
#' identical region semantics with Illumina array analyses.
#'
#' Results are cached in memory for the duration of the R session to avoid
#' repeated package loads. CpG island tracks downloaded from AnnotationHub
#' are also cached on disk in \code{tools::R_user_dir("SEMseeker", "cache")}.
#'
#' @param area_subarea Character scalar: area and subarea joined by \code{"_"}
#'   (e.g. \code{"GENE_BODY"}, \code{"ISLAND_N_SHORE"}, \code{"CHR_CYTOBAND"}).
#' @param genome_build Character scalar: reference assembly.  Defaults to
#'   \code{ssEnv$genome_build} (set by \code{\link{core_init_env}}), or
#'   \code{"hg19"} if the session is not initialised.
#'
#' @return A \code{GRanges} object.  \code{mcols(gr)$label} contains the
#'   subarea identifier used downstream to group CpGs (gene symbol, island
#'   coordinate, cytoband name, etc.).
#'
#' @section PROBE_WHOLE semantics by technology:
#' \code{PROBE_WHOLE} is \strong{not} handled by this function.  It is resolved
#' inline by \code{\link{anno_probe_features_get}}:
#' \itemize{
#'   \item \strong{Illumina}: one row per array probe (manufacturer ID,
#'     e.g. \code{cg00000029}).  Probe identity is meaningful and cross-study
#'     comparable for the same array platform.
#'   \item \strong{WGBS / LONGREAD}: treated as \code{POSITION_WHOLE} - one
#'     row per genomic position encoded as \code{"CHR\_START"} (e.g.
#'     \code{"1\_10000"}).  Cross-study comparisons require the same
#'     \code{genome_build}.
#' }
#'
#' @section Required packages:
#' Install via \code{BiocManager::install()}.
#' \describe{
#'   \item{GENE areas}{\code{TxDb.Hsapiens.UCSC.hg19.knownGene} (or hg38/mm10),
#'     \code{GenomicFeatures}, \code{GenomicRanges}, \code{IRanges}.
#'     \code{org.Hs.eg.db} is required in practice: without it the label of
#'     a gene region stays the Entrez id that TxDb keys on, while the
#'     Illumina path labels with gene symbols, so the two backends name the
#'     same gene two different ways. The code falls back rather than
#'     failing, but the result is not comparable across backends.}
#'   \item{ISLAND areas}{\code{AnnotationHub} (downloads track on first use,
#'     then caches locally).}
#'   \item{CHR / DMR areas}{no extra packages needed.}
#' }
#'
#' @importFrom GenomicRanges GRanges
anno_area_granges_build <- function(area_subarea, genome_build = NULL) {

  if (is.null(genome_build)) {
    ssEnv        <- tryCatch(core_get_session_info(), error = function(e) list())
    genome_build <- ssEnv$genome_build %||% "hg19"
  }

  if (!grepl("_", area_subarea))
    area_subarea <- paste0(area_subarea, "_WHOLE")

  cache_key <- paste0(area_subarea, "__", genome_build)
  if (exists(cache_key, envir = .area_gr_cache, inherits = FALSE))
    return(get(cache_key, envir = .area_gr_cache, inherits = FALSE))

  # Check for GenomicRanges (needed for all paths)
  for (pkg in c("GenomicRanges", "IRanges", "S4Vectors")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for anno_area_granges_build().\n",
           "Install with: BiocManager::install(c('GenomicRanges','IRanges','S4Vectors'))")
  }

  parts   <- strsplit(area_subarea, "_", fixed = TRUE)[[1]]
  area    <- parts[1]
  subarea <- paste(parts[-1], collapse = "_")

  gr <- switch(area,
    GENE = {
      if (!requireNamespace("GenomicFeatures", quietly = TRUE))
        stop("GenomicFeatures required for GENE areas.\n",
             "Install: BiocManager::install('GenomicFeatures')")
      txdb <- .anno_get_txdb(genome_build)
      .anno_build_gene_area(subarea, txdb)
    },
    ISLAND = {
      .anno_build_island_area(subarea, genome_build)
    },
    CHR = {
      txdb <- .anno_get_txdb(genome_build)
      .anno_build_chr_area(subarea, genome_build, txdb)
    },
    DMR = {
      .anno_build_dmr_area(subarea)
    },
    stop("Unknown area '", area, "' in area_subarea = '", area_subarea, "'.\n",
         "Supported areas: GENE, ISLAND, CHR, DMR")
  )

  assign(cache_key, gr, envir = .area_gr_cache)
  gr
}

# Null-coalescing operator (safe to define here if not defined elsewhere)
`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !is.na(x[1]) && nzchar(x[1])) x else y
