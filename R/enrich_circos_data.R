#' Gene-to-term pairs out of an enrichment report
#'
#' Reads the columns the rules row names and returns one row per (term, gene)
#' pair. Pure: it takes the report as a data.frame and knows nothing about where
#' it came from, which is what lets the pairs be asserted on directly.
#'
#' @section Why a declared column that is absent is a refusal:
#' A report whose gene column is missing produces no pairs, and no pairs draws a
#' circle with an ideogram and nothing on it - which is indistinguishable from a
#' study where the enrichment found nothing. The two have to be told apart, so a
#' missing column stops here and names itself.
#'
#' @param report data.frame. The enrichment report as read from disk.
#' @param rules one row of `ssEnv$key_enrichment_format`.
#' @param terms_selection character vector of term identifiers to keep, or
#'   `NULL` for all.
#' @param genes_separator character. Separator inside a gene column.
#'
#' @return data.frame with `TERM`, `TERM_LABEL` and `GENE`, or `NULL` when the
#'   selection leaves nothing.
#'
#' @keywords internal
#' @noRd
.enrich_circos_pairs_build <- function(report, rules, terms_selection = NULL,
                                       genes_separator = ",") {

  column_of_id <- as.character(rules$column_of_id)
  column_of_description <- as.character(rules$column_of_description)
  gene_columns <- util_split_and_clean(as.character(rules$column_of_genes))

  declared <- c(column_of_id, column_of_description, gene_columns)
  absent <- setdiff(declared, colnames(report))
  if (length(absent) > 0L)
    stop("the ", as.character(rules$label), " report does not have the ",
         "column(s) its format declares: ", paste(absent, collapse = ", "),
         ". Columns present: ", paste(colnames(report), collapse = ", "))

  if (!is.null(terms_selection))
    report <- report[as.character(report[[column_of_id]]) %in%
                       as.character(terms_selection), , drop = FALSE]
  if (nrow(report) == 0L)
    return(NULL)

  pairs <- do.call(rbind, lapply(seq_len(nrow(report)), function(i) {
    genes <- unlist(lapply(gene_columns, function(column) {
      cell <- report[i, column]
      if (length(cell) == 0L || is.na(cell) || !nzchar(trimws(cell)))
        return(character(0))
      trimws(strsplit(gsub(" ", "", as.character(cell)),
                      genes_separator, fixed = TRUE)[[1]])
    }))
    genes <- unique(genes[nzchar(genes)])
    if (length(genes) == 0L)
      return(NULL)
    data.frame(TERM = as.character(report[i, column_of_id]),
               TERM_LABEL = as.character(report[i, column_of_description]),
               GENE = genes, stringsAsFactors = FALSE)
  }))

  if (is.null(pairs) || nrow(pairs) == 0L) NULL else unique(pairs)
}

#' Lay the terms on a pseudo-chromosome and turn the pairs into links
#'
#' Each term gets an equal slice of the pseudo-sector, in the order the terms
#' first appear, and each pair becomes a link from the gene's interval to its
#' term's slice. Pure, so the arithmetic of the layout can be asserted without
#' a device.
#'
#' @section Equal slices, and the alternative that was rejected:
#' Slices could have been sized by the enrichment or by the gene count, which
#' would put more information in the picture. They are equal because a circos
#' already encodes a quantity in the position of every gene, and a second
#' quantity encoded as arc length reads as a genomic extent that does not exist.
#' A term is not a region.
#'
#' @param pairs data.frame from [.enrich_circos_pairs_build()].
#' @param gene_coordinates data.frame with `GENE`, `CHR`, `START`, `END`: one
#'   row per gene. A gene with no coordinates is dropped and counted, because a
#'   link needs two ends.
#' @param sector_name character. Name of the pseudo-chromosome.
#' @param sector_width numeric. Its length.
#'
#' @return A list of `labels`, `links_from` and `links_to`, or `NULL` when no
#'   pair has coordinates on both ends.
#'
#' @keywords internal
#' @noRd
.enrich_circos_layout_build <- function(pairs, gene_coordinates,
                                        sector_name = "term",
                                        sector_width = 2e8) {

  located <- merge(pairs, gene_coordinates, by = "GENE")
  if (nrow(located) == 0L)
    return(NULL)

  dropped <- length(setdiff(unique(pairs$GENE), unique(located$GENE)))
  if (dropped > 0L)
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " ", dropped, " of ", length(unique(pairs$GENE)),
              " genes have no coordinates in the annotation and carry no link.")

  # The order of the terms around the sector comes from the REPORT, not from the
  # merge above: a merge reorders, and an enrichment report is normally sorted
  # by significance, so taking its order makes the picture reproducible and puts
  # the strongest terms together instead of scattering them arbitrarily.
  terms <- unique(pairs$TERM)
  terms <- terms[terms %in% located$TERM]

  # The slice boundaries are computed as boundaries rather than as a start plus
  # a width: starting at 1 and adding the width n times ends at sector_width + 1,
  # one past the pseudo-chromosome, and circlize answers that with a warning
  # about regions past the end of the chromosome rather than an error.
  boundaries <- seq(from = 1, to = sector_width, length.out = length(terms) + 1L)
  term_layout <- data.frame(
    TERM = terms,
    TERM_START = boundaries[seq_len(length(terms))],
    TERM_END = boundaries[seq_len(length(terms)) + 1L],
    stringsAsFactors = FALSE)
  term_layout$TERM_LABEL <- located$TERM_LABEL[match(terms, located$TERM)]

  located <- merge(located, term_layout[, c("TERM", "TERM_START", "TERM_END")],
                   by = "TERM")

  gene_labels <- unique(located[, c("CHR", "START", "END", "GENE")])
  colnames(gene_labels) <- c("CHR", "START", "END", "LABEL")
  term_labels <- data.frame(CHR = sector_name,
                            START = term_layout$TERM_START,
                            END = term_layout$TERM_END,
                            LABEL = term_layout$TERM_LABEL,
                            stringsAsFactors = FALSE)

  list(
    labels = unique(rbind(gene_labels, term_labels)),
    links_from = data.frame(CHR = located$CHR, START = located$START,
                            END = located$END, stringsAsFactors = FALSE),
    links_to = data.frame(CHR = sector_name, START = located$TERM_START,
                          END = located$TERM_END, stringsAsFactors = FALSE))
}

#' From a report to everything the circos drawing needs
#'
#' Chains the three steps that have to happen in this order: the pairs out of
#' the report, the gene symbols normalised so they can be matched against the
#' annotation, and the coordinates of each gene. Each step is its own function
#' so that the one that can be wrong about the data is testable without a
#' device and without an annotation package.
#'
#' @param report data.frame. The enrichment report.
#' @param rules one row of `ssEnv$key_enrichment_format`.
#' @param terms_selection character vector or `NULL`.
#' @param tech character. Array technology, passed to the annotation builder.
#' @param sector_name,sector_width passed to [.enrich_circos_layout_build()].
#'
#' @return The list [.enrich_circos_layout_build()] returns, or `NULL`.
#'
#' @keywords internal
#' @noRd
.enrich_circos_data_build <- function(report, rules, terms_selection, tech,
                                      sector_name = "term",
                                      sector_width = 2e8) {

  pairs <- .enrich_circos_pairs_build(report, rules, terms_selection)
  if (is.null(pairs))
    return(NULL)

  pairs$GENE <- .enrich_official_symbols(pairs$GENE)
  pairs <- unique(pairs[!is.na(pairs$GENE), , drop = FALSE])
  if (nrow(pairs) == 0L)
    return(NULL)

  coordinates <- .enrich_gene_coordinates(unique(pairs$GENE), tech)
  if (is.null(coordinates))
    return(NULL)

  .enrich_circos_layout_build(pairs, coordinates, sector_name, sector_width)
}

#' Aliases to official symbols
#'
#' An enrichment report names genes however its source named them, and the probe
#' annotation names them officially. Without this step a gene whose report entry
#' is an alias has no coordinates and silently loses its link, which looks like a
#' gene the array does not cover.
#'
#' When the annotation packages are not installed the names are returned
#' unchanged rather than the chart refusing: the aliases are a minority, the
#' packages are suggested, and a chart of the genes that did match is still a
#' true chart of those genes. The count that could not be normalised is logged.
#'
#' @param genes character vector.
#' @return character vector of the same length, `NA` where nothing matched.
#'
#' @keywords internal
#' @noRd
.enrich_official_symbols <- function(genes) {

  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " AnnotationDbi or org.Hs.eg.db is not installed: gene aliases",
              " are not normalised to official symbols, so a gene named by an",
              " alias in the report will carry no link.")
    return(genes)
  }

  entrez <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db, keys = as.character(genes),
    column = "ENTREZID", keytype = "ALIAS", multiVals = "first"))

  resolvable <- !is.na(entrez)
  official <- as.character(genes)
  if (any(resolvable)) {
    symbols <- suppressMessages(AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db, keys = as.character(entrez[resolvable]),
      column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
    official[resolvable] <- as.character(symbols)
  }
  if (any(!resolvable))
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " ", sum(!resolvable), " of ", length(genes),
              " gene names did not resolve to an official symbol and are kept",
              " as the report wrote them.")
  official
}

#' The interval a gene occupies, from the probe annotation
#'
#' The widest span of the probes annotated to that gene: the smallest start and
#' the largest end. A gene is not one position and the annotation does not carry
#' gene bounds, so this is the extent the array actually interrogates, which is
#' also the honest thing for a chart drawn from array data to show.
#'
#' The `chr` prefix is added here. The annotation stores the chromosome without
#' it and the cytoband expects it; getting that wrong produces no error, only a
#' circos with every link missing.
#'
#' @param genes character vector of official symbols.
#' @param tech character. Array technology.
#' @return data.frame with `GENE`, `CHR`, `START`, `END`, or `NULL`.
#'
#' @keywords internal
#' @noRd
.enrich_gene_coordinates <- function(genes, tech) {

  probe_annotation <- anno_probe_annotation_build(tech)
  if (is.null(probe_annotation) || nrow(probe_annotation) == 0L)
    return(NULL)

  on_gene <- probe_annotation[
    as.character(probe_annotation$GENE_WHOLE) %in% as.character(genes),
    c("CHR", "START", "END", "GENE_WHOLE"), drop = FALSE]
  if (nrow(on_gene) == 0L)
    return(NULL)

  colnames(on_gene)[colnames(on_gene) == "GENE_WHOLE"] <- "GENE"
  starts <- stats::aggregate(START ~ GENE, data = on_gene, FUN = min)
  ends   <- stats::aggregate(END ~ GENE, data = on_gene, FUN = max)
  chromosome <- unique(on_gene[, c("GENE", "CHR")])
  # A gene appearing on more than one chromosome in the annotation cannot be one
  # interval; the first is taken and the rest counted, rather than producing two
  # links for one gene.
  duplicated_gene <- duplicated(chromosome$GENE)
  if (any(duplicated_gene))
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " ", sum(duplicated_gene), " gene(s) are annotated on more than",
              " one chromosome; the first is used.")
  chromosome <- chromosome[!duplicated_gene, , drop = FALSE]

  coordinates <- merge(merge(chromosome, starts, by = "GENE"), ends,
                       by = "GENE")
  coordinates$CHR <- paste0("chr", gsub("^chr", "", as.character(coordinates$CHR)))
  coordinates[, c("GENE", "CHR", "START", "END")]
}
