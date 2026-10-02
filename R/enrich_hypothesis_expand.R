#' The two biological questions an enrichment has to answer first (internal)
#'
#' An enrichment over genes is not one analysis, it is a family of them,
#' and which member you ran is decided by two answers. Until this release the
#' package answered both on the caller's behalf and said so nowhere:
#' `util_keys_create()` rewrote every `HYPER` and `HYPO` figure to `HYPER_HYPO`
#' and every region window other than `WHOLE` to `ALL_SUBAREAS`, so no argument
#' of any exported function could reach an enrichment restricted to
#' hypermethylation, or to promoters.
#'
#' The two answers are now the caller's, they have no defaults, and they travel
#' with the result.
#'
#' **Why they are one question and not two.** The direction of an epimutation
#' maps onto transcription only together with where it sits: hypermethylation at
#' a promoter silences, hypomethylation there de-represses, but hypermethylation
#' in the body of a gene is associated with *active* transcription - the same
#' sign means the opposite thing. Naming a direction without naming a region is
#' not yet a biological hypothesis, which is why neither may be omitted.
#'
#' @name enrich_hypothesis
#' @keywords internal
NULL

#' Region windows a `gene_region` answer stands for (internal)
#'
#' The aliases are named after the biology and expand to the windows of the
#' Illumina gene annotation. An alias is admissible here - where the package
#' otherwise refuses to bury a choice in a constant - only because its expansion
#' does not stay hidden: it is written in the log when it is expanded and
#' stamped into the enrichment result, so a reader of that result can see which
#' windows produced it and disagree. A caller who reads `PROMOTER` differently
#' passes the windows explicitly instead.
#'
#' \itemize{
#'   \item `PROMOTER` - `PROMOTER`. One window, not three. The promoter is a
#'     region class of its own now, masked over `TSS200`, `TSS1500` and the first
#'     exon while the analysis runs, so a gene has one measure of its promoter and
#'     one p-value. Naming the three windows here instead put the gene in three
#'     families, and a gene with three p-values has none.
#'   \item `GENE_BODY` - `BODY`. The transcribed region, where the association
#'     between methylation and expression runs the other way.
#'   \item `WHOLE_GENE` - `WHOLE`. Every position of the gene, one window: the
#'     hypothesis is that the locus is affected, without saying where.
#' }
#'
#' All three stand for a single window, so none of them is a union any more. A
#' caller can still create one by passing several windows explicitly, and what
#' that costs is stated where it is chosen.
#'
#' @param gene_region one alias, or an explicit vector of windows.
#' @return character vector of `SUBAREA` values, and the alias as an attribute.
#' @keywords internal
#' @noRd
util_gene_region_expand <- function(gene_region) {

  if (missing(gene_region) || is.null(gene_region) || !length(gene_region) ||
      !all(nzchar(as.character(gene_region))))
    stop(
      "gene_region has no default, because the answer changes what the result ",
      "means:\n",
      "  \"PROMOTER\"   - the promoter region class. Methylation here gates\n",
      "                 transcription, so the enrichment speaks about\n",
      "                 transcriptional dysregulation.\n",
      "  \"GENE_BODY\"  - BODY. The association between methylation and\n",
      "                 expression runs the other way here.\n",
      "  \"WHOLE_GENE\" - every position of the gene. The hypothesis is only\n",
      "                 that the locus is affected, not where.\n",
      "Or pass the windows explicitly, e.g. c(\"TSS200\", \"TSS1500\").",
      call. = FALSE)

  # Each alias stands for one window. PROMOTER named three until the promoter
  # became a region class of its own: the three are masked together while the
  # analysis runs, so the grouping that used to happen here, on the way out, is
  # not needed and is not wanted - it reached the same gene from three families.
  aliases <- list(PROMOTER   = "PROMOTER",
                  GENE_BODY  = "BODY",
                  WHOLE_GENE = "WHOLE")

  asked <- core_name_cleaning(as.character(gene_region))

  if (length(asked) == 1L && asked %in% names(aliases)) {
    windows <- aliases[[asked]]
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " gene_region = ", asked, " stands for ",
              paste(windows, collapse = " + "),
              ". This is stamped on the result; pass the windows explicitly to ",
              "use a different definition.")
    return(structure(windows, alias = asked))
  }

  if (any(asked %in% names(aliases)))
    stop("gene_region mixes an alias with explicit windows: ",
         paste(asked, collapse = ", "),
         ". Pass one alias, or the windows it should stand for - not both, or ",
         "the result cannot say which definition produced it.", call. = FALSE)

  structure(asked, alias = paste(asked, collapse = "+"))
}

#' Figures an `epimutation_direction` answer stands for (internal)
#'
#' `ANY` is the union of the two directions, and it is a legitimate hypothesis -
#' for a marker of instability the direction is how the loss of epigenetic
#' control shows itself, not what defines it. What it costs is stated where it
#' is chosen: over a union a gene reaches the list from more than one family, so
#' **a p-value per gene is not defined**, and combining the ones it has would be
#' the ex-post combination this package refuses elsewhere. `ANY` therefore
#' yields a set of genes and no p-value column, and a backend that needs one
#' refuses it by name.
#'
#' @param epimutation_direction `"HYPER"`, `"HYPO"` or `"ANY"`.
#' @return character vector of `FIGURE` values, and whether it is a union.
#' @keywords internal
#' @noRd
util_epimutation_direction_expand <- function(epimutation_direction) {

  if (missing(epimutation_direction) || is.null(epimutation_direction) ||
      length(epimutation_direction) != 1L ||
      !nzchar(as.character(epimutation_direction)))
    stop(
      "epimutation_direction has no default, because the answer changes what ",
      "you are claiming:\n",
      "  \"HYPER\" - hypermethylation. At a promoter this silences: the\n",
      "            hypothesis is loss of expression.\n",
      "  \"HYPO\"  - hypomethylation. At a promoter this de-represses: the\n",
      "            hypothesis is gain.\n",
      "  \"ANY\"   - neither: the hypothesis is only that the locus lost\n",
      "            epigenetic control. A p-value per gene is NOT defined on\n",
      "            this answer, and backends that need one refuse it.\n",
      "The sign is read together with gene_region: HYPER at a promoter\n",
      "silences, HYPER in the body of a gene goes with active transcription.",
      call. = FALSE)

  asked <- core_name_cleaning(as.character(epimutation_direction))

  figures <- switch(asked,
                    HYPER = "HYPER",
                    HYPO  = "HYPO",
                    ANY   = c("HYPER", "HYPO"),
                    stop("epimutation_direction '", asked,
                         "' is not one of HYPER, HYPO, ANY.", call. = FALSE))

  structure(figures, alias = asked, union = length(figures) > 1L)
}
