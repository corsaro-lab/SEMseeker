#' Circos plot of an enrichment result: the genes of each term, and where they sit
#'
#' Draws one circular plot per key of the study: the genome as an ideogram, the
#' genes an enrichment implicated as labels around it, the enriched terms laid
#' out on a pseudo-chromosome, and a link from each gene to the term it belongs
#' to. Which enrichment is read is a parameter, the way it is for
#' [enrichment_analysis()], so one door covers every enricher whose report
#' declares where its genes are.
#'
#' @section What this replaces, and why one door and not two:
#' Through 0.99.6 the circos drawing existed twice.
#' `enrich_pathfindR_circlize()` read a pathfindR report, built the gene-to-term
#' pairs and drew; `plot_area_plot_circlize()` held the same thirteen drawing
#' lines and nothing else - it was an extraction of that drawing left half done,
#' and the frame it drew its links from, `results`, was never assigned and was
#' never a parameter, so it raised `object 'results' not found` on the line that
#' chose the link colours. There was no second chart, only a second copy of one.
#'
#' Both also drew the ideogram of the wrong genome whenever a study was not
#' hg19: they called `circlize::read.cytoband()` with no species, and the
#' session has recorded `genome_build` all along.
#'
#' @section The enricher is dispatched through the rules table:
#' `ssEnv$key_enrichment_format` already says, per enricher, which column of its
#' report holds the identifier, the description, the adjusted p-value and the
#' enrichment. It did not say which column holds the **genes**, which is the one
#' thing a circos needs, because nothing in the package had ever read a gene
#' list out of any report except pathfindR's. The table now carries
#' `column_of_genes`, `+`-separated when an enricher splits its genes across
#' more than one column, as pathfindR does with `Up_regulated` and
#' `Down_regulated`.
#'
#' Six of the seven enrichers have that field empty, and this function
#' **refuses** them by name rather than drawing an empty circle. Filling one in
#' is a single table entry, and it should be filled by someone looking at that
#' enricher's report rather than guessed from its documentation: a chart drawn
#' from a column that does not exist is a chart of nothing, and an empty circos
#' looks exactly like a study with no enrichment.
#'
#' @param enricher character. The enricher to read, as its `label` appears in
#'   `ssEnv$key_enrichment_format`, e.g. `"pathfindR"`.
#' @param inference_details data.frame. One row per request, the same shape
#'   [association_analysis()] takes. One chart is drawn per row and per key.
#' @param terms_selection character vector of term identifiers to keep, or
#'   `NULL` for all of them. A whole-genome circos of every enriched term is
#'   unreadable, so naming the terms is the normal case rather than the
#'   exception.
#' @param statistic_parameter character. As in the enrichment run that produced
#'   the report, because it is part of the name of the file to read.
#' @param pvalue_column character. Likewise.
#' @param significance logical. Likewise.
#' @param overwrite logical. Redraw when the file exists.
#'
#' @return Invisibly a character vector of the files written, which is empty
#'   when no report of that enricher was found.
#'
#' @examples
#' # The chart reads an enrichment report, so it needs a result folder where an
#' # enrichment has already run. Its first argument is the enricher, named as
#' # the rules table names it:
#' inference_details <- data.frame(
#'   independent_variable = "Sample_Group",
#'   family_test          = "wilcoxon",
#'   transformation_y     = "none",
#'   aggregation          = "MEAN",
#'   scope                = "INSTANCE"
#' )
#' \dontrun{
#' enrich_circos_plot(
#'   enricher          = "pathfindR",
#'   inference_details = inference_details,
#'   terms_selection   = c("hsa04110", "hsa04151")
#' )
#' }
#' @export
enrich_circos_plot <- function(enricher, inference_details,
                               terms_selection = NULL,
                               statistic_parameter = "",
                               pvalue_column = "PVALUE_ADJ_ALL_BH",
                               significance = TRUE,
                               overwrite = FALSE) {

  ssEnv <- core_get_session_info()
  rules <- .enrich_rules_for(enricher, ssEnv)
  # Checked here, before any file is read: an enricher whose genes are not
  # declared cannot produce a link, and the refusal belongs at the door.
  .enrich_gene_columns_of(rules, ssEnv)

  keys <- unique(ssEnv$keys_for_pathway)
  written <- character(0)

  for (request_index in seq_len(nrow(inference_details))) {
    inference_detail <- inference_details[request_index, ]

    for (key_index in seq_len(nrow(keys))) {
      key <- keys[key_index, ]
      suffix <- if (identical(statistic_parameter, "")) "without_signal_" else ""
      analysis_name <- enrich_phenotype_analysis_name(
        inference_detail, key, prefix = "", suffix = suffix,
        pvalue_column, ssEnv$alpha, significance)

      report_path <- io_file_path_build(
        io_enrichment_folder(inference_detail, enricher), analysis_name, "csv")
      if (!file.exists(report_path))
        next

      chart_folder <- io_dir_check_and_create(ssEnv$result_folderChart,
                                              c("CIRCOS", enricher))
      chart_path <- io_file_path_build(chart_folder, analysis_name,
                                       ssEnv$plot_format)
      if (file.exists(chart_path) && !overwrite) {
        core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                  " chart already exists, skipped: ", chart_path)
        written <- c(written, chart_path)
        next
      }

      report <- utils::read.csv2(report_path, header = TRUE,
                                 stringsAsFactors = FALSE)
      circos_data <- .enrich_circos_data_build(report, rules, terms_selection,
                                               ssEnv$tech)
      if (is.null(circos_data)) {
        core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                  " no gene-to-term pair survived for ", enricher, " / ",
                  analysis_name, ": nothing drawn.")
        next
      }

      .plot_circos_links_build(
        labels = circos_data$labels,
        links_from = circos_data$links_from,
        links_to = circos_data$links_to,
        genome_build = ssEnv$genome_build,
        plot_path = chart_path, plot_format = ssEnv$plot_format,
        resolution_ppi = ssEnv$plot_resolution_ppi,
        extra_sector = .ENRICH_CIRCOS_TERM_SECTOR)
      written <- c(written, chart_path)
    }
  }

  invisible(written)
}

# The pseudo-chromosome the enriched terms are laid out on. A name rather than a
# literal in three places, because it has to be the same string in the cytoband
# and in the links or the links land on a sector that does not exist.
.ENRICH_CIRCOS_TERM_SECTOR <- "term"

# The rules row of one enricher, or a refusal naming what is available. Reading
# the enricher from the table is what makes one function serve all of them.
#
# This is the lookup ONLY. The requirement that an enricher declare its gene
# column belongs to the chart that draws links and not to every chart that reads
# a report: a lollipop of enriched terms needs no genes, so checking here would
# refuse six enrichers a chart that can draw all seven.
.enrich_rules_for <- function(enricher, ssEnv) {

  formats <- as.data.frame(ssEnv$key_enrichment_format)
  row_index <- which(tolower(as.character(formats$label)) == tolower(enricher))

  if (length(row_index) == 0L)
    stop("unknown enricher '", enricher, "'. Known: ",
         paste(formats$label, collapse = ", "))

  formats[row_index[1], ]
}

# The gene columns of an enricher, or a refusal. Only a chart that links genes
# to terms needs this, which is why it is not part of the lookup above.
.enrich_gene_columns_of <- function(rules, ssEnv) {

  gene_columns <- trimws(as.character(rules$column_of_genes))

  if (length(gene_columns) == 0L || is.na(gene_columns) ||
      !nzchar(gene_columns)) {
    formats <- as.data.frame(ssEnv$key_enrichment_format)
    declared <- formats$label[nzchar(as.character(formats$column_of_genes))]
    stop("enricher '", as.character(rules$label), "' does not declare which ",
         "column of its report holds the genes, so there is nothing to link. ",
         "Add column_of_genes for it in the enrichment format table. ",
         "Enrichers that declare it: ", paste(declared, collapse = ", "))
  }

  util_split_and_clean(gene_columns)
}
