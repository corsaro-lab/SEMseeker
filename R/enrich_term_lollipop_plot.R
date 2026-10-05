#' Lollipop plot of the enriched terms of a study
#'
#' One row per enriched term, the term's significance on the x axis, the size of
#' each point its fold enrichment, and its colour and shape a column of the
#' report chosen by the caller. One chart per category of term, because a GO
#' biological process and a KEGG pathway are not comparable on one axis.
#'
#' It is the second door of the same family as [enrich_circos_plot()]: both read
#' the report the enrichment run wrote, both take the enricher as a parameter and
#' resolve its column names through `ssEnv$key_enrichment_format`, so neither
#' knows anything about a particular enricher.
#'
#' @section Through 0.99.6 this was `enrich_lollipop_plot()`, and it could not run:
#' Measured on 2026-10-05, against a report of the shape the pipeline actually
#' writes, three defects, each sufficient on its own:
#'
#' * `performance_category` was used six times - for the colour, the shape, both
#'   legend titles, a factor conversion and a level count - and was never
#'   assigned, was not a parameter, and exists nowhere in the package. The call
#'   ended in `object 'performance_category' not found`. It is now
#'   `group_column`, a parameter, validated against the report.
#' * the `ggplot` object was never assigned. The whole chart was built as an
#'   expression statement inside a `for` loop, where a value is not printed, and
#'   then `ggsave()` was called **without** a `plot` argument, so it saved
#'   `last_plot()` - whatever had been displayed most recently, which inside that
#'   loop is nothing or a chart belonging to something else.
#' * the descriptions were starred with `data_to_plot$by_keyword`, a column
#'   produced only by `enrich_taxonomies_armonyser()`, which has no callers. The
#'   report the live pipeline writes has no such column, so after fixing the
#'   first defect this one assigned a zero-length vector to a column of n rows.
#'   The star is now drawn when the column is there and skipped when it is not,
#'   because it is an optional annotation from a chain that is currently dead.
#'
#' @section A p-value of exactly zero:
#' An enricher that reports `0` is reporting "below what I can represent", not
#' "impossible". On a `-log10` axis that is an infinity, which drops the point
#' and silently removes the most significant term from the chart. Zeroes are
#' replaced by the smallest positive value in the same report, so the term stays
#' on the chart at the edge of what that run could resolve. The substitution is
#' a floor and not a measurement: it says the term is at least this significant.
#'
#' @param enricher character. The enricher whose report to read, as its `label`
#'   appears in `ssEnv$key_enrichment_format`.
#' @param inference_details data.frame. One row per request. One chart is drawn
#'   per row, per key and per category of term.
#' @param group_column character. Column of the report that gives each point its
#'   colour and shape, e.g. `"MARKER"`, `"SS_CATEGORY"`, `"PHENOTYPE"`. A column
#'   the report does not have is refused, with the available ones named.
#' @param sort_column character. Column the terms are ordered by before the top
#'   ones are kept, e.g. `"SS_RANK"`.
#' @param top integer. How many distinct terms to keep.
#' @param statistic_parameter,pvalue_column,significance As in the enrichment
#'   run that produced the report, because they are part of the name of the file
#'   to read.
#' @param overwrite logical. Redraw when the file exists.
#'
#' @return Invisibly a character vector of the files written, empty when no
#'   report of that enricher was found.
#'
#' @examples
#' # The chart reads an enrichment report, so it needs a result folder where an
#' # enrichment has already run. The enricher and the column that groups the
#' # points are its two choices:
#' inference_details <- data.frame(
#'   independent_variable = "Sample_Group",
#'   family_test          = "wilcoxon",
#'   transformation_y     = "none",
#'   aggregation          = "MEAN",
#'   scope                = "INSTANCE"
#' )
#' \dontrun{
#' enrich_term_lollipop_plot(
#'   enricher          = "pathfindR",
#'   inference_details = inference_details,
#'   group_column      = "MARKER"
#' )
#' }
#' @export
enrich_term_lollipop_plot <- function(enricher, inference_details,
                                      group_column = "MARKER",
                                      sort_column = "SS_RANK",
                                      top = 50L,
                                      statistic_parameter = "",
                                      pvalue_column = "PVALUE_ADJ_ALL_BH",
                                      significance = TRUE,
                                      overwrite = FALSE) {

  ssEnv <- core_get_session_info()
  rules <- .enrich_rules_for(enricher, ssEnv)

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

      report <- utils::read.csv2(report_path, header = TRUE,
                                 stringsAsFactors = FALSE)
      canonical <- .enrich_report_canonical(report, rules)
      if (is.null(canonical))
        next

      chart_folder <- io_dir_check_and_create(ssEnv$result_folderChart,
                                              c("TERM_LOLLIPOP", enricher))

      for (category in unique(as.character(canonical$SS_CATEGORY))) {
        of_category <- canonical[
          as.character(canonical$SS_CATEGORY) == category, , drop = FALSE]
        if (nrow(of_category) == 0L)
          next

        chart_path <- io_file_path_build(
          chart_folder,
          c(analysis_name, core_name_cleaning(category)), ssEnv$plot_format)
        if (file.exists(chart_path) && !overwrite) {
          core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                    " chart already exists, skipped: ", chart_path)
          written <- c(written, chart_path)
          next
        }

        plot_object <- .enrich_lollipop_plot_build(
          of_category, group_column = group_column,
          sort_column = sort_column, top = top,
          palette = ssEnv$color_palette, alpha = as.numeric(ssEnv$alpha))

        ggplot2::ggsave(chart_path, plot = plot_object, width = 16, height = 9,
                        units = "in",
                        dpi = as.numeric(ssEnv$plot_resolution_ppi))
        written <- c(written, chart_path)
      }
    }
  }

  invisible(written)
}

# The report with the four columns the chart needs under fixed names, taken from
# the rules row of the enricher. The report on disk carries the enricher's own
# names: `highest_p` for pathfindR, `FDR` for WebGestalt, `Score` for
# phenolyzer. Renaming once here is what lets the drawing below know nothing
# about which enricher it is drawing.
.enrich_report_canonical <- function(report, rules) {

  wanted <- c(ID = as.character(rules$column_of_id),
              Description = as.character(rules$column_of_description),
              P_Value = as.character(rules$column_of_adj_pvalue),
              Enrichment = as.character(rules$column_of_enrichment))

  absent <- wanted[!(wanted %in% colnames(report))]
  if (length(absent) > 0L)
    stop("the ", as.character(rules$label), " report does not have the ",
         "column(s) its format declares: ",
         paste(sprintf("%s (for %s)", absent, names(absent)),
               collapse = ", "),
         ". Columns present: ", paste(colnames(report), collapse = ", "))

  for (canonical_name in names(wanted))
    colnames(report)[colnames(report) == wanted[[canonical_name]]] <-
      canonical_name

  if (nrow(report) == 0L)
    return(NULL)

  # A category is needed to split the charts. A report without one is one
  # category, which is true and keeps the chart drawable.
  if (!("SS_CATEGORY" %in% colnames(report)))
    report$SS_CATEGORY <- "ALL"

  report
}

# The drawing, separated from the door. It returns the object, which is the
# whole point: the chart this replaces was built and dropped, and ggsave() was
# left to guess which plot to write from whatever had last been displayed.
.enrich_lollipop_plot_build <- function(data, group_column, sort_column,
                                        top = 50L, palette,
                                        alpha = 0.05) {

  for (column in c(group_column, sort_column))
    if (!(column %in% colnames(data)))
      stop("'", column, "' is not a column of the report. Available: ",
           paste(colnames(data), collapse = ", "))

  data <- as.data.frame(data)
  data$P_Value <- as.numeric(data$P_Value)
  data$Enrichment <- as.numeric(data$Enrichment)

  # Zero means "below what this run can represent", and -log10(0) is an
  # infinity that drops the most significant term off the chart. The floor is
  # the smallest positive value the same report reached.
  positive <- data$P_Value[is.finite(data$P_Value) & data$P_Value > 0]
  if (length(positive) > 0L)
    data$P_Value[!is.na(data$P_Value) & data$P_Value == 0] <- min(positive)
  data$log_fdr <- -log10(data$P_Value)

  data <- data[order(data[[sort_column]], decreasing = FALSE), , drop = FALSE]
  kept <- utils::head(unique(as.character(data$Description)), top)
  data <- data[as.character(data$Description) %in% kept, , drop = FALSE]

  data[[group_column]] <- as.factor(data[[group_column]])

  # The star marks a term matched by a keyword. The column comes from a chain
  # that currently has no callers, so it is drawn when present and skipped when
  # absent rather than being assumed.
  if ("by_keyword" %in% colnames(data))
    data$Description <- ifelse(as.logical(data$by_keyword),
                               paste(data$Description, "*"),
                               as.character(data$Description))

  # Most significant at the top: the y axis is discrete, so the order is the
  # order of the levels, and ggplot2 draws the first level at the bottom.
  data$Description <- factor(as.character(data$Description),
                             levels = rev(unique(as.character(data$Description))))

  group_count <- length(levels(data[[group_column]]))

  ggplot2::ggplot(
      data,
      ggplot2::aes(x = .data$log_fdr, y = .data$Description,
                   size = .data$Enrichment,
                   colour = .data[[group_column]],
                   shape = .data[[group_column]])) +
    ggplot2::geom_point(fill = NA, stroke = 1.5) +
    ggplot2::scale_shape_manual(values = seq_len(group_count) + 20L) +
    ggplot2::scale_size_continuous(range = c(3, 10)) +
    ggplot2::scale_colour_manual(
      values = grDevices::colorRampPalette(palette)(group_count)) +
    ggplot2::geom_vline(xintercept = -log10(alpha), linetype = "dashed",
                        colour = palette[1], linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = -log10(alpha / 5), linetype = "dashed",
                        colour = palette[2], linewidth = 0.5) +
    ggplot2::labs(x = "-log10(adjusted p-value)", y = NULL,
                  colour = group_column, shape = group_column,
                  size = "Enrichment") +
    ggplot2::guides(
      shape = ggplot2::guide_legend(override.aes = list(size = 5))) +
    ggplot2::theme_minimal()
}
