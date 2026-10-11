# ---------------------------------------------------------------------------
# .plot_comparison_pvalue_label()
#
# Formatted overall p-value for the box/violin group comparison. Replaces
# ggpubr::stat_compare_means(label = "p.format"), whose internal
# after_stat(create_p_label()) expression fails under ggplot2 >= 4.0 when
# ggpubr is loaded via :: but not attached - the CI ERROR observed in
# test-plot_box_plot.R ("could not find function 'create_p_label'").
#
# The p-value is computed directly from the stats package, so the result no
# longer depends on ggpubr/ggplot2 internal evaluation scoping. Returns
# NA_character_ for families that were not previously annotated, or when the
# test cannot be computed.
# ---------------------------------------------------------------------------
.plot_comparison_pvalue_label <- function(dataFrameToPlot, independent_variable,
                                          dependent_variable, family_test) {
  fml <- stats::as.formula(sprintf("`%s` ~ `%s`", dependent_variable, independent_variable))
  n_groups <- length(unique(stats::na.omit(dataFrameToPlot[[independent_variable]])))

  pval <- tryCatch({
    if (family_test == "t.test" && n_groups == 2) {
      stats::t.test(fml, data = dataFrameToPlot)$p.value
    } else if (family_test == "wilcox.test" && n_groups == 2) {
      stats::wilcox.test(fml, data = dataFrameToPlot)$p.value
    } else if (family_test == "kruskal.test") {
      stats::kruskal.test(fml, data = dataFrameToPlot)$p.value
    } else if (family_test == "anova") {
      stats::anova(stats::lm(fml, data = dataFrameToPlot))[["Pr(>F)"]][1]
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)

  if (is.null(pval) || length(pval) != 1 || is.na(pval)) return(NA_character_)
  paste0("p = ", format.pval(pval, digits = 2, eps = 1e-3))
}

plot_box_plot <- function (dataFrameToPlot, independent_variable,dependent_variable, transformation_y, family_test, samples_sql_condition="",key)
{
  if (!assoc_is_family_dicotomic(family_test))
    return()

  area <- as.character(key$AREA)
  subarea <- as.character(key$SUBAREA)
  marker <- as.character(key$MARKER)
  figure <- as.character(key$FIGURE)

  ssEnv <- core_get_session_info()
  chartFolder <- io_dir_check_and_create(ssEnv$result_folderChart,c("COMPARISON",core_name_cleaning(as.character(samples_sql_condition))))
  filename  <-  io_file_path_build(chartFolder,toupper(c("BOX_PLOT",family_test,as.character(transformation_y), independent_variable,"Vs", dependent_variable,area, subarea, marker, figure)),ssEnv$plot_format)
  if(!file.exists(filename))
  {
    num_boxplots <- length(unique(dataFrameToPlot[, independent_variable]))

    # Overall comparison p-value, computed once and reused for both plots.
    p_value_label <- .plot_comparison_pvalue_label(dataFrameToPlot, independent_variable,
      dependent_variable, family_test)

    p <- ggplot2::ggplot(dataFrameToPlot, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
      ggplot2::geom_boxplot(ggplot2::aes(fill = !!ggplot2::sym(independent_variable)), outlier.size = 1.5) +
      ggplot2::coord_cartesian(ylim = c(min(dataFrameToPlot[[dependent_variable]], na.rm=TRUE)*0.95,
        max(dataFrameToPlot[[dependent_variable]], na.rm=TRUE)*1.05)) +  # Zoom in to remove extremes
      ggplot2::theme_bw(base_size = 15) +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
      )

    if (num_boxplots <= length(ssEnv$color_palette)) {
      p <- p + ggplot2::scale_fill_manual(values = ssEnv$color_palette[seq_len(num_boxplots)])
    }

    ## add p-value annotation (manual; replaces ggpubr::stat_compare_means)
    if (!is.na(p_value_label)) {
      p <- p + ggplot2::annotate("text", x = (num_boxplots + 1) / 2, y = 0.5,
        label = p_value_label, size = 3)
    }

    ggplot2::ggsave(filename = filename, plot = p, width = 7, height = 7, dpi = as.numeric(ssEnv$plot_resolution_ppi), units = "in")


    #### VIOLIN PLOT #####
    filename  <-  io_file_path_build(chartFolder,
      toupper(c("VIOLIN_PLOT",family_test,as.character(transformation_y), independent_variable,"Vs", dependent_variable,area, subarea, marker, figure)),ssEnv$plot_format)

    p <- ggplot2::ggplot(dataFrameToPlot, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
      ggplot2::geom_violin(ggplot2::aes(fill = !!ggplot2::sym(independent_variable))) +
      ggplot2::stat_summary(fun.data = ggplot2::mean_sdl, fun.args = list(mult = 1), geom = "crossbar", width = 0.2, color = "black") +
      ggplot2::stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
      ggplot2::coord_cartesian(ylim = c(min(dataFrameToPlot[[dependent_variable]], na.rm=TRUE)*0.95,
        max(dataFrameToPlot[[dependent_variable]], na.rm=TRUE)*1.05)) +  # Zoom in to remove extremes
      ggplot2::theme_bw(base_size = 15) +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
      )

    if (num_boxplots <= length(ssEnv$color_palette)) {
      p <- p + ggplot2::scale_fill_manual(values = ssEnv$color_palette[seq_len(num_boxplots)])
    }

    ## add p-value annotation (manual; replaces ggpubr::stat_compare_means)
    if (!is.na(p_value_label)) {
      p <- p + ggplot2::annotate("text", x = (num_boxplots + 1) / 2, y = 0.5,
        label = p_value_label, size = 3)
    }

    ggplot2::ggsave(filename = filename, plot = p, width = 7, height = 7, dpi = as.numeric(ssEnv$plot_resolution_ppi), units = "in")

  }
}
