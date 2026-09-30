sem_cluster_analysis <- function(cluster_variables,ellipsis=TRUE, sql_sample_selection="", result_folder, maxResources = 90,
  parallel_strategy  = "multisession",start_fresh = FALSE, ...)
{
  # AI-061+ (2026-06-09): propagate start_fresh to core_init_env (was previously
  # hardcoded to FALSE). Same name → same semantics: TRUE wipes the whole
  # result_folder (Chart, Enrichment, Phenotype, Pivots, Data, Log, Inference);
  # FALSE preserves everything. The redundant unlink of result_folderInference
  # below is kept as a no-op when start_fresh=TRUE (core_init_env already wiped)
  # and as the SOLE reset when start_fresh=FALSE-but-caller-wants-inference-only
  # - but that latter case never happens with this signature.
  ssEnv <- core_init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = start_fresh, ...)

  core_log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), " SEMseeker will perform the cluster analysys for project \n in ", ssEnv$result_folderData)

  if(start_fresh)
    unlink(ssEnv$result_folderInference, recursive = TRUE)
  io_dir_check_and_create(ssEnv$result_folderInference,c())
  localKeys <- ssEnv$keys_areas_subareas_markers_figures

  localKeys <- localKeys[localKeys$AREA != "POSITION", ]

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = seq_len(nrow(localKeys)*length(cluster_variables)))
  else
    progress_bar <- ""


  for (i in seq_along(cluster_variables))
  {
    cluster_variable_name <- cluster_variables[i]
    if (is.na(cluster_variable_name) || cluster_variable_name=="")
    {
      core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " sem_cluster_analysis: cluster_variable is empty")
      stop()
    }

    study_summary <-   sem_study_summary_get(sql_sample_selection)
    if (!(cluster_variable_name %in% colnames(study_summary)))
    {
      core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " sem_cluster_analysis: cluster_variable is not in study_summary")
      stop()
    }

    for ( k in seq_len(nrow(localKeys)))
    {
      key <- localKeys[k,]
      chart_folder <- io_dir_check_and_create(ssEnv$result_folderChart, c("CLUSTER_ANALYSIS", core_name_cleaning(sql_sample_selection)))
      plot_filename <- io_file_path_build(baseFolder = chart_folder,
        detailsFilename = c(key$MARKER, key$FIGURE, key$AREA, key$SUBAREA, cluster_variable_name), extension = ssEnv$plot_format )

      if(file.exists(plot_filename))
        next

      # AI-027: read via unified dispatcher.
      pivot_lazy <- io_read_pivot(key$MARKER, key$FIGURE, key$AREA, key$SUBAREA)
      if (is.null(pivot_lazy))
        next
      pivot_data <- as.data.frame(pivot_lazy$collect())
      pivot_data <- pivot_data[, study_summary$Sample_ID]
      pivot_data[is.na(pivot_data)] <- 0

      tsne <- M3C::tsne(pivot_data,labels=as.factor( study_summary[, cluster_variable_name]))
      # get the plot
      plot <- ggplot2::ggplot(
        tsne$data,
        ggplot2::aes(x = .data$X1, y = .data$X2, color = as.factor(study_summary[, cluster_variable_name]))
      ) +
        ggplot2::geom_point() +
        ggplot2::labs(
          title = paste("t-SNE plot for", key$MARKER, key$FIGURE, key$AREA, key$SUBAREA),
          x = "t-SNE 1",
          y = "t-SNE 2"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::scale_color_discrete(name = cluster_variable_name)

      if(ellipsis)
        plot <- plot + ggplot2::stat_ellipse(ggplot2::aes(group = as.factor(study_summary[, cluster_variable_name])))


      # Save the plot
      plot_filename <- io_file_path_build(baseFolder = chart_folder,
        detailsFilename = c(key$MARKER, key$FIGURE, key$AREA, key$SUBAREA, cluster_variable_name), extension = "png")
      ggplot2::ggsave(plot_filename, plot = plot, width = 8, height = 6, dpi = ssEnv$plot_resolution_ppi, units = "in")

      if(ssEnv$showprogress)
        progress_bar(sprintf("Cluster analysis of area: %s, subarea: %s, marker: %s and figure: %s",key$AREA, key$SUBAREA, key$MARKER, key$FIGURE ))
    }
  }
  core_log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), " SEMseeker finished the cluster analysis for project \n in ", ssEnv$result_folderData)
  return(invisible())
}
