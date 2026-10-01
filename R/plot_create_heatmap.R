#' plot_create_heatmap load the multiple bed resulting from
#' analysis organized into files and folders per marker and produce a pivot
#'
#' @return nothing
#' @importFrom doRNG %dorng%
plot_create_heatmap <- function() {

  ssEnv <- core_get_session_info()
  sample_group_comb <- utils::combn(ssEnv$keys_sample_groups, 2)
  for (g in seq_len(sample_group_comb))
  {
    localKeys <- reshape::expand.grid.df(as.data.frame(ssEnv$keys_areas_subareas_markers_figures), sample_group_comb )
    variables_to_export <- c("markers", "annotatedData", "sample_group_comb", "chartFolder", "figures",
      "%dorng%", "j", "iter", "RNGseed", "checkRNGversion", "getRNG", "%||%", ".getDoParName",
      "getDoParName", "getDoBackend", "setDoBackend", "RNGtype", "showRNG", "doRNGversion",
      ".getRNG", ".getRNGattribute", "hasRNG", "isNumber", "isReal", "isInteger",
      "nextRNG", ".foreachGlobals", "RNGkind", "setRNG", "RNGprovider", ".RNGkind_length",
      "tail", "RNGstr", "localKeys","doRNGseq", "%dopar%", "getDoPar","variables_to_export_nested","get_probes")

    variables_to_export_nested <- c(variables_to_export, "markers", "annotatedData", "sample_group_comb", "chartFolder","g","figures","variables_to_export")
    # g <- 0
    i <- 0
    j <- 0

    foreach::foreach(j = seq_len(nrow(localKeys)), .export = variables_to_export_nested) %dorng%
    # for(j in 1:nrow(localKeys))
    {
      core_update_session_info(ssEnv)
      area <- localKeys[j,"AREA"]
      subarea <- localKeys[j,"SUBAREA"]
      marker <- localKeys[j,"MARKER"]
      figure <- localKeys[j,"FIGURE"]
      chartFolder <- io_dir_check_and_create(ssEnv$result_folderChart, area)

      annotatedData <-  read_annotated_bed(figure,marker,area,subarea)
      annotatedData <- subset(annotatedData, annotatedData$SAMPLE_GROUP %in% sample_group_comb)
      annotatedData <- subset(annotatedData, annotatedData$VALUE != 0 )
      annotatedData$SAMPLEID <- as.factor(annotatedData$SAMPLEID)
      annotatedData$FIGURE <- as.factor(annotatedData$FIGURE)
      annotatedData$MARKER <- as.factor(annotatedData$MARKER)
      annotatedData$SAMPLE_GROUP <- as.factor(annotatedData$SAMPLE_GROUP)
      annotatedData$VALUE <- as.numeric(annotatedData$VALUE)
      annotatedData$SUBAREA <- subarea

      if (is.null(annotatedData) | plyr::empty(annotatedData) | nrow(annotatedData)==0 || nrow(annotatedData) < 10)
        next

      if(!plyr::empty(annotatedData))
        if(nrow(annotatedData)>2)
        {
          annotatedData <- reshape2::dcast(data = annotatedData, SAMPLEID + SAMPLE_GROUP ~ SUBAREA, value.var = "VALUE", sum)
          row.names(annotatedData) <- annotatedData$SAMPLEID

          mainTitle <- paste0( paste0( sample_group_comb, collapse ="_Vs_")," ",marker, sep="")
          if(ncol(annotatedData)>1000)
          {
            annotatedData <- annotatedData[,seq_len(min(ncol(annotatedData),1005))]
            mainTitle <- paste0( mainTitle," (first 1000)",  sep="")
          }

          # col<- colorRampPalette(c("violet","white",ssEnv$color_palette[1]))(1024)
          # skip heatmap if no enough data are available
          tt <- as.data.frame(annotatedData[,3:ncol(annotatedData)])
          if (!plyr::empty(tt))
            if(nrow(tt) > 2 & ncol(tt) > 2)
            {

              # Prepare the filename
              filename <- paste0(chartFolder, "/", paste0(sample_group_comb, collapse = "_Vs_"), "_", marker, "_", figure, ".",ssEnv$plot_format)

              # Convert the data to a long format for ggplot2
              heatmap_data <- reshape2::melt(annotatedData, id.vars = "SAMPLE_GROUP")

              # Create the ggplot heatmap
              p <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = .data$variable, y = .data$SAMPLE_GROUP, fill = .data$value)) +
                ggplot2::geom_tile() +
                ggplot2::scale_fill_gradientn(colours = grDevices::cm.colors(256)) +
                ggplot2::theme_minimal(base_size = 15) +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1), plot.margin = unit(c(2.5, 2.5, 2.5, 2.5), "cm")) +
                ggplot2::labs(title = mainTitle, x = "", y = "")

              # Save the plot
              ggplot2::ggsave(filename = filename, plot = p, width = 2480/ssEnv$plot_resolution, height = 2480/ssEnv$plot_resolution, dpi = as.numeric(ssEnv$plot_resolution_ppi))


              # filename = paste0( chartFolder,"/",paste0( sample_group_comb, collapse ="_Vs_"),"_",marker,"_",figure, ".",ssEnv$plot_format,sep="")
              # grDevices::png(file= filename, width=2480, height = 2480, pointsize = 15, res = ssEnv$plot_resolution)
              # stats::heatmap(as.matrix(annotatedData[,3:ncol(annotatedData)]),
              #   col = grDevices::cm.colors(256),
              #   scale = "column",
              #   RowSideColors =as.vector(annotatedData$SAMPLE_GROUP),
              #   margins = c(25, 25),
              #   main = mainTitle
              # )
              # grDevices::dev.off()
            }

      }
  }
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Heatmap created." )

  }
}
