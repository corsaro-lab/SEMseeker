enrich_pathfindR_circlize <- function(
    pathways_selection,
  statistic_parameter="", adjust_per_area = FALSE, adjust_globally = FALSE,adjustment_method = "BH",
  pvalue_column="PVALUE_ADJ_ALL_BH",
  inference_details, significance = TRUE,
  result_folder, maxResources = 90, parallel_strategy  = "multisession", ...)
{

  ssEnv <- core_init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)
  keys <- unique(ssEnv$keys_for_pathway)
  total_progress <- nrow(keys)*nrow(inference_details)
  progress <- 0

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = seq_len(total_progress))
  else
    progress_bar <- ""

  for(id in seq_len(nrow(inference_details)))
  {
    inference_detail <- inference_details[id,]
    for(i in seq_len(nrow(keys)))
    {
      progress <- progress + 1
      if(ssEnv$showprogress)
        progress_bar(sprintf("Creating circle for pathway for %s",keys[i,]$COMBINED))
      key <- paste(keys[i,]$AREA,keys[i,]$SUBAREA,keys[i,]$MARKER,keys[i,]$FIGURE, sep="_")
      suffix <- ""
      if(statistic_parameter=="")
        suffix = "without_signal_"

      enrich_phenotype_analysis_name <- enrich_phenotype_analysis_name(inference_detail, keys[i,],prefix ="", suffix= suffix , pvalue_column, ssEnv$alpha, significance)
      path <- io_dir_check_and_create(ssEnv$result_folderEnrichment,c("pathfindR",core_name_cleaning(inference_detail$areas_sql_condition),core_name_cleaning(inference_detail$samples_sql_condition), core_name_cleaning(inference_detail$association_results_sql_condition)))
      pathway_report_path <- io_file_path_build(path,enrich_phenotype_analysis_name,"csv")
      message("Pathway report path: ", pathway_report_path)
      if(file.exists(pathway_report_path))
      {
        pathways <- utils::read.csv2(pathway_report_path, header = TRUE)
        if(length(pathways_selection)>0)
          pathways <- subset(pathways, ID %in% pathways_selection)
        results <- data.frame()
        pathways <- as.data.frame(pathways)
        for ( n in seq_len(nrow(pathways)))
        {
          Up_regulated <- gsub(" ","",pathways[n,]$Up_regulated)
          if(!is.na(Up_regulated))
            Up_regulated <- strsplit(Up_regulated,",")[[1]]
          else
            Up_regulated <- c()
          Down_regulated <- gsub(" ","",pathways[n,]$Down_regulated)
          if(!is.na(Down_regulated))
            Down_regulated <- strsplit(Down_regulated,",")[[1]]
          else
            Down_regulated <- c()
          genes <- c(Down_regulated, Up_regulated)
          # remove spaces
          genes <- gsub(" ", "", genes)
          tmp <- expand.grid(pathways[n,"ID"], genes)
          results <- plyr::rbind.fill(tmp, results)
        }
        results <- unique(results)
        colnames(results) <- c("PATHWAY","GENE")
        # table(results$GENE)


        entrez_ids <- AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db,
          keys = as.vector(results$GENE),
          column = "ENTREZID",
          keytype = "ALIAS",
          multiVals = "first"
        )
        table(is.na(entrez_ids))

        official_symbols <- AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db,
          keys = entrez_ids,
          column = "SYMBOL",
          keytype = "ENTREZID",
          multiVals = "first"
        )
        table(is.na(official_symbols))
        results$GENE <- official_symbols

        results$PATHWAY_CHR <- "pathway"
        ssEnv_local <- core_get_session_info()
        probe_anno  <- anno_probe_annotation_build(ssEnv_local$tech)
        genes_impacted_cytoband <- unique(
          probe_anno[probe_anno$GENE_WHOLE %in% results$GENE,
                     c("CHR", "START", "END", "CHR_CYTOBAND", "GENE_WHOLE")])

        results <- merge(results,genes_impacted_cytoband, by.x="GENE", by.y="GENE_WHOLE")
        #reduce gene getting the start as min of start and end as max of end
        gene_start <- aggregate(START ~ GENE, data = results, min)
        gene_end <- aggregate(END ~ GENE, data = results, max)
        # remove START and END from results
        results <- results[,c("CHR","CHR_CYTOBAND","GENE","PATHWAY","PATHWAY_CHR")]
        results <- merge(results, gene_start, by.x="GENE", by.y="GENE")
        results <- merge(results, gene_end, by.x="GENE", by.y="GENE")
        results <- unique(results)
        results <- results[,c("CHR","START","END","GENE","PATHWAY","PATHWAY_CHR")]
        results$CHR <- paste0("chr",results$CHR)

        pathway_start <- aggregate(START ~ PATHWAY, data = results, min)
        colnames(pathway_start) <- c("PATHWAY","PATHWAY_START")
        # create progressive start
        pathway_start$PATHWAY_START <- seq(1, by = 2e8/nrow(pathway_start), length.out = nrow(pathway_start))
        pathway_start$PATHWAY_END <- pathway_start$PATHWAY_START + 2e8/nrow(pathway_start)

        results <- merge(results, pathway_start, by.x="PATHWAY", by.y="PATHWAY")
        results <- results[,c("CHR","START","END","GENE","PATHWAY","PATHWAY_CHR","PATHWAY_START","PATHWAY_END")]


        cytoband <- circlize::read.cytoband()$df
        cytoband <- rbind(cytoband,
          data.frame(V1 = "pathway", V2 = 1,  V3 = 2e8, V4 = "", V5 = "")
        )
        # circos.initializeWithIdeogram(cytoband)
        # circos.genomicLink(results[, 1:3], results[, 6:8])
        # circos.clear()
        ####

        # circos.par(gap.after = c(rep(1, 23), 5, 5))
        # circos.genomicInitialize(cytoband, plotType = NULL)
        #
        # # the labels track
        # label_df = rbind(results[, 1:4], stats::setNames(results[, c(6:8, 5)], colnames(results)[1:4]))
        # label_df = unique(label_df)
        # label_df = label_df[,c("CHR","START","END","GENE")]
        # circos.genomicLabels(label_df, labels.column = 4, side = "outside", cex = 0.6)
        #
        # # the chromosome names track
        # circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
        #   circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        #     niceFacing = TRUE, adj = c(0.5, 0), cex = 0.8)
        # }, track.height = graphics::strheight("fj", cex = 0.8)*1.2, bg.border = NA, cell.padding = c(0, 0, 0, 0))
        #
        # # ideogram track
        # circos.genomicIdeogram(cytoband)
        #
        # # genomic links
        # my_col <- grDevices::rainbow(nrow(results))
        # circos.genomicLink(results[, 1:3], results[, 6:8], col = my_col)


        # Variables to control the font size, link length, chr font size, and link start height
        label_font_size <- 0.3
        chr_font_size <- 0.3

        label_connection_height <- 0.03

        link_start_height <- 2
        link_length <- 2

        enrich_phenotype_analysis_name <- enrich_phenotype_analysis_name(inference_detail, keys[i,],prefix ="", suffix= suffix , pvalue_column, ssEnv$alpha, significance)
        path <- io_dir_check_and_create(ssEnv$result_folderChart,"pathfindR")
        pathway_report_path <- io_file_path_build(path,enrich_phenotype_analysis_name,"png")
        # Save the plot as a PNG file
        if(ssEnv$plot_format == "png")
          grDevices::png(file =  pathway_report_path, width = 2480,height = 2480, pointsize  =  15, res = as.numeric(ssEnv$plot_resolution_ppi))
        if(ssEnv$plot_format == "eps")
          grDevices::postscript(file =  pathway_report_path, width = 2480,height = 2480, pointsize  =  15, res = as.numeric(ssEnv$plot_resolution_ppi))
        # grDevices::png(pathway_report_path, width = 2048, height = 2048, res = as.numeric(ssEnv$plot_resolution_ppi))

        # Set circos parameters to adjust the gap and overall appearance
        circlize::circos.par(gap.after = c(rep(1, 23), 5, 5), track.height = link_start_height)

        # Initialize the circos plot without plotting anything yet
        circlize::circos.genomicInitialize(cytoband, plotType = NULL)

        # Prepare label data frame
        label_df <- rbind(results[, 1:4], stats::setNames(results[, c(6:8, 5)], colnames(results)[1:4]))
        label_df <- unique(label_df)
        label_df <- label_df[, c("CHR", "START", "END", "GENE")]

        # Add genomic labels with adjusted size and position
        circlize::circos.genomicLabels(label_df, labels.column = 4, side = "outside", cex = label_font_size, connection_height = label_connection_height)

        # Add chromosome names with adjusted size and positioning
        circlize::circos.track(track.index = circlize::get.current.track.index(), panel.fun = function(x, y) {
          circlize::circos.text(circlize::CELL_META$xcenter, circlize::CELL_META$ylim[1], circlize::CELL_META$sector.index,
            niceFacing = TRUE, adj = c(0.5, 0), cex = chr_font_size)
        }, track.height = graphics::strheight("fj", cex = chr_font_size) * 2, bg.border = "grey", bg.col= "white", cell.padding = c(0, 0, 0, 0))

        # Add the ideogram track
        circlize::circos.genomicIdeogram(cytoband)

        # Add genomic links with unique colors
        my_col <- grDevices::rainbow(nrow(results))
        circlize::circos.genomicLink(results[, 1:3], results[, 6:8], col = my_col)

        # Clear the circos plot after drawing
        circlize::circos.clear()

        # Close the PNG device
        grDevices::dev.off()
      }
    }
  }
}

