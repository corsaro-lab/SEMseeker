enrich_WebGestalt <- function(study,
  types=c("BP","MF"),  enrich_methods = c("ORA"),
  adjust_per_area = FALSE, adjust_globally = FALSE,adjustment_method = "BH", pvalue_column="PVALUE_ADJ_ALL_BH",
  inference_detail,significance)
{

  # AI-311: what an enrichment can be about is declared once, not written
  # out again here. A pathway is a set of genes, so the input is one row
  # per gene (SCOPE = INSTANCE) of the GENE region class.
  .enrich_in <- enrich_input_invariant()

  #
  # start_fresh <- FALSE
  # ssEnv <- core_init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = start_fresh, ...)
  ssEnv <- core_get_session_info()
  pvalue_column <- core_name_cleaning(pvalue_column)
  keys <- unique(ssEnv$keys_for_pathway)
  path <- io_dir_check_and_create(ssEnv$result_folderEnrichment,c("WebGestalt",core_name_cleaning(inference_detail$areas_sql_condition),core_name_cleaning(inference_detail$samples_sql_condition), core_name_cleaning(inference_detail$association_results_sql_condition)))
  tmp <- tempdir()
  tempFolder <- io_dir_check_and_create(tmp,c("/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]")))


  #check if optional package is installed
  if(!requireNamespace("WebGestaltR", quietly = TRUE))
  {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " WebGestaltR package is not installed. Please install WebGestaltR package to use this function")
    return()
  }

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = seq_len(nrow(keys)*length(types)*length(enrich_methods)))
  else
    progress_bar <- ""

  for (i in seq_len(nrow(keys)))
  {
    for( t in seq_along(types))
    {
      type <- types[t]
      for ( em in seq_along(enrich_methods))
      {
        projectName <- enrich_phenotype_analysis_name( inference_detail = inference_detail,key = keys[i,], prefix="",suffix=""  ,
          pvalue_column=pvalue_column, as.numeric(ssEnv$alpha), significance)
        filenameResult <- io_file_path_build(path,projectName,"csv")
        # if(file.exists(filenameResult))
        #   next

        if(file.exists(filenameResult))
        {
          pp <- utils::read.csv2(filenameResult,stringsAsFactors = FALSE)
          if(nrow(pp)==0)
            next
          enrich_result_save(pp, filenameResult, "WebGestalt")
          next
        }

        #

        enrich_method <- enrich_methods[em]
        if(ssEnv$showprogress)
          progress_bar(sprintf("Searching for disease using WebGestalt: %s with %s and %s",keys[i,]$COMBINED,enrich_method,type))
        key <- paste(keys[i,]$FIGURE,keys[i,]$MARKER,keys[i,]$AREA,keys[i,]$SUBAREA, sep="_")

        results_inference <- assoc_results_get(
          inference_detail =  inference_detail,
          # AI-257: enrichment happens for genes and nothing else — a pathway is a set
          # of genes. And it needs a p-value PER gene, so a collapsed artefact (one
          # number per sample) has nothing to list. Two coordinates, both invariant.
          area  = .enrich_in$area,
          scope = .enrich_in$scope,
          marker = keys[i,"MARKER"],
          pvalue_column=  pvalue_column,
          adjustment_method= adjustment_method,
          significance = TRUE)

        #

        if (nrow(results_inference)==0)
          next

        if(keys[i,]$SUBAREA=="ALL_SUBAREAS")
          gene_set <- results_inference[results_inference$SUBAREA!="WHOLE",]
        else
          gene_set <- results_inference[results_inference$SUBAREA==keys[i,]$SUBAREA,]

        if(keys[i,]$FIGURE=="HYPER_HYPO")
          gene_set <- results_inference[results_inference$FIGURE=="HYPER" | results_inference$FIGURE=="HYPO",]
        else
          gene_set <- results_inference[results_inference$FIGURE==keys[i,]$FIGURE,]

        if(nrow(gene_set)==0)
        {
          core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " No genes found for the key: ", keys[i,])
          next
        }
        # remove duplicates
        gene_set <- aggregate(gene_set[,pvalue_column], by = list(gene_set$AREA_OF_TEST), max)

        colnames(gene_set) <- c("AREA_OF_TEST",pvalue_column)
        gene_set <- gene_set[!duplicated(gene_set$AREA_OF_TEST),]
        gene_set <- gene_set[order(gene_set[,pvalue_column]),]
        gene_set <- na.omit(gene_set)

        core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Number of genes in the gene set: ",nrow(gene_set), " key: ", keys[i,])
        if (nrow(gene_set)==0)
          next
        projectName <- enrich_phenotype_analysis_name( inference_detail = inference_detail,key = keys[i,], prefix="",suffix= paste(type,enrich_method, sep="_") , pvalue_column=pvalue_column, as.numeric(ssEnv$alpha), significance)

        # entrez <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = as.vector(gene_set$AREA_OF_TEST),column = "ENTREZID", keytype = "SYMBOL")
        # entrez <- unique(entrez)
        # entrez <- na.omit(entrez)

        # core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Number of genes in the gene set (ENTREZ): ",nrow(gene_set), " key: ", keys[i,])
        if (length(gene_set$AREA_OF_TEST)==0)
          next

        geneFile <- file.path(system.file(package="WebGestaltR"),"extdata/interestingGenes.txt")
        write.table(unique(gene_set$AREA_OF_TEST),geneFile,row.names=FALSE,col.names = FALSE,quote =FALSE)

        enrichDataBase <- switch(
          type,
          "BP"="geneontology_Biological_Process",
          "MF"="geneontology_Molecular_Function"
        )

        if(nrow(gene_set)<2)
          next


        enrichResult <- tryCatch({
          WebGestaltR::WebGestaltR(
            enrichMethod=enrich_method,
            organism="hsapiens",
            enrichDatabase= enrichDataBase,
            interestGeneFile=geneFile,
            interestGeneType="genesymbol",
            referenceSet= "genome",
            isOutput=FALSE,
            outputDirectory= tempFolder,
            projectName= projectName,
            sigMethod ="top",
            fdrMethod = "BH",
            fdrThr = 1,
            topThr = 200,
            maxNum = 2000,
            minNum = 1
          )
        },
          catch = function(e) {
            core_log_event(paste("Error in WebGestaltR: ",e))
            NULL
          },
          finally = {
            # core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " WebGestaltR done for key: ", keys[i,])
          }
        )

        if(is.null(enrichResult))
          next

        if(nrow(enrichResult)==0)
          next
        enrichResult$ALPHA <- as.numeric(ssEnv$alpha)
        enrichResult$pvalue_column <- pvalue_column
        enrichResult$type <- type
        enrichResult$enrich_method <- enrich_method
        enrichResult$key <- key
        enrichResult$MARKER <- keys[i,]$MARKER
        enrichResult$FIGURE <- keys[i,]$FIGURE
        enrichResult$AREA <- keys[i,]$AREA
        enrichResult$SUBAREA <- keys[i,]$SUBAREA

        # format FDR to use scientific notation

        enrichResultToPlot <- enrichResult[1:15,]
        # # browser
        # if(min(enrichResultToPlot$FDR) < as.numeric((ssEnv$alpha)))
        # {
        plotFileName <- io_file_path_build(path,projectName,ssEnv$plot_format)
        # grDevices::png(file= plotFileName, width=2048,height=2048, bg = "transparent")
        enrichResultToPlot <- as.data.frame(enrichResultToPlot)
        enrichResultToPlot <- enrichResultToPlot[order(enrichResultToPlot$expect),]


        # i think plot os too much
        # # Create a new column for bar colors based on the FDR condition
        # enrichResultToPlot <- enrichResultToPlot %>%
        #   plyr::mutate(bar_colors = ifelse(FDR < as.numeric(ssEnv$alpha), ssEnv$color_palette[1], ssEnv$color_palette[2]))
        #
        #
        # # Generate the ggplot
        # gg <- ggplot2::ggplot(enrichResultToPlot, ggplot2::aes(x = expect, y = reorder(description, expect), fill = bar_colors)) +
        #   ggplot2::geom_bar(stat = "identity", color = "black") +
        #   ggplot2::scale_fill_identity() +
        #   ggplot2::labs(
        #     title = paste0(keys[i, "AREA"], " ", keys[i, "SUBAREA"], " ", keys[i, "MARKER"], " ", keys[i, "FIGURE"], " ", pvalue_column),
        #     x = "Expected",
        #     y = ""
        #   ) +
        #   ggplot2::theme_minimal() +
        #   ggplot2::theme(
        #     axis.text.y = ggplot2::element_text(size = 12),
        #     axis.title.x = ggplot2::element_text(size = 12),
        #     plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
        #     legend.position = "none"
        #   )
        #
        #
        # # Save the plot
        # ggplot2::ggsave(plotFileName, gg, width = 4096, height = 2048, units = "px", dpi = as.numeric(ssEnv$plot_resolution_ppi))

        file.remove(geneFile)
        unlink(geneFile)

        if(exists("enrichResult"))
          if(exists("enrichResultFinal"))
            enrichResultFinal <- rbind(enrichResultFinal,enrichResult)
        else
          enrichResultFinal <- enrichResult

      }
    }

    if(exists("enrichResultFinal"))
    {
      projectName <- enrich_phenotype_analysis_name( inference_detail = inference_detail,key = keys[i,], prefix="",suffix=""  , pvalue_column=pvalue_column, as.numeric(ssEnv$alpha), significance)
      filenameResult <- io_file_path_build(path,projectName,"csv")
      enrich_result_save(enrichResultFinal, filenameResult, "WebGestalt")
      rm(enrichResultFinal)
    }
  }
}
