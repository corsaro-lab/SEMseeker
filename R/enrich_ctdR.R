#' @importFrom doRNG %dorng%
#' @importFrom doFuture %dofuture%
enrich_ctdR <- function(study,
  statistic_parameter="",
  adjust_per_area = FALSE, adjust_globally = FALSE,adjustment_method = "BH", pvalue_column="PVALUE_ADJ_ALL_BH",
  inference_detail, significance = TRUE)
{

  # AI-311: what an enrichment can be about is declared once, not written
  # out again here. A pathway is a set of genes, so the input is one row
  # per gene (SCOPE = INSTANCE) of the GENE region class.
  .enrich_in <- enrich_input_invariant()

  #
  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  pvalue_column <- core_name_cleaning(pvalue_column)

  # start_fresh <- FALSE
  # ssEnv <- core_init_env(result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = start_fresh, ...)
  ssEnv <- core_get_session_info()
  keys <- unique(ssEnv$keys_for_pathway)

  #check if optional package is installed
  if(!requireNamespace("ctdR", quietly = TRUE))
  {
    core_log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y"),"ctdR package is not installed. Please install ctdR package to use this function.")
    return()
  }

  total_progress <- nrow(keys)*nrow(inference_detail)
  progress <- 0

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = seq_len(total_progress))
  else
    progress_bar <- ""


  for(id in seq_len(nrow(inference_detail)))
  {
    inference_detail <- inference_detail[id,]
    # foreach::foreach(i = 1:nrow(keys)) %dorng%
    for(i in seq_len(nrow(keys)))
    {
      seq <- 0
      progress <- progress + 1

      if(exists("existing_db"))
        existing_keys <- unique(pathway_report[pathway_report$source==db,"key"])

      if(ssEnv$showprogress)
        progress_bar(sprintf("Searching for pathway using ctdR for %s",keys[i,]$COMBINED))

      key <- paste(keys[i,]$AREA,keys[i,]$SUBAREA,keys[i,]$MARKER,keys[i,]$FIGURE, sep="_")

      suffix <- ""
      if(statistic_parameter=="")
        suffix <- "without_signal_"

      enrich_phenotype_analysis_name <- enrich_phenotype_analysis_name(inference_detail, keys[i,],prefix ="", suffix= suffix , pvalue_column, ssEnv$alpha, significance)
      path <- io_dir_check_and_create(ssEnv$result_folderEnrichment,c("ctdR",core_name_cleaning(inference_detail$areas_sql_condition),core_name_cleaning(inference_detail$samples_sql_condition), core_name_cleaning(inference_detail$association_results_sql_condition)))
      pathway_report_path <- io_file_path_build(path,enrich_phenotype_analysis_name,"csv")

      # if(file.exists(pathway_report_path))
      #   next
      if(file.exists(pathway_report_path))
      {
        pp <- utils::read.csv2(pathway_report_path,stringsAsFactors = FALSE)
        if(nrow(pp)==0)
          next
        enrich_result_save(pp, pathway_report_path, "ctdR")
        next
      }

      results_inference <- assoc_results_get(
        inference_detail =  inference_detail,
        # AI-257: enrichment happens for genes and nothing else — a pathway is a set
        # of genes. And it needs a p-value PER gene, so a collapsed artefact (one
        # number per sample) has nothing to list. Two coordinates, both invariant.
        area  = .enrich_in$area,
        scope = .enrich_in$scope,
        marker = keys[i,"MARKER"],
        pvalue_column=  pvalue_column,
        adjustment_method = adjustment_method,
        significance = TRUE)

      results_inference$AREA_OF_TEST <- results_inference$ENTREZID

      # remove rows wher AREA_OF_TEST is na
      results_inference <- results_inference[!is.na(results_inference$AREA_OF_TEST),]

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


      # gene_set <- na.omit(gene_set)
      if (statistic_parameter=="")
      {
        gene_set <- gene_set[,c("AREA_OF_TEST",pvalue_column)]
        gene_set <- aggregate(gene_set[,pvalue_column], by = list(gene_set$AREA_OF_TEST), mean)
        colnames(gene_set) <- c("AREA_OF_TEST",pvalue_column)
        gene_set$statistic_parameter <- NA
        gene_set <- gene_set[,c(pvalue_column,"statistic_parameter","AREA_OF_TEST")]
      }
      else
      {
        if(!(statistic_parameter %in% colnames(gene_set)))
          core_log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y"),  "this column", statistic_parameter, " doesn't exists.")
        #
        if(nrow(gene_set)<2)
          next
        gene_set <- gene_set[,c("AREA_OF_TEST",statistic_parameter,pvalue_column)]
        # aggregate the gene set using the max of pvalue_columns
        gene_set_pval <- aggregate(gene_set[,pvalue_column], by = list(gene_set$AREA_OF_TEST), mean)
        colnames(gene_set_pval) <- c("AREA_OF_TEST",pvalue_column)
        gene_set_stat <- aggregate(gene_set[,statistic_parameter], by = list(gene_set$AREA_OF_TEST), mean)
        colnames(gene_set_stat) <- c("AREA_OF_TEST",statistic_parameter)
        gene_set <- merge(gene_set_pval,gene_set_stat,by="AREA_OF_TEST")
        gene_set <- gene_set[,c(pvalue_column,statistic_parameter,"AREA_OF_TEST")]
        gene_set <- na.omit(gene_set)
      }

      # rename column to  pvalue, logfc and gene

      gene_set <- unique(gene_set)
      gene_set <- gene_set[order(gene_set[,pvalue_column],decreasing = FALSE),]
      key <- paste(keys[i,]$AREA,keys[i,]$SUBAREA,keys[i,]$MARKER,keys[i,]$FIGURE, sep="_")

      core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Number of genes in the gene set: ",nrow(gene_set), " key: ", keys[i,])

      if(nrow(gene_set)<5)
        next

      gene_set <- gene_set[, c("AREA_OF_TEST",pvalue_column)]
      colnames(gene_set) <- c("entrez_ids", "pvalue")
      result_pathway <- ctdR::enrichment_CTD(gene_set)
      result_pathway <- as.data.frame(result_pathway)

      # remove leadingEdge
      result_pathway$leadingEdge <- NULL

      core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Started ctdR analysis")


      if(nrow(result_pathway)==0)
        next

      # Assume total background number)
      result_pathway$key <- key
      result_pathway$source <- "ctdR"

      core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Done ctdR analysis")


      if(exists("result_pathway"))
      {
        enrich_result_save(result_pathway, pathway_report_path, "ctdR")      }
    }
  }
}
