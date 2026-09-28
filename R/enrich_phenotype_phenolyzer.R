enrich_phenotype_phenolyzer <- function(study,
  disease,phenolyzer_folder_bin,minimum_score = 0.5, statistic_parameter = "",
  adjust_per_area = FALSE, adjust_globally = FALSE,adjustment_method = "BH", pvalue_column="PVALUE_ADJ_ALL_BH",
  inference_detail, significance = TRUE)
{

  # AI-311: what an enrichment can be about is declared once, not written
  # out again here. A pathway is a set of genes, so the input is one row
  # per gene (SCOPE = INSTANCE) of the GENE region class.
  .enrich_in <- enrich_input_invariant()

  # start_fresh <- FALSE
  # ssEnv <- core_init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = start_fresh, ...)
  ssEnv <- core_get_session_info()
  pvalue_column <- core_name_cleaning(pvalue_column)

  # check existence of phenolyzer
  if (!file.exists(phenolyzer_folder_bin))
  {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Phenolyzer not found in ",phenolyzer_folder_bin)
    return()
  }

  tmp <- tempdir()

  #
  diseases <- enrich_disease_get(disease =  disease,vocabulary =  "HPO")
  disease_label <- gsub(":","_",disease)

  if (is.null(nrow(diseases)))
  {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " No disease found for ",disease_label)
    return()
  }

  keys <- unique(ssEnv$keys_for_pathway)

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = seq_len(nrow(keys)))
  else
    progress_bar <- ""

  diseases_summary <- data.frame()
  seq <-0
  for(i in seq_len(nrow(keys)))
  {
    random_string <- stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]")
    tempFolder <- io_dir_check_and_create(tmp,c("/semseeker/",random_string))

    if(ssEnv$showprogress)
      progress_bar(sprintf("Searching for disease using phenolyzer: %s",keys[i,]$COMBINED))
    key <- paste(keys[i,]$FIGURE,keys[i,]$MARKER,keys[i,]$AREA,keys[i,]$SUBAREA, sep="_")

    base_path <- io_dir_check_and_create(ssEnv$result_folderPhenotype,c("phenolyzer",core_name_cleaning(inference_detail$areas_sql_condition)))
    path <- io_dir_check_and_create(baseFolder = base_path, subFolders = "scores")
    annotated_gene_file <- file.path(tempFolder,"ex8.annotated_gene_scores")
    enrich_phenotype_analysis_name <- enrich_phenotype_analysis_name( inference_detail = inference_detail,key = keys[i,], prefix="",suffix=paste("_",disease_label,"_gene_scores",sep=""), pvalue_column=pvalue_column, ssEnv$alpha, significance)
    phenotype_report_path <- io_file_path_build(path,enrich_phenotype_analysis_name,"csv")
    #
    if(file.exists(phenotype_report_path))
      next

    results_inference <- assoc_results_get(
      inference_detail =  inference_detail,
      marker = keys[i,"MARKER"],
      # AI-257: neither coordinate was declared here, so this read took the
      # GENE default and every scope — the collapsed rows included. A pathway
      # needs a p-value per gene; a per-sample burden has no genes to list.
      area  = .enrich_in$area,
      scope = .enrich_in$scope,
      pvalue_column=  pvalue_column,
      adjustment_method= adjustment_method,
      significance = TRUE)

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

    colnames(gene_set) <- c("AREA_OF_TEST",pvalue_column)
    gene_set <- gene_set[order(gene_set[,pvalue_column]),]

    if (nrow(results_inference)==0)
      next

    enrich_phenotype_analysis_name <- enrich_phenotype_analysis_name( inference_detail = inference_detail,key = keys[i,], prefix="",suffix=disease_label, pvalue_column=pvalue_column, ssEnv$alpha, significance)
    phenotype_report_path <- io_file_path_build(base_path,enrich_phenotype_analysis_name,"csv")

    if(statistic_parameter!="")
      gene_set <- results_inference[,c("AREA_OF_TEST",statistic_parameter,pvalue_column)]
    else
      gene_set <- results_inference[,c("AREA_OF_TEST",pvalue_column)]
    gene_set <- na.omit(gene_set)

    if(nrow(gene_set)<2)
      next

    gene_set <- aggregate(gene_set[,pvalue_column], by = list(gene_set$AREA_OF_TEST), mean)
    colnames(gene_set) <- c("AREA_OF_TEST",pvalue_column)
    # reduce gene set to unique genes and taking the maximum pvalue
    gene_set <- gene_set[order(gene_set[,pvalue_column],decreasing = FALSE),]
    gene_set <- gene_set[!duplicated(gene_set$AREA_OF_TEST),]

    core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Number of genes in the gene set: ",nrow(gene_set), " key: ", keys[i,])

    file_term <- file.path(tempFolder, paste0("term_",random_string,".txt"))
    file_genes <- file.path(tempFolder, paste0("genes_",random_string,".txt"))
    write.table(unique(gene_set[,"AREA_OF_TEST"]), file_genes, quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(diseases, file_term, quote = FALSE, row.names = FALSE, col.names = FALSE)

    #
    nCore <-  ssEnv$parallel$nCore
    # perl disease_annotation.pl alzheimer -p -ph -logistic -out ex8 -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25
    pcommand <-
      paste0(
        "cd ",
        tempFolder ,
        " && perl ", phenolyzer_folder_bin ,"disease_annotation.pl " ,
        " -p " ,file_term,
        " -f --gene " ,file_genes,
        " -wordcloud ",
        " -prediction -phenotype -logistic -out ex8 ",
        " -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25 ",
        " -nproc ", nCore ,
        sep = ""
      )
    message(pcommand)
    try(
      {
        # sudo cpan App::cpanminus IO::String Bio::OntologyIO Graph::Directed JSON Parallel::ForkManager
        # biotools:phenolyzer
        # cpan Module::Name
        # cd /usr/local/lib/
        # git clone https://github.com/WGLab/phenolyzer.git
        # sudo cpan Bio::OntologyIO
        # sudo cpan Graph::Directed

        # perl /usr/local/lib/phenolyzer/disease_annotation.pl -p ./hypertension.txt  -f --gene ./genes_list.txt  -prediction -phenotype -logistic -out hypertension -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25 -nproc 15

        null_device <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
        message(pcommand)
        sink(null_device)
        # pcommand uses shell features ('cd ... && perl ...'); system2() does not
        # honour these directly, so route through the shell explicitly. Bioconductor
        # discourages bare system() in package code.
        shell_bin <- if (.Platform$OS.type == "windows") "cmd" else "bash"
        shell_arg <- if (.Platform$OS.type == "windows") "/c" else "-c"
        system2(shell_bin, c(shell_arg, pcommand))
        sink()
      }
    )

    annotated_gene_file <- file.path(tempFolder,"ex8.annotated_gene_list")
    if(!file.exists(annotated_gene_file))
    {
      core_log_event("INFO:" , format(Sys.time(), "%a %b %d %X %Y"), "No annotated gene file found !")
      next
    }

    annotated_genes <-  utils::read.csv2(annotated_gene_file, sep = "\t")
    annotated_genes_temp <- annotated_genes[annotated_genes$Status=="SeedGene",]

    if(nrow(annotated_genes_temp)==0)
      annotated_genes_temp <- annotated_genes[annotated_genes$Score>=minimum_score,]

    if(nrow(annotated_genes_temp)==0)
      next
    seq <- seq + 1
    if(nrow(annotated_genes)==0)
      next()

    annotated_genes$key <- key
    annotated_genes$seq <- seq
    annotated_genes$gene_count <- nrow(gene_set)
    annotated_genes$order <- seq_len(nrow(annotated_genes))
    annotated_genes$diseases <- paste(diseases$DISEASE, collapse = ";")

    enrich_phenotype_analysis_name <- enrich_phenotype_analysis_name( inference_detail = inference_detail,key = keys[i,], prefix="",suffix=paste("_", disease_label,"_report",sep=""), pvalue_column=pvalue_column, ssEnv$alpha, significance)
    path <- io_dir_check_and_create(baseFolder = base_path, subFolders = "summary")
    phenotype_report_path <- io_file_path_build(path,enrich_phenotype_analysis_name,"csv")
    utils::write.csv2(annotated_genes, phenotype_report_path)


    path <- io_dir_check_and_create(baseFolder = base_path, subFolders = "scores")
    annotated_gene_file <- file.path(tempFolder,"ex8.annotated_gene_scores")
    enrich_phenotype_analysis_name <- enrich_phenotype_analysis_name( inference_detail = inference_detail,key = keys[i,], prefix="",suffix=paste("_", disease_label,"_gene_scores",sep=""), pvalue_column=pvalue_column, ssEnv$alpha, significance)
    phenotype_report_path <- io_file_path_build(path,enrich_phenotype_analysis_name,"csv")
    file.copy(annotated_gene_file, phenotype_report_path, overwrite = TRUE)


    path <- io_dir_check_and_create(baseFolder = base_path, subFolders = "worldcloud")
    # copy all files ending with _wordcloud.png to the worldcloud folder
    wordcloud_files <- list.files(tempFolder, pattern = "_wordcloud.png", full.names = TRUE)
    for (w in seq_along(wordcloud_files))
    {
      wordcloud_file <- wordcloud_files[w]
      enrich_phenotype_analysis_name <- enrich_phenotype_analysis_name( inference_detail = inference_detail,key = keys[i,], prefix="",suffix=paste(disease_label,"_wordcloud",sep=""), pvalue_column=pvalue_column, ssEnv$alpha, significance)
      phenotype_report_path <- io_file_path_build(path,enrich_phenotype_analysis_name,"png")
      file.copy(wordcloud_file, phenotype_report_path, overwrite = TRUE)
    }

    # path <- io_dir_check_and_create(baseFolder = base_path, subFolders = "disease")
    # copy all files ending with _diseases.png to the worldcloud folder
    diseases_files <- list.files(tempFolder, pattern = "_diseases", full.names = TRUE)

    for (w in seq_along(diseases_files))
    {

      diseases_file <- diseases_files[w]
      disease_temp <- utils::read.csv2(diseases_file, sep = "\t", header = FALSE)
      colnames(disease_temp) <- c("Description","Score")
      disease_temp$key <- key
      disease_temp$HPO <- disease
      diseases_summary <- plyr::rbind.fill(diseases_summary, disease_temp)

      projectName <- enrich_phenotype_analysis_name( inference_detail = inference_detail,key = keys[i,], prefix="",suffix= disease_label  , pvalue_column=pvalue_column, as.numeric(ssEnv$alpha), significance)
      filenameResult <- io_file_path_build(base_path,projectName,"csv")
      utils::write.csv2(disease_temp, filenameResult)
    }

    # enrich_phenotype_analysis_name <- enrich_phenotype_analysis_name( inference_detail = inference_detail,key = keys[i,], prefix="",suffix=paste(disease,"_gene_cloud",sep=""), pvalue_column=pvalue_column, ssEnv$alpha)
    # cloud_gene_file <- file.path(tempFolder,"out.annotated_gene_scores")
    # file.copy(cloud_gene_file, phenotype_report_path, overwrite = TRUE)

    unlink(tempFolder, recursive = TRUE)
  }

  # family_test <- inference_detail$family_test
  # transformation_y <- inference_detail$transformation_y
  # significance_label <- ifelse(significance, "significant", "non_significant")
  # analysis_name <- paste("summary",pvalue_column, significance_label, ssEnv$alpha, as.character(transformation_y), as.character(family_test), sep="_")
  # phenotype_report_path <- io_file_path_build(base_path,analysis_name,"csv")
  # if(file.exists(phenotype_report_path))
  # {
  #   disease_temp <- utils::read.csv2(phenotype_report_path)
  #   diseases_summary <- plyr::rbind.fill(diseases_summary, disease_temp)
  # }
  # diseases_summary <- unique(diseases_summary)
  # write.csv2(diseases_summary, phenotype_report_path)
}
