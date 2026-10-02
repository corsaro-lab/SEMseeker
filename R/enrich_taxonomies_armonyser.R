enrich_taxonomies_armonyser <- function(inference_details, result_folder, pvalue_column="PVALUE_ADJ_ALL_BH",
  significance = TRUE, disease_description,keywords,stop_keywords,alphas,disease_hpo, ...)
{
  if(length(disease_hpo)>0)
    disease_original <- gsub("[:]","_",disease_hpo)
  else
    disease_original <- ""

  ssEnv <- core_init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  keys <- unique(ssEnv$keys_for_pathway)
  inference_details <- as.data.frame(inference_details)
  pvalue_column <- core_name_cleaning(pvalue_column)
  key_enrichment_format <- ssEnv$key_enrichment_format
  aggregated_patwhay_result_total <- data.frame()

  for (a in alphas)
  {
    ssEnv$alpha <- a
    core_update_session_info(ssEnv)
    for (pt in seq_len(nrow(key_enrichment_format)))
    {
      if(key_enrichment_format[pt,"label"]=="phenolyzer")
        disease <- disease_original
      else
        disease <- ""

      column_of_id <- key_enrichment_format[pt,"column_of_id"]
      column_of_adj_pvalue <- key_enrichment_format[pt,"column_of_adj_pvalue"]
      column_of_description <- key_enrichment_format[pt,"column_of_description"]
      column_of_enrichment <- key_enrichment_format[pt,"column_of_enrichment"]

      for (id in seq_len(nrow(inference_details)))
      {
        # id <- 1
        # id <- 2
        inference_detail <- inference_details[id,]
        if(key_enrichment_format[pt,"type"]=="Pathway")
          path <- io_enrichment_folder(inference_detail, key_enrichment_format[pt,"label"])
        else
          path <- io_enrichment_folder(inference_detail, key_enrichment_format[pt,"label"], base = ssEnv$result_folderPhenotype)
        family_test <- inference_detail$family_test
        transformation_y <- as.character(inference_detail$transformation_y)
        independent_variable <- inference_detail$independent_variable
        covariates <- paste(inference_detail$covariates, collapse="_")

        file_prfx <- paste(independent_variable, transformation_y,family_test,covariates, pvalue_column,ssEnv$alpha, sep="_")

        missed_keys <- data.frame()
        for (i in seq_len(nrow(keys)))
        {
          file_name <- enrich_phenotype_analysis_name(inference_detail = inference_detail,key = keys[i,], prefix="",
            suffix=ifelse(disease=="","",paste("_",disease,sep="")),
            pvalue_column=pvalue_column, alpha = ssEnv$alpha, significance = significance)
          file_name <- io_file_path_build(path,file_name,"csv")
          # chech agin with disease without underscore
          if(!file.exists(file_name))
          {
            file_name <- enrich_phenotype_analysis_name(inference_detail = inference_detail,key = keys[i,], prefix="",
              suffix=ifelse(disease=="","",paste(disease,sep="")),
              pvalue_column=pvalue_column, alpha = ssEnv$alpha, significance = significance)
            file_name <- io_file_path_build(path,file_name,"csv")
          }

          # check with without_signal at the end for pathfindR
          if(!file.exists(file_name))
          {
            file_name <- enrich_phenotype_analysis_name(inference_detail = inference_detail,key = keys[i,], prefix="",
              suffix=ifelse(disease=="","",paste(disease,sep="")),
              pvalue_column=pvalue_column, alpha = ssEnv$alpha, significance = significance)
            file_name <- io_file_path_build(path,c(file_name,"without_signal"),"csv")
          }


          # #
          if (!file.exists(file_name))
          {
            # add a row if the file is missed it means non pathway was found !
            missed_keys <- plyr::rbind.fill(missed_keys, keys[i,])
            core_log_event("WARNING:", format(Sys.time(), "%a %b %d %X %Y") , " pathway_result ", file_name, " is missed !")
            next
          }

          if(key_enrichment_format[pt,"label"]=="phenolyzer")
          {
            # read the pathway_result
            pathway_result <- utils::read.csv2(file_name, dec=".")
          }
          else
            pathway_result <- utils::read.csv2(file_name)

          pathway_result <- enrich_analysy_add_category(data =  pathway_result, source =  key_enrichment_format[pt,"label"])
          cols_to_check <- c(column_of_id,column_of_enrichment,column_of_description, column_of_adj_pvalue)
          # check column names contain the required columns
          if(!all(cols_to_check %in% colnames(pathway_result)))
          {

            core_log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y") , " column_of_id ", column_of_id, " is missed in ", file_name)
            next
          }

          pathway_result$MARKER <- keys[i,"MARKER"]
          pathway_result$FIGURE <- keys[i,"FIGURE"]
          # pathway_result$key <- paste(pathway_result$MARKER,pathway_result$FIGURE,pathway_result$AREA,pathway_result$SUBAREA,sep="_")

          if(key_enrichment_format[pt,"label"]=="pathfindR")
          {
            # adjust pvalue
            pathway_result$PVALUE_ADJ_ALL_FDR <-stats::p.adjust(pathway_result[,"highest_p"], method = "fdr")

          }
          pathway_result$by_keyword <- FALSE
          if(length(keywords)>0)
            pathway_result[,"by_keyword"] <- vapply(pathway_result[,column_of_description], function(x) any(grepl(paste(keywords, collapse = "|"), x, ignore.case = TRUE)), logical(1))
          if(length(stop_keywords)>0)
            # remove item with the stop kewywords where by_keyword is TRUE
            pathway_result[pathway_result$by_keyword,"by_keyword"] <- vapply(pathway_result[pathway_result$by_keyword,column_of_description], function(x) !any(grepl(paste(stop_keywords, collapse = "|"), x, ignore.case = TRUE)), logical(1))

          # assign the same column name
          colnames(pathway_result)[colnames(pathway_result)==column_of_id] <- "ID"
          colnames(pathway_result)[colnames(pathway_result)==column_of_adj_pvalue] <- "FDR"
          colnames(pathway_result)[colnames(pathway_result)==column_of_description] <- "DESCRIPTION"
          colnames(pathway_result)[colnames(pathway_result)==column_of_enrichment] <- "ENRICHMENT"

          # take only IS,FDR,DESCRIPTION,ENRICHMENT
          pathway_result <- pathway_result[,c("ID","FDR","DESCRIPTION","ENRICHMENT","by_keyword","SS_CATEGORY")]
          pathway_result$MARKER <- keys[i,"MARKER"]
          pathway_result$alpha <- a
          pathway_result$key <- keys[i,"COMBINED"]
          pathway_result$independent_variable <- independent_variable
          pathway_result$covariates <- covariates
          pathway_result$transformation_y <- transformation_y
          pathway_result$family_test <- family_test
          pathway_result$areas_sql_condition <- inference_detail$areas_sql_condition
          pathway_result$samples_sql_condition <- inference_detail$samples_sql_condition
          pathway_result$association_results_sql_condition <- inference_detail$association_results_sql_condition
          pathway_result$combined_key <- paste(
            (pathway_result$alpha),
            (pathway_result$key),
            (pathway_result$family_test),
            (pathway_result$areas_sql_condition),
            (inference_detail$samples_sql_condition),
            sep="_")


          aggregated_patwhay_result_total <- plyr::rbind.fill(aggregated_patwhay_result_total, pathway_result)
        }
      }
    }
  }
  aggregated_patwhay_result_total$combined_key_number <- as.numeric(as.factor(aggregated_patwhay_result_total$combined_key))
  return(aggregated_patwhay_result_total)
}
