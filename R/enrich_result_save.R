enrich_result_save <- function(result_pathway, pathway_report_path, pathway_package, study)
{
  ssEnv <- core_get_session_info()
  if(nrow(result_pathway)!=0)
  {
    result_pathway <- enrich_analysy_add_category(pathway_package,result_pathway)
    pvalue_column_adj <- ssEnv$key_enrichment_format[ssEnv$key_enrichment_format$label==pathway_package,"column_of_pvalue"]
    description_column <- ssEnv$key_enrichment_format[ssEnv$key_enrichment_format$label==pathway_package,"column_of_description"]
    # result_pathway <- result_pathway[, !(grepl("PVALUE_ADJ_", colnames(result_pathway))]
    col_p <- core_name_cleaning(paste0("PVALUE_ADJ_", ssEnv$multiple_test_adj))
    tryCatch({
      if(ssEnv$multiple_test_adj=="q")
        result_pathway[,col_p] <- qvalue::qvalue(result_pathway[,pvalue_column_adj], fdr.level = ssEnv$alpha, pi0.method="bootstrap", na.rm=TRUE)$qvalues
      else
        result_pathway[,col_p] <- stats::p.adjust(result_pathway[,pvalue_column_adj],method  =  ssEnv$multiple_test_adj)
      # sort by col_p
      result_pathway <- result_pathway[order(result_pathway[,col_p]),]
    }, error = function(e) {

    })
    # study arrives as an argument. It used to be read as a bare name, which
    # only works under dynamic scoping: every caller happens to have a
    # parameter called study, but R resolves names lexically, so this line
    # failed with "object 'study' not found" on every single call.
    result_pathway$PHENOTYPE <- grepl(study,result_pathway[,description_column], ignore.case = TRUE)
    # REMOVE COLUMNS with NAMES X, X.1 and X.2
    result_pathway <- result_pathway[,!grepl("^X$|^X\\.[0-9]+$", colnames(result_pathway))]
    utils::write.csv2(result_pathway, pathway_report_path, row.names = FALSE)
    rm(result_pathway)
  }
}
