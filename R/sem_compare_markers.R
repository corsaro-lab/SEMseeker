sem_compare_markers <- function(inference_detail, result_folder, studies, adjust_per_area = FALSE, adjust_globally = FALSE, pvalue_column = "PVALUE_ADJ_ALL_BH", adjustment_method = "BH", ...){

  ssEnv <- core_init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  color_palette <- ssEnv$color_palette
  localKeys <- ssEnv$keys_areas_subareas_markers_figures

  aggregated_markers_results <- data.frame()
  for (a in seq_len(nrow(localKeys)))
  {
    MARKER <- localKeys[a, "MARKER"]
    AREA <- localKeys[a, "AREA"]
    # for each study in studies
    for (s in seq_len(nrow(studies)))
    {
      # get the inference details for the study
      # AI-257: SCOPE declared. Without it this read mixed the collapsed row of a
      # region class with its per-instance rows, which are different quantities.
      temp_res <- assoc_results_get(inference_detail = inference_detail, marker = MARKER,
         area = AREA, scope = "INSTANCE",
         pvalue_column= pvalue_column,
        adjustment_method = adjustment_method,significance = TRUE)
      if(nrow(temp_res) != 0)
        temp_res$STUDY <- studies[s,"STUDY"]
      aggregated_markers_results <- plyr::rbind.fill(aggregated_markers_results, temp_res)
    }
  }
}
