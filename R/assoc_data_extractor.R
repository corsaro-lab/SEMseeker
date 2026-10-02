assoc_data_extractor <- function(inference_details,destination_folder="", result_folder="", ...)
{

  arguments <- list(...)
  areas_selection <- c()
  if(!is.null(arguments[["areas_selection"]]))
  {
    areas_selection <-arguments$areas_selection
    arguments[["areas_selection"]] <- NULL
  }

  if(result_folder!="")
    ssEnv <- core_init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  else
    ssEnv <- core_get_session_info()

  localKeys <- ssEnv$keys_markers_figures
  markers <- unique(localKeys$MARKER)
  final_results <- data.frame()
  for(z in seq_len(nrow(inference_details)))
  {
    start_time <- Sys.time()
    inference_detail <- inference_details[z,]
    for (a in seq_along(markers) )
    {
      if (exists("results_inference"))
        rm(results_inference)

      keys <- localKeys[localKeys$MARKER==markers[a],]
      keys <- unique(keys)
      fileNameResults <- io_inference_file_name(inference_detail, markers[a], ssEnv$result_folderInference)
      if (file.exists(fileNameResults))
      {
        results_inference <- utils::read.csv2(fileNameResults, header  =  TRUE)
        results_inference <- subset(results_inference, FIGURE %in% keys$FIGURE)
        results_inference <- assoc_filter_sql(inference_detail$areas_sql_condition, results_inference)

        # remove any column with name containg SAMPLES_SQL_CONDITION
        results_inference <- results_inference[,!grepl("SAMPLES_SQL_CONDITION", colnames(results_inference))]

        # AI-257: this function used to write the filtered frame back over
        # `fileNameResults` - the very file it had just read. Extracting was
        # therefore destructive: the three filters above (the figures of the
        # marker, the areas_sql_condition of the request, the dropped
        # SAMPLES_SQL_CONDITION columns) were applied to the canonical inference
        # result and the rows that did not pass were gone from disk. Two
        # consequences, both silent: a narrow `areas_sql_condition` permanently
        # reduced the result of the run, and the next iteration of this same loop
        # read a file that a previous iteration had already truncated.
        #
        # A consumer reads. `destination_folder` below is where an extraction is
        # allowed to write, and it writes a copy under a name of its own.

        if(inference_detail$samples_sql_condition!="")
          results_inference$SAMPLES_SQL_CONDITION <- inference_detail$samples_sql_condition

        final_results <- plyr::rbind.fill(final_results, results_inference)
        # core_log_event("DEBUG: ",format(Sys.time(), "%a %b %d %X %Y")," sql executed:", areas_sql_condition)
        filePathResults <- file.path(destination_folder, paste0(inference_detail$family_test, "_", markers[a], ".csv"))
        if(destination_folder != "")
          utils::write.csv2(results_inference, filePathResults, row.names = FALSE)
      }
    }
  }
  return(final_results)
}


