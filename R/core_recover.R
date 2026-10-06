core_recover <- function(result_folder, maxResources = 90,
  parallel_strategy  = "multisession",start_fresh = FALSE, ...)
{

  arguments <- list(...)
  areas_selection <- c()
  if(!is.null(arguments[["areas_selection"]]))
  {
    areas_selection <-arguments$areas_selection
    arguments[["areas_selection"]] <- NULL
  }

  ssEnv <- core_init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)
  core_log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), " SEMseeker will try to core_recover the analysys for project \n in ", ssEnv$result_folderData)

  study_summary <-   sem_study_summary_get()

  anno_create_position_pivots(study_summary,ssEnv$keys_markers_figures)
  sem_deltaX_get()
  sem_study_summary_total()
  anno_annotate_position_pivots()

  # Single point of sidecar materialisation.
  core_ensure_sidecars(ssEnv$result_folderData)
}
