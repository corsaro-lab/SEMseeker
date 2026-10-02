io_plot_file_name <- function(inference_detail, folder,file_extension="png", prefix = "", suffix="", key=c())
{
  ssEnv <- core_get_session_info()
  area <- as.character(key$AREA)
  subarea <- as.character(key$SUBAREA)
  marker <- as.character(key$MARKER)
  figure <- as.character(key$FIGURE)

  covariates <- inference_detail$covariates
  covariates <- if(length(covariates) !=  0 && !is.null(covariates)) unlist(t(strsplit( gsub(" ","",covariates),split  =  "+", fixed  =  TRUE)))
  covariates_dummy <- inference_detail$covariates_dummy
  covariates_dummy <- if(length(covariates_dummy) !=  0 && !is.null(covariates_dummy)) unlist(t(strsplit( gsub(" ","",covariates_dummy),split  =  "+", fixed  =  TRUE)))
  covariates_pca <- ifelse(is.null(inference_detail$covariates_pca),FALSE,inference_detail$covariates_pca)

  family_test <- inference_detail$family_test
  transformation_y <- inference_detail$transformation_y
  independent_variable <- gsub(" ","", inference_detail$independent_variable)

  file_suffix <- ""

  if(!is.null(covariates) & (length(covariates) !=0))
  {
    long_covariates <- length(covariates) > 2
    # split each covariates by _
    if (long_covariates)
    {
      covariates <- unlist(t(strsplit( gsub(" ","",covariates),split  =  "_", fixed  =  TRUE)))
      covariates <- unique(covariates)
    }
    file_suffix <- paste(covariates, collapse = "_")
  }

  if (!is.null(covariates_dummy) && (length(covariates_dummy)  !=  0))
  {
    long_covariates <- length(covariates_dummy) > 2
    # split each covariates_dummy by _
    if (long_covariates)
    {
      covariates_dummy <- unlist(t(strsplit( gsub(" ","",covariates_dummy),split  =  "_", fixed  =  TRUE)))
      covariates_dummy <- unique(covariates_dummy)
    }
    file_suffix <- c(file_suffix, paste(covariates_dummy, collapse = "_"))
  }

  if(length(covariates_dummy))
    file_suffix <- c(file_suffix, paste("dummy", sep = "_"))

  if(covariates_pca)
    file_suffix <- c(file_suffix, paste("pca", sep = "_"))

  file_suffix <- paste(file_suffix, suffix, sep = "_")

  file_result_prefix <- paste(as.character(family_test),as.character(independent_variable),as.character(transformation_y),sep="_")
  file_result_prefix <- paste(prefix, file_result_prefix, sep = "_")

  if(area != "")
    file_result_prefix <- paste(file_result_prefix, area, sep = "_")

  if(subarea != "")
    file_result_prefix <- paste(file_result_prefix, subarea, sep = "_")

  if(length(inference_detail$samples_sql_condition)>0)
    folder <- io_dir_check_and_create(folder,core_name_cleaning(inference_detail$samples_sql_condition))

  fileNameResults <- io_file_path_build(folder,c(as.character(marker),as.character(figure),file_result_prefix ,file_suffix) , file_extension)

  return(fileNameResults)
}


