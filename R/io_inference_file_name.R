io_inference_file_name <- function(inference_detail, marker, folder,file_extension="csv", prefix = "", suffix="",
                                  skip_dir_create = FALSE)
{
  ssEnv <- core_get_session_info()

  covariates <- util_split_and_clean(inference_detail$covariates, "\\+")
  covariates_dummy <- util_split_and_clean(inference_detail$covariates_dummy, "\\+")
  covariates_pca <- util_boolean_check(inference_detail$covariates_pca)
  family_test <- util_split_and_clean(inference_detail$family_test)
  transformation_y <- util_split_and_clean(inference_detail$transformation_y)
  inference_detail$independent_variable <- as.character(inference_detail$independent_variable)
  inference_detail$independent_variable <- as.character(gsub("_SCALED","",inference_detail$independent_variable))
  if(inference_detail$transformation_x=="scale")
    inference_detail$independent_variable <- paste0(inference_detail$independent_variable, "_SCALED")
  independent_variable <- util_split_and_clean(inference_detail$independent_variable)

  file_suffix <- ""
  file_result_prefix <- paste(as.character(independent_variable),as.character(transformation_y),as.character(family_test),sep="_")

  if(length(covariates)>0)
    file_suffix <- paste(covariates, collapse = "_")

  if(length(covariates_dummy)>0)
    file_suffix <- c(file_suffix, paste(covariates_dummy, collapse = "_"))

  if(length(covariates_dummy)>0)
    file_suffix <- c(file_suffix, paste("dummy", sep = "_"))

  if(covariates_pca)
    file_suffix <- c(file_suffix, paste("pca", sep = "_"))

  file_result_prefix <- paste(prefix, file_result_prefix, sep = "_")
  file_suffix <- paste(file_suffix, suffix, sep = "_")

  # skip_dir_create is for callers that want the NAME and must not touch the
  # filesystem: checking whether another study holds the results of a request
  # would otherwise create folders inside that study's tree. Same idiom as
  # io_bed_file_name().
  if(length(inference_detail$samples_sql_condition)>0)
  {
    if (skip_dir_create)
      folder <- file.path(folder, core_name_cleaning(inference_detail$samples_sql_condition))
    else
      folder <- io_dir_check_and_create(folder,core_name_cleaning(inference_detail$samples_sql_condition))
  }


  fileNameResults <- io_file_path_build(folder,c(as.character(marker),file_result_prefix ,file_suffix) , file_extension)

  return(fileNameResults)
}


