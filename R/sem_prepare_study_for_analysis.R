#' Prepare study_summary, covariates, and sample_names for one inference job
#'
#' Extracted from association_analysis() (was inline at lines 143-214).
#' Wraps assoc_covariates_model() + independent_variable validation + factor
#' conversion for dichotomous families + sample_names construction with
#' complete_cases filtering.
#'
#' @param inference_detail single-row data.frame.
#' @param study_summary data.frame returned by sem_study_summary_get().
#' @param family_test character. Already util_split_and_clean'd.
#' @return list(study_summary, covariates, sample_names, independent_variable,
#'   transformation_y, inference_detail, file_result_prefix)
#'   OR NULL if the job must be skipped (logged via core_log_event).
#' @keywords internal
sem_prepare_study_for_analysis <- function(inference_detail, study_summary, family_test) {
  transformation_y <- inference_detail$transformation_y
  if (is.null(transformation_y) || length(transformation_y) == 0)
    transformation_y <- NULL

  res_model_covariates <- assoc_covariates_model(inference_detail, study_summary)
  study_summary    <- res_model_covariates$study_summary
  covariates       <- res_model_covariates$covariates
  inference_detail <- res_model_covariates$inference_detail

  independent_variable <- gsub(" ", "", inference_detail$independent_variable)

  if (independent_variable %in% covariates) {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
      " The independent variable is also present as covariate!")
    return(NULL)
  }

  if (is.null(independent_variable) || length(independent_variable) == 0) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
      " One indipendent variable is missed! Skipped.")
    return(NULL)
  }

  # transform independent variable as factor for dichotomous families
  if (family_test == "binomial" || family_test == "wilcoxon" || family_test == "t.test")
    study_summary[, independent_variable] <- as.factor(study_summary[, independent_variable])

  file_result_prefix <- as.character(independent_variable)

  if (!(independent_variable %in% colnames(study_summary))) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
      " This indipendent variabile:", independent_variable, " is missed! Skipping")
    return(NULL)
  }

  # build sample_names with optional covariates
  has_covariates <- !is.null(covariates) && length(covariates) != 0
  if (has_covariates) {
    sample_names <- data.frame(study_summary[, c("Sample_ID", independent_variable, covariates)])
    sample_names <- sample_names[, c("Sample_ID", independent_variable, covariates)]
    colnames(sample_names) <- c("Sample_ID", independent_variable, covariates)
  } else {
    sample_names <- data.frame(study_summary[, c("Sample_ID", independent_variable)])
    sample_names <- unique(sample_names[, c("Sample_ID", independent_variable)])
    colnames(sample_names) <- c("Sample_ID", independent_variable)
  }

  # remove samples with missing values
  sample_names <- sample_names[complete.cases(sample_names), ]
  if (nrow(sample_names) == 0) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
      " No samples with complete data for the analysis! Skipped.")
    return(NULL)
  }

  list(
    study_summary        = study_summary,
    covariates           = covariates,
    sample_names         = sample_names,
    independent_variable = independent_variable,
    transformation_y     = transformation_y,
    inference_detail     = inference_detail,
    file_result_prefix   = file_result_prefix
  )
}
