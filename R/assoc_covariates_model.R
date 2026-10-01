assoc_covariates_model <- function(inference_detail, study_summary)
{
  ssEnv <- core_get_session_info()
  collinearity_check <- util_boolean_check(inference_detail$collinearity_check)
  covariates_dummy <- util_split_and_clean(inference_detail$covariates_dummy)
  covariates_pca <- util_boolean_check(inference_detail$covariates_pca)
  covariates <- util_split_and_clean(inference_detail$covariates)
  independent_variable <- as.character(inference_detail$independent_variable)
  transformation_x <- as.character(inference_detail$transformation_x)

  # remove all "_SCALED" and "_DUMMY" from study_summary
  study_summary <- study_summary[, !grepl("_SCALED|_DUMMY", colnames(study_summary))]

  prev_columns <- colnames(study_summary)

  if(transformation_x=="scale")
  {
    scaled_cov <- c()
    for(cc in seq_along(covariates))
    {
      cname <- covariates[cc]
      if(is.numeric(study_summary[,cname]))
      {
        study_summary[,paste0(cname,"_SCALED")] <- scale(study_summary[,cname], center = TRUE, scale = TRUE)
        scaled_cov <- c(scaled_cov, cname)
      }
      # replace convariate name
      covariates[cc] <- paste0(cname,"_SCALED")
    }
    core_log_event("JOURNAL: Scaling and centering applied on covariate: ", scaled_cov)
    study_summary[,paste0(independent_variable,"_SCALED")] <- scale(study_summary[,independent_variable], center = TRUE, scale = TRUE)
    core_log_event("JOURNAL: Scaling and centering applied on independent variable: ", independent_variable)
    inference_detail$independent_variable <- paste0(inference_detail$independent_variable,"_SCALED")
  }

  # dummify covariates dummy
  if(length(covariates_dummy) > 0)
    for(i in seq_along(covariates_dummy))
    {
      encoded_covariate <- NULL
      covariate_dummy <- as.character(covariates_dummy[i])
      encoded_covariate <- fastDummies::dummy_cols(study_summary, select_columns = covariate_dummy, remove_first_dummy = TRUE)
      encoded_covariate <- encoded_covariate[, !(colnames(encoded_covariate) %in% colnames(study_summary))]
      # AI-068: when the dummy column is constant within the subset filtered
      # by samples_sql_condition (e.g. Tumour_Locus is always 'Breast' once
      # samples are restricted to Tissue=='Breast'), dummy_cols + remove_first
      # leaves zero columns. Without this guard the colnames<- below fails with
      # 'names attribute [N] must be the same length as the vector [0]'.
      n_enc <- if (is.null(dim(encoded_covariate))) length(encoded_covariate) else ncol(encoded_covariate)
      if (is.null(n_enc) || n_enc == 0L)
      {
        core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                  " dummy expansion of '", covariate_dummy,
                  "' produced 0 columns (covariate constant within sample subset?) - skipping.")
        next
      }
      if(is.null(dim(encoded_covariate)))
      {
        encoded_covariate <- data.frame(encoded_covariate)
        colnames(encoded_covariate) <- paste0(covariate_dummy,"_DUMMY")
      }
      else
      {
        encoded_covariate <- encoded_covariate[,!(colnames(encoded_covariate) %in% c("X"))]
        encoded_columns <- core_name_cleaning(colnames(encoded_covariate))
        colnames(encoded_covariate) <- encoded_columns
      }
      for (cc in seq_along(colnames(encoded_covariate)))
      {
        cname <- colnames(encoded_covariate)[cc]
        prev_columns <- prev_columns[!grepl(paste0("^",cname), prev_columns)]
        # remove column starting with the same name of encoded covariate
        study_summary <- study_summary[, !grepl(paste0("^",cname), colnames(study_summary))]
        study_summary[,cname] <- as.numeric(encoded_covariate[,cname])
      }
      covariates <- unique(c(covariates, colnames(encoded_covariate)))
    }


  if(covariates_pca)
    # check if covariates are numeric
    if(length(covariates)  !=  0)
    {
      # check if all covariates are numeric
      if(!all(vapply(study_summary[,covariates], is.numeric, logical(1))))
        core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " Not all covariates are numeric! Skipped.")


      zero_var_cols <- vapply(study_summary[, covariates], function(x) sd(x, na.rm = TRUE) == 0, logical(1))
      # Keep only non-constant columns
      filtered_covariates <- covariates[!zero_var_cols]

      # AI-070: PCA on < 3 covariates is pointless - 1-col PCA = identity,
      # 2-col PCA = orthogonal rotation that loses interpretability without
      # reducing dimensionality. Skip PCA in those cases and use the raw
      # (non-constant) covariates directly. This also subsumes the AI-069
      # Kaiser-Guttman edge case (sdev^2 == 1 on single scaled dummy).
      if (length(filtered_covariates) < 3L) {
        core_log_event("JOURNAL: PCA skipped - only ", length(filtered_covariates),
                  " non-constant covariate(s), using raw: ",
                  paste(filtered_covariates, collapse = ", "))
        covariates <- filtered_covariates
      } else {
      # Then run PCA
      pca_result <- prcomp(study_summary[, filtered_covariates], center = TRUE, scale. = TRUE)

      core_log_event("JOURNAL: PCA,scaling and centering, applied on covariates: ", paste(covariates, collapse = ", "))
      # preserve components with eigenvalue (sdev^2) above 1 - Kaiser-Guttman
      # criterion. If NO PC passes this filter (e.g. when only a single dummy
      # is left after subset filtering: scale=TRUE forces sdev^2 == 1 exactly,
      # which fails the strict '> 1' inequality) we fall back to keeping the
      # first PC so the design matrix has at least one degree of freedom.
      keep <- which(pca_result$sdev^2 > 1)
      if (length(keep) == 0L) {
        core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                  " No PC passed Kaiser-Guttman (sdev^2 > 1) filter - keeping PC1 as fallback.")
        keep <- 1L
      }
      pca_result <- pca_result$x[, keep, drop = FALSE]
      pca_result <- as.data.frame(pca_result)
      colnames(pca_result) <- paste0("PC", seq_len(ncol(pca_result)))
      covariates <- colnames(pca_result)
      for(cc in seq_along(colnames(pca_result)))
      {
        study_summary <- study_summary[, !grepl(paste0("^",colnames(pca_result)[cc],"."), colnames(study_summary))]
        cname <- colnames(pca_result)[cc]
        prev_columns <- prev_columns[!grepl(paste0("^",cname), prev_columns)]
        study_summary[,cname] <- pca_result[,cname]
      }
      }  # close AI-070 else (PCA branch)
    }

  covariates_to_remove <- c()
  if(collinearity_check && length(covariates) > 0)
    covariates_to_remove <- assoc_calculate_collinearity_score(study_summary[,c(independent_variable, covariates)])

  # preserve independent variable
  covariates_to_remove <- covariates_to_remove[!covariates_to_remove %in% inference_detail$independent_variable]

  if(length(covariates_to_remove) > 0)
  {
    core_log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), " The following covariates are collinear and will be removed: ", paste(covariates_to_remove, collapse = ", "))
    core_log_event("JOURNAL: The following covariates are collinear and will be removed: ", paste(covariates_to_remove, collapse = ", "))
    covariates <- setdiff(covariates, covariates_to_remove)
    if(length(covariates) == 0)
    {
      inference_detail <- t(as.data.frame(inference_detail))
      core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " No covariates left after collinearity check. Please check your data. ",
        "Inference detail: ", paste(inference_detail, collapse = ", "))
      stop("No covariates left after collinearity check. Please check your data.")
    }
  }

  #transform in factor not numeric covariate
  if(length(covariates)>0)
    for(cc in seq_along(covariates))
    {
      cname <- covariates[cc]
      if(!is.numeric(stats::na.omit(study_summary[,cname])))
      {
        study_summary[,cname] <- as.factor(study_summary[,cname])
      }
    }


  post_columns <- colnames(study_summary)
  columns_to_save <- post_columns[!post_columns %in% prev_columns]
  inf_file_name <- io_inference_file_name(inference_detail,"",ssEnv$result_folderInference,"csv",  suffix = "covariates_model_")
  write.csv2(study_summary[,c("Sample_ID",columns_to_save)], file = inf_file_name, row.names = FALSE)
  # file_path <- io_file_path_build(ssEnv$result_folderData, "sample_sheet_result" ,"csv")
  # write.csv2(study_summary, file_path, row.names = FALSE)

  return (list(covariates = covariates, study_summary = study_summary, inference_detail = inference_detail))

}
