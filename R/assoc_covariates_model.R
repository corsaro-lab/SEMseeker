#' Prepare the covariates of one request
#'
#' Turns the covariate fields of an \code{inference_detail} into columns that
#' exist in \code{study_summary} and a vector of names the formula can use. Four
#' things happen here, in this order, and each one can rename a covariate.
#'
#' @section Every transformation renames its column, and that is the point:
#' A transformed covariate is written as a NEW column, \code{<NAME>_<SUFFIX>},
#' and the name in the returned \code{covariates} vector is replaced with it. The
#' model is therefore fitted on the new column and its coefficient is reported
#' under the name of the quantity it belongs to.
#'
#' This is not bookkeeping. A covariate transformed in place would be fitted
#' under its original name, so \code{AGE_ESTIMATE} would be the coefficient of
#' \code{exp(AGE)} - the same defect that had a polynomial estimate carrying the
#' literal string \code{INDEPENDENT_VARIABLE} instead of the variable's name.
#' \code{scale} already worked this way, with \code{_SCALED}; the rest follow it.
#'
#' @section The four steps:
#' \describe{
#'   \item{columns from a previous run}{Removed first, so running twice on the
#'     same \code{study_summary} is the same as running once. The two original
#'     patterns are unanchored as they were; the suffixes added later are anchored
#'     at the end, because a covariate name can contain one of those words and
#'     only the tail is a suffix this function wrote.}
#'   \item{\code{transformation_x == "scale"}}{Scales the INDEPENDENT variable
#'     and renames it. Through 0.99.5 it also scaled every covariate, which is a
#'     transformation asked for one variable and applied to several: a request for
#'     \code{scale} on the predictor scaled age and body mass index with it. It
#'     no longer reaches the covariates.}
#'   \item{\code{covariates_transformation}}{One transformation per covariate,
#'     paired with \code{covariates} by POSITION. The lengths have to agree and a
#'     mismatch is refused naming both, because R recycles a short vector without
#'     saying so and would apply one covariate's transformation to another's
#'     values. The vocabulary is \code{io_transform_apply()}'s and the suffix is
#'     \code{io_transform_suffix()}'s. A covariate the study summary does not
#'     carry is skipped with a warning rather than stopping the run, because a
#'     samples filter can legitimately remove the column.}
#'   \item{\code{covariates_dummy}, then \code{covariates_pca}}{Categorical
#'     encoding and then, optionally, principal components. \code{"factor"} is
#'     deliberately NOT a covariate transformation: a categorical covariate goes
#'     through \code{covariates_dummy}, which encodes it, where relabelling it as
#'     a factor does not.}
#' }
#'
#' @param inference_detail One row of the request specification.
#' @param study_summary The per-sample table the covariates live in.
#' @return A list with \code{covariates}, the possibly renamed names,
#'   \code{study_summary} carrying the new columns, and the possibly renamed
#'   \code{inference_detail}.
#'
#' @keywords internal
#' @noRd
assoc_covariates_model <- function(inference_detail, study_summary)
{
  ssEnv <- core_get_session_info()
  collinearity_check <- util_boolean_check(inference_detail$collinearity_check)
  covariates_dummy <- util_split_and_clean(inference_detail$covariates_dummy)
  covariates_pca <- util_boolean_check(inference_detail$covariates_pca)
  covariates <- util_split_and_clean(inference_detail$covariates)
  # Positional: paired with covariates one to one, so duplicates are kept.
  covariates_transformation <- util_split_and_clean(
    inference_detail$covariates_transformation, unique_values = FALSE)
  independent_variable <- as.character(inference_detail$independent_variable)
  transformation_x <- as.character(inference_detail$transformation_x)

  # Columns this function made on a previous run, removed so that running twice
  # on the same study_summary is the same as running once. The two original
  # patterns are left unanchored, as they were; the suffixes added later are
  # anchored at the end, because a covariate name can carry one of those words
  # inside it and only the tail is a suffix this function wrote.
  study_summary <- study_summary[, !grepl(
    "_SCALED|_DUMMY|_LOG$|_LOG2$|_LOG10$|_EXP$|_POW[0-9.]+$|_QUANTILE[0-9]+$|_FACTOR$",
    colnames(study_summary))]

  prev_columns <- colnames(study_summary)

  # transformation_x applies to the independent variable and to nothing else.
  # It used to scale every covariate as well, which is a transformation asked
  # for one variable and applied to several: asking for log10 of x would have
  # meant log10 of age and of body mass index too, which is not a request anyone
  # would write on purpose. A covariate now carries its own transformation, in
  # covariates_transformation, one per covariate.
  if(transformation_x=="scale")
  {
    study_summary[,paste0(independent_variable,"_SCALED")] <- scale(study_summary[,independent_variable], center = TRUE, scale = TRUE)
    core_log_event("JOURNAL: Scaling and centering applied on independent variable: ", independent_variable)
    inference_detail$independent_variable <- paste0(inference_detail$independent_variable,"_SCALED")
  }

  # One transformation per covariate, positionally paired with covariates.
  if (length(covariates_transformation) > 0L && any(nzchar(covariates_transformation)))
  {
    # R recycles a short vector without saying so, which would apply one
    # covariate's transformation to another's values and name the result for the
    # column it came from. The two lengths have to agree, and the refusal names
    # both so the request can be corrected without guessing which is short.
    if (length(covariates_transformation) != length(covariates))
      stop("covariates_transformation has ", length(covariates_transformation),
           " entries and covariates has ", length(covariates), ".\n",
           "  covariates:                ", paste(covariates, collapse = ", "), "\n",
           "  covariates_transformation: ", paste(covariates_transformation, collapse = ", "), "\n",
           "  One per covariate, in the same order, \"none\" for a covariate that ",
           "is used as it is.", call. = FALSE)

    for (cc in seq_along(covariates))
    {
      suffix <- io_transform_suffix(covariates_transformation[cc])
      if (is.na(suffix)) next

      cname <- covariates[cc]
      if (!cname %in% colnames(study_summary))
      {
        core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                  " covariates_transformation names ", covariates_transformation[cc],
                  " for ", cname, ", which the study summary does not carry: skipped.")
        next
      }

      new_name <- paste0(cname, "_", suffix)
      study_summary[, new_name] <- io_transform_apply(study_summary[, cname],
                                                      covariates_transformation[cc])
      # The name moves with the values. The model is fitted on the new column, so
      # its coefficient is reported as <COVARIATE>_<TRANSFORMATION>_ESTIMATE and
      # says which quantity it belongs to.
      covariates[cc] <- new_name
      core_log_event("JOURNAL: ", covariates_transformation[cc],
                " applied on covariate ", cname, " as ", new_name)
    }
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
