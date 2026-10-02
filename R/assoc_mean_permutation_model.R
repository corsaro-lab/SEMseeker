#' Compute mean delta permutation (CPU)
#'
#' @param sig.formula formula to apply
#' @param df dataframe to use
#' @param shuffle logical; if TRUE, permute the independent variable before fitting
#'
#' @return A numeric scalar: the observed mean difference between the two groups
#'   (used as the test statistic in the permutation test).
#'
assoc_compute_mean_delta_permutation_cpu <- function(sig.formula,df, shuffle = FALSE)
{
  # #
  cols <- colnames(df)
  indepVar <- as.character(all.vars(sig.formula)[2])
  burden <- as.character(all.vars(sig.formula)[1])
  tempDataFrame <- df
  idx <- which(colnames(tempDataFrame) == indepVar)
  tempDataFrame[,idx] <- as.factor(tempDataFrame[,idx])
  if (shuffle == TRUE)
    tempDataFrame[ ,indepVar] <- sample(tempDataFrame[,indepVar])
  tempDataFrame <- as.data.frame(tempDataFrame)
  colnames(tempDataFrame) <- cols
  # calculate mean difference based on indepVar
  mean_difference <- tapply(tempDataFrame[,burden], tempDataFrame[,indepVar], mean)
  # statistic_parameter <- mean(tempDataFrame[tempDataFrame[,indepVar]==1,burden]) - mean(tempDataFrame[tempDataFrame[,indepVar]==0,burden])
  statistic_parameter <- diff(mean_difference)
  statistic_parameter <- -1 * as.numeric(statistic_parameter[[1]])
  return(statistic_parameter)
}

assoc_compute_mean_delta_permutation_gpu <- function(sig.formula,df, n_permutations, shuffle=FALSE)
{
  ssEnv <- core_get_session_info()

  # burden_var <- as.character(all.vars(sig.formula)[1])
  independent_var <- as.character(all.vars(sig.formula)[2])
  burden <- df[,-which(colnames(df) %in% c(independent_var))]
  independent_var <- df[,independent_var]
  n_samples <- length(burden)
  n_areas <- ncol(burden)

  if (is.null(ssEnv$independent_var_permutated))
  {
    # create a independent_var_permutated of n_permutations x n_samples where each n_permutation is a permutation of independent_var
    independent_var_permutated <- matrix(rep(independent_var, n_permutations), nrow = n_permutations, byrow = TRUE)
    independent_var_permutated <- t(apply(independent_var_permutated, 1, function(x) sample(x, length(x))))
    ssEnv$independent_var_permutated <- independent_var_permutated
  }
  else
  {
    independent_var_permutated <- ssEnv$independent_var_permutated
  }
  # core_log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y"), " Permuted independent_var!")

  context <- OpenCL::oclContext(precision = "single", device="gpu")
  # core_log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y"), " OpenCL Context Created")

  compute_permutation_diff_kernel_obj <- OpenCL::oclSimpleKernel(context, "compute_row_means_diff", compute_permutation_diff_kernel, "single")

  ocl_run <- function(compute_permutation_diff_kernel_obj, n_permutations,burden, independent_var_permutated,n_samples)
    as.numeric(
      OpenCL::oclRun(
        kernel = compute_permutation_diff_kernel_obj,
        dim = c(n_permutations * n_areas), # First argument: number of n_permutations (vector)
        as.integer(n_permutations), # Second argument: number of n_permutations (scalar)
        as.integer(n_samples), # Fifth argument: number of columns (scalar)
        OpenCL::as.clBuffer(as.single(burden), context, mode = "single"), # Third argument: buffer for the matrix
        OpenCL::as.clBuffer(as.integer(independent_var_permutated), context, mode = "int") # Fourth argument: buffer for the phenotype vector
      )
    )

  # measure time
  # system.time(res <- ocl_run(compute_permutation_diff_kernel_obj_ok, n_permutations,burden_permutated_buf, phenotype,n_samples))
  output_diffs_buf <- OpenCL::as.clBuffer(as.single(rep(0,n_permutations)), context, mode = "single")
  res <- ocl_run(compute_permutation_diff_kernel_obj, n_permutations,burden, independent_var_permutated,n_samples)
  # core_log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y"), " Execution OPENCL done !")
  # system.time(res_ok <- apply(burden_permutated, 1, function(x) mean(x[phenotype == 0]) - mean(x[phenotype == 1])))
  # all.equal(res, res_ok)
  #
  return(res)
}


#' Mean permutation model
#'
#' @param family_test family test string encoding model type and permutation parameters
#' @param sig.formula formula of the model
#' @param tempDataFrame data
#' @param independent_variable name of regressor
#' @param plot logical; if TRUE, generate diagnostic plots
#' @param samples_sql_condition SQL condition string used to filter samples (used for plot file naming)
#' @param key named list with AREA, SUBAREA, MARKER and FIGURE identifiers for this test
#'
#' @return A numeric p-value from the permutation-based mean-difference test.
#'
assoc_mean_permutation <- function(family_test, sig.formula, tempDataFrame, independent_variable,plot, samples_sql_condition=samples_sql_condition, key)
{
  area <- as.character(key$AREA)
  subarea <- as.character(key$SUBAREA)
  marker <- as.character(key$MARKER)
  figure <- as.character(key$FIGURE)

  ssEnv <- core_get_session_info()
  n_permutations <- NA
  n_permutations_test <- NA
  ci.lower <- NA
  ci.upper <- NA
  aic_value <- NA
  residuals <- NA
  shapiro_pvalue <- NA
  std.error <- NA
  statistic_parameter <- NA
  pvalue <- NA

  mean_params <- unlist(strsplit(as.character(family_test),"_"))

  # mean_params template mean + first_round_of_permutations + second_round_of_permutations + confidence_interval_of_beta
  # apply permutation to obtain signal using mean
  # Define function to compute delta mean regression coefficient
  n_permutations_test <- as.numeric(mean_params[2])
  n_permutations <- as.numeric(mean_params[3])
  conf.level <- as.numeric(mean_params[4])



  # Compute signal and p-value for n_permutations replications
  statistic_parameter <-  assoc_compute_mean_delta_permutation_cpu(sig.formula=sig.formula, df=tempDataFrame, shuffle = FALSE)
  # if (ssEnv$opencl)
  #   permutation_vector <- assoc_compute_mean_delta_permutation_gpu(sig.formula=sig.formula, df=tempDataFrame, shuffle=FALSE, n_permutations=1)
  # else
  #   permutation_vector <- replicate(1, assoc_compute_mean_delta_permutation_cpu(sig.formula=sig.formula, df=tempDataFrame, shuffle=FALSE))

  if (ssEnv$opencl)
    permutation_vector <- assoc_compute_mean_delta_permutation_gpu(sig.formula=sig.formula, df=tempDataFrame, shuffle=TRUE, n_permutations=n_permutations_test)
  else
    permutation_vector <- replicate(n_permutations_test, assoc_compute_mean_delta_permutation_cpu(sig.formula=sig.formula, df=tempDataFrame, shuffle=TRUE))

  summary_results <- assoc_exact_pvalue(permutation_vector, statistic_parameter, conf.level = conf.level)
  pvalue <- summary_results[3]
  pvalue_limit <- summary_results[4]


  # Compute average signal and p-value
  if ((pvalue < pvalue_limit) && (n_permutations_test < n_permutations))
    if (ssEnv$opencl)
      permutation_vector <- assoc_compute_mean_delta_permutation_gpu(sig.formula=sig.formula, df=tempDataFrame, shuffle=TRUE, n_permutations=n_permutations)
    else
      permutation_vector <- replicate(n_permutations, assoc_compute_mean_delta_permutation_cpu(sig.formula=sig.formula, df=tempDataFrame, shuffle=TRUE))


  summary_results <- assoc_exact_pvalue(permutation_vector, statistic_parameter, conf.level = conf.level)
  pvalue <- summary_results[3]

  r_model <- "assoc_mean_permutation"
  n_permutations <- length(permutation_vector)

  ci.lower <- summary_results[1]
  ci.upper <- summary_results[2]


  return (data.frame(ci.lower,ci.upper, pvalue, statistic_parameter,aic_value,residuals,shapiro_pvalue,r_model,std.error,n_permutations))
}
