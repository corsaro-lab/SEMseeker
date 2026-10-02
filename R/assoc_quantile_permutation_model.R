#' Compute quantile delta permutation (CPU)
#'
#' @param sig.formula formula to apply
#' @param df dataframe to use
#' @param shuffle logical; if TRUE, permute the independent variable before fitting
#' @param quantile quantile level (0-1) at which to compute the group difference; default 0.5 (median)
#'
#' @return A numeric scalar: the observed quantile difference between groups
#'   (used as the test statistic in the permutation test).
#'
assoc_compute_quantile_delta_permutation <- function(sig.formula,df, shuffle = FALSE, quantile = 0.5)
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
  # calculate quantile difference based on indepVar
  quantile_difference <- tapply(tempDataFrame[,burden], tempDataFrame[,indepVar], function(x) stats::quantile(x, probs=quantile))
  statistic_parameter <- diff(quantile_difference)
  statistic_parameter <- as.numeric(statistic_parameter[[1]])
  return(statistic_parameter)
}

#' assoc_quantile_permutation_model calculate differences between the same quantile of two distribution
#'
#' @param family_test family quantile
#' @param sig.formula formula of the model
#' @param tempDataFrame data
#' @param independent_variable name of regressor
#' @param transformation_y transformation to apply to the dependent variable
#' @param plot logical; if TRUE, generate diagnostic plots
#' @param samples_sql_condition SQL condition string used to filter samples (used for plot file naming)
#' @param key named list with AREA, SUBAREA, MARKER and FIGURE identifiers for this test
#'
#' @return A numeric p-value from the permutation-based quantile-difference test.
#'
assoc_quantile_permutation_model <- function(family_test, sig.formula, tempDataFrame, independent_variable, transformation_y, plot =plot, samples_sql_condition=samples_sql_condition, key)
{
  quantile_params <- unlist(strsplit(as.character(family_test),"_"))
  area <- as.character(key$AREA)
  subarea <- as.character(key$SUBAREA)
  marker <- as.character(key$MARKER)
  figure <- as.character(key$FIGURE)

  # quantile_params template quantile + first_round_of_permutations + second_round_of_permutations + confidence_interval_of_beta
  # apply permutation to obtain signal using quantile
  # Define function to compute delta quantile regression coefficient

  n_permutations_test <- as.numeric(quantile_params[2])
  res <- data.frame(n_permutations = n_permutations_test)

  n_permutations <- as.numeric(quantile_params[3])
  conf.level <- as.numeric(quantile_params[4])
  res$conf.level <- conf.level
  tau <- as.numeric(quantile_params[5])
  res$tau <- tau

  res$independent_variable <- as.character(all.vars(sig.formula)[2])

  pvalue_limit <- 1 - conf.level
  pvalue_limit_inf <- (pvalue_limit/2)
  pvalue_limit_sup <- 1 - (pvalue_limit/2)

  # Compute signal and p-value for n_permutations replications
  # The statistic and the p-value are kept in local variables and then
  # written to the result columns. They used to be assigned straight into
  # res$... and read back as bare names three lines later, which are two
  # different things: the bare names never existed and the function died
  # here every time it was called.
  statistic_parameter <- assoc_compute_quantile_delta_permutation(sig.formula = sig.formula, df = tempDataFrame, shuffle = FALSE)
  res$statistic_parameter <- statistic_parameter
  permutation_vector <- replicate(n_permutations_test, assoc_compute_quantile_delta_permutation(sig.formula=sig.formula, df=tempDataFrame, shuffle = TRUE, quantile = tau))
  pvalue <- mean(abs(permutation_vector) >= abs(statistic_parameter))
  if (pvalue > 1)
    pvalue <- 1
  res$pvalue <- pvalue
  # Compute average signal and p-value
  if ((pvalue < pvalue_limit) && (n_permutations_test < n_permutations))
    permutation_vector <- replicate(n_permutations, assoc_compute_quantile_delta_permutation(sig.formula=sig.formula, df=tempDataFrame, shuffle=TRUE, quantile = tau))
  res$r_model <- "assoc_quantile_permutation_model"
  if(length(permutation_vector) == n_permutations_test)
    res$n_permutations <- n_permutations_test
  pvalue <- mean(abs(permutation_vector) >= abs(statistic_parameter))
  res$pvalue <- pvalue
  ci <- stats::quantile(permutation_vector, probs = c(pvalue_limit_inf, pvalue_limit_sup))
  res$ci.lower <- ci[1]
  res$ci.upper <- ci[2]

  # independent_variable <- as.character(all.vars(sig.formula)[2])
  # dependent_variable <- as.character(all.vars(sig.formula)[1])

  # predicted_values <- lqmm::predict.lqm(model, tempDataFrame[,independent_variable])
  # expected_values <- quantile(tempDataFrame[,independent_variable], tau)
  # res <- assoc_quantreg_metrics(predicted_values = predicted_values, expected_values = expected_values,
  #   tau = tau, res = res, family_test = family_test, independent_variable = independent_variable, transformation_y = transformation_y,
  #   dependent_variable = dependent_variable, permutation_vector = permutation_vector, plot = plot)

  return(res)
}
