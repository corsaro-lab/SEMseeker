#' Mediation analysis using linear models
#'
#' Fits a causal mediation model via \code{\link[mediation]{mediate}}, testing
#' whether the effect of \code{treatment} on \code{outcome} methods::is (partially)
#' mediated by a \code{mediator} variable.  The formula must follow the
#' convention \code{mediator ~ outcome + treatment [+ covariates]}.
#'
#' The two-round permutation scheme mirrors the approach used in
#' \code{assoc_mean_permutation} and \code{assoc_spearman_permutation}: a fast first round
#' (\code{permutations_test}) establishes significance; if significant, a
#' slower second round (\code{permutations}) refines the estimate.
#'
#' @param family_test Character string encoding the model type and permutation
#'   counts, formatted as \code{"mediation_<permutations_test>_<permutations>"}.
#' @param tempDataFrame \code{data.frame} containing all variables referenced
#'   in \code{sig.formula}.
#' @param sig.formula Formula of the form
#'   \code{mediator ~ outcome + treatment [+ covariates]}.
#' @param transformation_y Character scalar: transformation applied to the
#'   dependent variable (e.g. \code{"none"}, \code{"log"}).
#' @param plot Logical; if \code{TRUE}, generate diagnostic plots.
#' @param samples_sql_condition Character scalar: SQL \code{WHERE} clause used
#'   to subset samples (propagated to file-naming helpers).
#' @param key Named list with elements \code{AREA}, \code{SUBAREA},
#'   \code{MARKER}, and \code{FIGURE} identifying this test instance.
#'
#' @return A \code{data.frame} with one row containing mediation results:
#'   ACME (average causal mediation effect), ADE (average direct effect),
#'   total effect, proportion mediated, their confidence intervals and
#'   p-values, plus the number of permutations used.
#'
# (family_test, tempDataFrame, sig.formula, transformation_y, plot)
assoc_mediation_linear_model <- function(family_test,tempDataFrame, sig.formula, transformation_y, plot, samples_sql_condition=samples_sql_condition, key)
{

  area <- as.character(key$AREA)
  subarea <- as.character(key$SUBAREA)
  marker <- as.character(key$MARKER)
  figure <- as.character(key$FIGURE)

  model_params <- unlist(strsplit(as.character(family_test),"_"))
  permutations_test <- as.numeric(model_params[2])
  res <- data.frame(permutations_test = permutations_test)
  permutations <- as.numeric(model_params[3])
  res$permutations <- permutations

  # decompose the formula, to move the mediator
  # the expected formula is "mediator ~ outcome + treatment + covariates"
  vars <- assoc_sig_formula_vars(sig.formula)
  mediator <- vars[[1]] # burdnValue
  outcome <- vars[[2]] # dependent variable
  treatment <- vars[[3]][1]

  res$mediator <- mediator
  res$outcome <- outcome
  res$treatment <- treatment

  if(length(vars[[3]])>1)
  {
    covariates <- vars[[3]][-c(1)]
    # transform covariates in a vector of string
    covariates <- vapply(covariates, as.character, character(1))
  } else
    covariates <- c()

  df_mediator <- data.frame(tempDataFrame)
  df_mediator$mediator <- df_mediator[,mediator]
  df_mediator$treatment <- df_mediator[,treatment]
  if(length(covariates) != 0)
    df_mediator <- df_mediator[,c("treatment", "mediator",covariates)]
  else
    df_mediator <- df_mediator[,c("treatment", "mediator")]

  df_outcome <- data.frame(data)
  df_outcome$outcome <- df_outcome[,outcome]
  df_outcome$mediator <- df_outcome[,mediator]
  df_outcome$treatment <- df_outcome[,treatment]
  if(length(covariates) != 0)
    df_outcome <- df_outcome[,c("outcome", "mediator", "treatment",covariates)]
  else
    df_outcome <- df_outcome[,c("outcome", "mediator", "treatment")]

  perms <- sort(unique(c(permutations_test, permutations)))

  for (i in seq_along(perms)) {

    env <- new.env()
    env$mediate_env <- mediation::mediate
    perm <- perms[i]
    mediation_result <- eval(expr = quote(mediate_env(
      model.m = lm(as.formula(paste("mediator", "~", ".")), data = df_mediator),
      model.y = lm(as.formula(paste("outcome", "~", ".")), data = df_outcome),
      treat = "treatment",
      mediator = "mediator",
      boot = TRUE,
      sims = perm
    )), envir = env)
    res$pvalue <- max(mediation_result[1,4], mediation_result[2,4], mediation_result[3,4], mediation_result[4,4])
    res$permutations <- perm
    # Extract key components from the summary
    res$acme <- mediation_result[1,1]
    res$acme_ci_start <- mediation_result[1,2]
    res$acme_ci_end <- mediation_result[1,3]
    res$acme_pvalue <- mediation_result[1,4]

    res$ade <- mediation_result[2,1]
    res$ade_ci_start <- mediation_result[2,2]
    res$ade_ci_end <- mediation_result[2,3]
    res$ade_pvalue <- mediation_result[2,4]

    res$total_effect <- mediation_result[3,1]
    res$total_effect_ci_start <- mediation_result[3,2]
    res$total_effect_ci_end <- mediation_result[3,3]
    res$total_effect_pvalue <- mediation_result[3,4]

    res$prop_mediated <- mediation_result[4,1]
    res$prop_mediated_ci_start <- mediation_result[4,2]
    res$prop_mediated_ci_end <- mediation_result[4,3]
    res$r_model <- family_test
    if(res$pvalue < as.numeric(ssEnv$alpha))
      next
  }

  # remove environment
  rm(env)


  # Return the summary of the mediation analysis
  return(res)
}

