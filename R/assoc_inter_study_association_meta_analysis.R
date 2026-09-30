#' Cross-study meta-analysis of association results
#'
#' Combines inference results from multiple studies using a random-effects
#' meta-analysis model (\code{\link[meta]{metagen}}).  For each unique
#' combination of FIGURE, SUBAREA and AREA_OF_TEST, the function pools
#' effect sizes (BETA) and standard errors across studies and reports
#' fixed-effect and random-effect estimates with heterogeneity statistics.
#'
#' Requires at least two studies per stratum; strata with fewer studies are
#' silently skipped.
#'
#' @param inference_details \code{data.frame} describing the inference
#'   configuration (same format as used by \code{association_analysis}).
#' @param statistic_parameter Character scalar: column name of the effect
#'   size estimate in the inference results (default \code{"BETA"}).
#' @param pvalue_column Character scalar: column name of the adjusted p-value
#'   (default \code{"PVALUE_ADJ_ALL_BH"}).
#' @param studies Character vector: study identifiers to include.
#' @param studies_base_folder Character scalar: base directory containing
#'   per-study result folders.
#' @param result_folder Character scalar: output directory for the
#'   meta-analysis results.
#'
#' @return Invisibly returns a \code{data.frame} with one row per stratum
#'   containing pooled effect estimates, confidence intervals, p-values, and
#'   heterogeneity statistics (\eqn{\tau^2}, Q-test p-value).
#'
#' @export
#' @examples
#' # Stub: see vignette('imprinting-disorders', package = 'SEMseeker') for a
#' # runnable Beckwith-Wiedemann workflow on the GSE133774 subset.
#' invisible(NULL)
assoc_inter_study_association_meta_analysis <- function(inference_details,statistic_parameter="BETA", pvalue_column="PVALUE_ADJ_ALL_BH",studies,
  studies_base_folder, result_folder)
{

  if (!requireNamespace("meta", quietly = TRUE))
    stop("Package 'meta' is required for cross-study meta-analysis. Install it with install.packages('meta').")
  if (!requireNamespace("tidyverse", quietly = TRUE))
    stop("Package 'tidyverse' is required for cross-study meta-analysis. Install it with install.packages('tidyverse').")
  keys <- na.omit(unique(results_inference[,c("FIGURE","SUBAREA")]))
  for (k in seq_len(nrow(keys)))
  {
    # k <- 1
    results_inference_for <- subset(results_inference, FIGURE==keys[k,"FIGURE"] & SUBAREA==keys[k,"SUBAREA"])
    results_inference_for <- na.omit(results_inference_for[,c("BETA","STD.ERROR","STUDY",pvalue_column,"AREA_OF_TEST")])
    areas <- na.omit(unique(results_inference[, "AREA_OF_TEST"]))
    for (g in seq_along(areas))
    {
      # g <- 1
      first_area <- areas[g]
      if(is.na(first_area))

      results_inference_subset <- subset(results_inference_for, AREA_OF_TEST==first_area)
      # if( is.null(results_inference_subset[,"STUDY"]))
      #   next
      studies_count <- length(unique(results_inference_subset[,"STUDY"]))
      if(studies_count<2)
        next
      meta_model <- meta::metagen(
        TE = results_inference_subset[,statistic_parameter],   # effect size for each study
        seTE = results_inference_subset$STD.ERROR,  # standard error of the effect size for each study
        studlab = results_inference_subset$STUDY,  # study label
        data = results_inference_subset,  # data frame containing all data
        comb.fixed = FALSE, # use random-effects model
        hakn = FALSE, # do not apply Hartung-Knapp correction
        TE.targ = 1, # target effect size is the regression coefficient
        pval = results_inference_subset[,pvalue_column],  # p-value column in the data frame
        ncpus = 9
      )
      meta_analysis_results <- summary(meta_model)
      # pval <- meta_analysis_results$pval
      beta <- meta_analysis_results$random$TE

      common.ci.lower <- meta_analysis_results$common$lower
      common.ci.upper <- meta_analysis_results$common$upper
      common.pval <- meta_analysis_results$pval.common

      random.ci.lower <- meta_analysis_results$random$lower
      random.ci.upper <- meta_analysis_results$random$upper
      random.pval <- meta_analysis_results$pval.random

      # test of heterogeneity
      pval.Q <- meta_analysis_results$pval.Q
      pval.fixed <- meta_analysis_results$pval.fixed
      k.study <- meta_analysis_results$k.study
      tau2 <- meta_analysis_results$tau2

      studies <- paste(sort(meta_analysis_results$studlab), collapse = " ")
      result <- data.frame(tau2,beta,common.ci.lower,common.ci.upper, common.pval, random.ci.lower,random.ci.upper, random.pval, pval.Q, pval.fixed,k.study, studies, first_area, studies_count)
      if (exists("final_result"))
        final_result <- rbind(final_result, result)
      else
        final_result <- result
    }
  }
}
