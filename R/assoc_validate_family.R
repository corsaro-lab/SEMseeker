assoc_validate_family_test <- function(family_test){

  if( is.null(family_test) || length(family_test)  ==  0 || is.na(family_test))
  {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " One test family_test is missed! Skipped.", family_test)
    return(FALSE)
  }

  core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " family_test: " , as.character(family_test))


  if(family_test=="multinomial" | family_test=="binomial" | family_test=="binomial_bulk" | family_test=="wilcoxon" | family_test=="jsd" | family_test=="t.test" | family_test=="poisson" |
      family_test=="chisq.test" | family_test=="fisher.test" | family_test=="kruskal.test" | family_test=="pearson" | family_test=="kendall" | family_test=="spearman" |
      family_test=="wilcoxon" | family_test=="gaussian")
    return(TRUE)

  if(grepl("mean-permutation",family_test))
  {
    mean_params <- unlist(strsplit(as.character(family_test),"_"))
    if (length(mean_params) != 4)
      core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), "mean-permutation family_test must have been with the follwing syntax mean-permutation_n.permutations.test_n.permutations_conf.level")
    else
      return(TRUE)
  }

  if(grepl("wilcoxon.paired",family_test))
  {
    mean_params <- unlist(strsplit(as.character(family_test),"@"))
    if (length(mean_params) != 2)
      core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), "wilcoxon.paired family_test must have been with the following syntax wilcoxon.paired@pairing_variable")
    else
      return(TRUE)
  }

  if(grepl("t.test.paired",family_test))
  {
    mean_params <- unlist(strsplit(as.character(family_test),"@"))
    if (length(mean_params) != 2)
      core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), "t.test.paired family_test must have been with the following syntax t.test.paired@pairing_variable")
    else
      return(TRUE)
  }

  if (grepl("quantile-permutation",family_test))
  {
    quantile_params <- unlist(strsplit(as.character(family_test),"_"))
    if (length(quantile_params) == 5)
      return(TRUE)
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), "quantile-permutation family_test must have been with the follwing syntax quantile-permutation_quantile_n.permutations.test_n.permutations_conf.level")
    return(FALSE)
  }

  if (grepl("quantreg-permutation",family_test))
  {
    quantile_params <- unlist(strsplit(as.character(family_test),"_"))
    if (length(quantile_params) != 5)
      core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), "quantreg-permutation family_test must have been with the follwing syntax quantile-permutation_tau_n.permutations.test_n.permutations_conf.level")
    else
      return(TRUE)
  }

  if (grepl("mediation-ridge", family_test))
    return(TRUE)

  if (grepl("mediation-linear", family_test))
    return(TRUE)

  if (grepl("spearman-permutation",family_test))
    return(TRUE)

  if (grepl("quantreg", family_test))
    return(TRUE)

  if(grepl("polynomial_",family_test))
    return(TRUE)

  # AI-040: limma_<degree>[_<partition>] and voom_<degree>[_<partition>].
  # Same parser shape as polynomial; the actual guard against missing
  # limma installation lives at the dispatch point (assoc_apply_stat_model
  # for the batch path, assoc_execute_model for the per-area path) following
  # the AI-038 dispatch=guard convention.
  if (grepl("^(limma|voom)_", family_test))
    return(TRUE)

  if (grepl("exp_",family_test))
    return(TRUE)

  if (grepl("log_",family_test))
    return(TRUE)

  if (grepl("log10_",family_test))
    return(TRUE)

  if (grepl("pow10_",family_test))
    return(TRUE)

  core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), "family_test is not recognized: ", family_test)
  return(FALSE)
}

#' Refuse a request whose family test cannot be run
#'
#' AI-309. `assoc_validate_family_test()` answers a question: is this string a
#' family this package knows how to fit, and answers it for one value. This is
#' the *consequence* of that answer at the door of `association_analysis()`,
#' alongside [assoc_validate_scope()] and [assoc_validate_aggregation()].
#'
#' It used to have no consequence worth the name. The per-row loop did
#' `if (!assoc_validate_family_test(family_test)) next`: the row vanished, the
#' log gained a line, and the run carried on to write an inference CSV that
#' looks exactly like one where every requested test had been fitted. A reader
#' of that file cannot tell a model that found nothing from a model that was
#' never run, and the whole point of naming the six coordinates in the output
#' is to make that distinction impossible to lose.
#'
#' Skipping is right for an *environmental* condition, where the request was
#' well formed and something outside it was missing. A family test that is
#' absent or is not in the vocabulary is a malformed request: no rewriting of it
#' can be trusted to mean what its author intended, so it stops the run before
#' any result exists.
#'
#' @param inference_details validated data.frame of requests.
#' @return `inference_details` unchanged; called for the refusal.
#' @keywords internal
#' @noRd
assoc_validate_family <- function(inference_details) {

  if (is.null(inference_details) || nrow(inference_details) == 0)
    return(inference_details)

  for (z in seq_len(nrow(inference_details))) {
    # Same cleaning the loop applies, so the door judges the value the loop
    # would have used and not a different one.
    family_test <- util_split_and_clean(inference_details$family_test[z])

    if (length(family_test) == 0)
      stop("inference_details row ", z, ": 'family_test' is required. Name the ",
           "model to fit: a group test, a GLM family, a correlation, a ",
           "quantile regression. The row used to be dropped with a line in the ",
           "log, which left a result file that could not be told apart from ",
           "one where the test had run and found nothing.", call. = FALSE)

    # The vocabulary lives in assoc_validate_family_test() and stays there; it
    # also logs the value it refused, which is worth keeping.
    if (!assoc_validate_family_test(family_test))
      stop("inference_details row ", z, ": family_test = '",
           paste(family_test, collapse = ", "),
           "' is not a model this package can fit. See ?association_analysis ",
           "for the families and their parameterised forms.", call. = FALSE)
  }

  inference_details
}
