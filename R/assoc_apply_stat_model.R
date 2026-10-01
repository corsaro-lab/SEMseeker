#' Title
#'
#' @param tempDataFrame data frame to apply association
#' @param g_start index of starting data
#' @param family_test family of test to run
#' @param covariates vector of covariates
#' @param key key to identify file to elaborate
#' @param transformation_y transformation_y to apply to covariates, burden and independent variable
#' @param session_folder where to save log file
#' @param independent_variable independent variable name
#' @param samples_sql_condition SQL condition string to filter samples
#' @param inference_detail one row of the inference specification, carrying the
#'   per-request fields (scope, aggregation, transformation_x) that the model
#'   call needs; NULL falls back to the defaults
#' @param ... extra parameters
#'
#' @return A data.frame with one row per tested genomic area, including columns
#'   for p-value, adjusted p-value, test statistic, AIC, residuals, and model
#'   metadata; returns \code{NULL} if no results could be computed.
#'
#' @importFrom doRNG %dorng%
#' @importFrom doFuture %dofuture%
#'
#'
assoc_apply_stat_model <- function(tempDataFrame, g_start, family_test, covariates = NULL, key, transformation_y,
  session_folder, independent_variable, samples_sql_condition,
  inference_detail = NULL, ...)
{
  arguments <- list(...)

  # AI-040 Fase 2+3: limma_<N> and voom_<N> bypass the per-area foreach.
  # limma needs the full M x N response matrix for eBayes shrinkage to
  # mean anything; voom literally cannot estimate its mean-variance
  # trend on a 1-row matrix. assoc_apply_stat_model_batch() fits once on the
  # whole chunk and returns N rows with the same schema the foreach
  # would have produced. The guard mirrors the AI-038 dispatch=guard
  # pattern: fail fast with install hint if limma is missing.
  if (grepl("^(limma|voom)_", family_test)) {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("family_test='", family_test,
           "' requires the 'limma' package. Install it with:\n",
           "  BiocManager::install('limma')")
    }
    return(assoc_apply_stat_model_batch(
      tempDataFrame        = tempDataFrame,
      g_start              = g_start,
      family_test          = family_test,
      covariates           = covariates,
      key                  = key,
      transformation_y     = transformation_y,
      session_folder       = session_folder,
      independent_variable = independent_variable,
      samples_sql_condition = samples_sql_condition,
      inference_detail     = inference_detail,
      ...
    ))
  }

  # AI-044 (2026-06-08): bulk path for logistic regression. Same
  # dispatch=guard pattern as limma/voom - guard against missing
  # Rfast lives inside assoc_glm_model_bulk(). Returns one row per probe
  # with the legacy schema (per-coef PVALUE/ESTIMATE + top-level
  # PVALUE/PVALUE_ADJ) so downstream CSV machinery doesn't change.
  if (family_test == "binomial_bulk") {
    return(assoc_glm_model_bulk(
      tempDataFrame        = tempDataFrame,
      g_start              = g_start,
      family_test          = family_test,
      covariates           = covariates,
      key                  = key,
      transformation_y     = transformation_y,
      session_folder       = session_folder,
      independent_variable = independent_variable,
      samples_sql_condition = samples_sql_condition,
        inference_detail     = inference_detail,
      ...
    ))
  }

  # Session info is only needed for the per-area (foreach) path -
  # the batch path above is intentionally session-independent so it
  # can be exercised in unit tests without a materialised session.
  ssEnv <- core_get_session_info()

  g_end <- ncol(tempDataFrame)
  transformation_x_local <- if (!is.null(inference_detail$transformation_x)) as.character(inference_detail$transformation_x) else "none"
  prepared_data <- io_data_preparation(family_test,transformation_y,tempDataFrame, independent_variable, g_start, g_end, covariates, key, transformation_x = transformation_x_local)
  # if(ncol(prepared_data$tempDataFrame) != ncol(tempDataFrame))
  #   return(NULL)

  tempDataFrame <- prepared_data$tempDataFrame
  independent_variable1stLevel <- prepared_data$independent_variableLevels[1]
  independent_variable2ndLevel <- prepared_data$independent_variableLevels[2]

  # AI-106 (2026-06-09): name memoisation. R formula identifiers cannot
  # contain '-', ':', '/', etc., but the upstream annotation may carry
  # raw names like "HLA-A" or "chr10:100028204-100028508". We sanitise
  # the colnames to a R-safe form BEFORE the foreach loop, keep a
  # safe→real mapping in scope, and reverse it when assigning
  # AREA_OF_TEST so the result CSV preserves the original names.
  # Enrichment downstream (WebGestalt / STRINGdb / pathfindR) then
  # resolves HGNC symbols correctly; resume match is exact.
  real_cols <- colnames(tempDataFrame)
  safe_cols <- gsub("[^A-Za-z0-9_.]", "_", real_cols)
  # Disambiguate if sanitisation collapses distinct real names onto the
  # same safe form (rare in practice for HGNC, but possible for noisy
  # annotations): append a numeric suffix to the duplicates.
  if (anyDuplicated(safe_cols)) {
    safe_cols <- make.unique(safe_cols, sep = "_")
  }
  safe_to_real <- stats::setNames(real_cols, safe_cols)
  colnames(tempDataFrame) <- safe_cols

  cols <- colnames(tempDataFrame)
  g_end <- length(cols)
  g <- 0

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = g_start:g_end)
  else
    progress_bar <- ""

  to_export <- c("cols", "family_test", "covariates", "independent_variable", "tempDataFrame",
    "independent_variable1stLevel", "independent_variable2ndLevel",
    "key", "transformation_y","assoc_exact_pvalue","g_end",
    "io_data_preparation","apply_stat_model_sig.formula","assoc_quantreg_permutation_model",
    "assoc_apply_stat_model_sig_formula", "data_distribution_info", "assoc_glm_model", "assoc_test_model", "assoc_test_model_paired", "Breusch_Pagan_pvalue",
    "progress_bar","progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace","signal_values","ssEnv","g_start",
    "assoc_execute_model", "assoc_is_family_dicotomic", "core_log_event","mediate","mediation","core_get_session_info", "samples_sql_condition",
    # AI-106 (2026-06-09): safe_to_real mapping must reach each foreach worker
    "safe_to_real")

  result_columns <- c("MARKER", "FIGURE", "AREA", "SUBAREA", "AREA_OF_TEST", "CI.LOWER", "CI.UPPER", "PVALUE", "STATISTIC_PARAMETER", "AIC_VALUE", "RESIDUALS", "SHAPIRO_PVALUE", "R_MODEL", "STD.ERROR", "N_PERMUTATIONS", "N_PERMUTATIONS_TEST")
  core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),  " Starting foreach with: ", g_end - g_start, " items")

  core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " I'll perform:",g_end - g_start," tests." )
  result_temp <- data.frame()
  # .packages loads SEMseeker in each worker so SEMseeker::: lookups resolve.
  # Internal helpers are prefixed with SEMseeker::: because they live in the
  # namespace (not in the caller's frame) and .export does not cover them.
  result_temp <- foreach::foreach(
    g = g_start:g_end,
    .combine = plyr::rbind.fill,
    .export = to_export,
    .packages = "SEMseeker"
  ) %dorng%
  # for(g in g_start:g_end)
  tryCatch({
    # NOTE: this tryCatch is intentional. doFuture internally wraps the foreach
    # body in tryCatch(error = identity), which returns the error *condition*
    # object as the task result. A simpleError has length 2 (message + call),
    # causing doFuture to throw "parsing result not of length one, but 2".
    # By catching errors ourselves and returning NULL, we prevent the condition
    # object from reaching doFuture's result-combination logic.
    # plyr::rbind.fill silently ignores NULL results.
    # AI-041: in-memory only; saveRDS happens at end-of-batch in the caller,
    # not per-gene (was the hot-path culprit causing ~5-7x slowdown).
    SEMseeker:::core_update_session_info(ssEnv, save_to_disk = FALSE)
    ssEnv <- SEMseeker:::core_get_session_info()

    burdenValue <- cols[g]
    if(ssEnv$showprogress)
      progress_bar(sprintf("doing genomic area: %s", stringr::str_pad(burdenValue, 10, pad = " ")))

    if(!is.null(tempDataFrame[,burdenValue]) & length(unique(tempDataFrame[,burdenValue]))>=2){


      #
      sig.formula <- SEMseeker:::assoc_apply_stat_model_sig_formula(family_test, burdenValue, independent_variable, covariates)
      model_result <- SEMseeker:::assoc_execute_model(family_test, tempDataFrame, sig.formula, burdenValue, independent_variable, transformation_y, (g_end - g_start < 10), samples_sql_condition, key)

      #
      local_result <- data.frame("INDIPENDENT_VARIABLE" = independent_variable)
      local_result$MARKER <- as.character(key$MARKER)
      local_result$FIGURE <-  as.character(key$FIGURE)
      local_result$AREA <-  as.character(key$AREA)
      local_result$SUBAREA <-  as.character(key$SUBAREA)
      # AI-248: which operator reduced the positions to this number. Absent for
      # the depths that do not aggregate, present wherever the caller declared
      # it - and part of the row identity, so two aggregations of the same
      # scope never collapse into one.
      if (!is.null(key$AGGREGATION))
        local_result$AGGREGATION <- as.character(key$AGGREGATION)
      # AI-255: SCOPE completes the identity of the row. Without it a burden over
      # the whole sample and a burden per gene are two rows that differ only in
      # AREA, which reads as two region classes rather than as two extents.
      if (!is.null(key$SCOPE))
        local_result$SCOPE <- as.character(key$SCOPE)
      # AI-106 (2026-06-09): reverse-map back to the upstream raw name
      # (HLA-A, chr10:100028204-100028508, ...) so the CSV preserves it
      # for enrichment / resume match. Fallback to burdenValue itself if
      # the mapping is missing (defensive - should not happen).
      local_result$AREA_OF_TEST <- if (burdenValue %in% names(safe_to_real)) {
        safe_to_real[[burdenValue]]
      } else {
        burdenValue
      }
      local_result$FAMILY_TEST <- family_test
      local_result$transformation_y <- transformation_y
      local_result$COVARIATES <- ifelse(length(covariates)>0,paste0(covariates,collapse=" "),NA)
      # local_result$bartlett.pvalue <- data_distribution_info(family_test, tempDataFrame, burdenValue, independent_variable)

      if (SEMseeker:::assoc_is_family_dicotomic(family_test))
      {
        #
        selector <- tempDataFrame[, independent_variable]==independent_variable1stLevel
        independent_variableData1stLevel <- stats::na.omit(tempDataFrame[selector,burdenValue])
        selector <- tempDataFrame[, independent_variable]==independent_variable2ndLevel
        independent_variableData2ndLevel <- stats::na.omit(tempDataFrame[selector,burdenValue])

        if(length(stats::na.omit(independent_variableData2ndLevel))==0 | length(stats::na.omit(independent_variableData1stLevel))==0)
        {
          SEMseeker:::core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " I skip this test because one of the two groups is empty." )
          colnames(local_result) <- toupper(colnames(local_result))
          local_result$PVALUE <- NA
        }
      }

      if (!SEMseeker:::assoc_is_family_dicotomic(family_test))
      {
        dependentVariableData <- as.numeric(stats::na.omit(tempDataFrame[!is.na(tempDataFrame[,independent_variable]),burdenValue]))
        independent_variableData <- as.numeric(stats::na.omit(tempDataFrame[  ,independent_variable]))

        if(sum(is.na(dependentVariableData)>0) | sum(is.na(independent_variableData)))
        {
          SEMseeker:::core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), "The submitted data are not factorial or numeric.")
          stop()
        }
      }

      if(nrow(model_result)>0)
        local_result <- cbind(local_result, model_result)

      colnames(local_result) <- toupper(colnames(local_result))

      # local_result
      # result_temp <- plyr::rbind.fill(result_temp, local_result)
      local_result
    }
  }, error = function(e) {
    SEMseeker:::core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " Skipping area '", if(exists("burdenValue")) burdenValue else "?",
              "': ", conditionMessage(e))
    NULL
  })


  # Report the actual count of non-NULL fits, not the chunk upper bound.
  # Before this change the log said "I performed: 1512 tests" even on chunks
  # where every gene was filtered out via area_to_remove or fitted to NULL,
  # which was confusing during resume runs that legitimately had 0 work left.
  n_done <- if (is.null(result_temp)) 0L else nrow(result_temp)
  if (n_done > 0L) {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " I performed:", n_done, " tests." )
  } else {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " No new tests to perform in this chunk (already cached)." )
  }

  # AI-041: end-of-foreach disk snapshot (workers used save_to_disk=FALSE
  # inside the per-gene loop; here we persist the session exactly once after
  # the parallel section closes).
  core_update_session_info(ssEnv, save_to_disk = TRUE)

  # & !is.null(result_temp)
  if(exists("result_temp") & !is.null(result_temp))
  {

    # result_temp <- unique(result_temp)
    result_temp <- result_temp %>% dplyr::distinct()

    # AI-257: no adjustment here any more. What this function holds is one
    # chunk - assoc_run_marker() splits a pivot at ceiling(6e6 / ncol)
    # rows - so any family it could form is a memory parameter, not a
    # statistical choice. assoc_analysis_save_results() is the only place where
    # every row of a family is together, and it computes all three levels there.
    #
    # The block removed here also split on `grepl("TOTAL", AREA_OF_TEST)`, and
    # TOTAL went with AI-255: it labelled rows built by composing aggregates
    # over a partition that is not disjoint. The predicate had been false on
    # every row since, so the two branches had quietly become one.
    colnames(result_temp) <- core_name_cleaning(colnames(result_temp))
    return(result_temp)
  }
  return(NULL)
}
