#' @title Enrichment analysis (pathway + phenotype)
#' @description Dispatcher per l'analisi di arricchimento di pathway e
#'   phenotype sui geni con epimutazioni significative identificate da
#'   `association_analysis`. Supporta multipli backend opt-in:
#'   WebGestalt, STRINGdb, pathfindR, CTDbase (ctdR), e Phenolyzer.
#'
#'   \strong{Phenolyzer requirement:} the \code{phenolyzer} and
#'   \code{Phenolyzer_STRINGdb} / \code{Phenolyzer_WebGestalt} backends require
#'   Phenolyzer to be installed and configured at the operating-system level
#'   (see \url{https://phenolyzer.wglab.org}). Supply the path to the
#'   Phenolyzer binary directory via \code{phenolyzer_folder_bin}.
#'
#' @param inference_details data.frame. Inference parameter table (must contain
#'   the per-gene artefacts (SCOPE = INSTANCE, AREA = GENE) are
#'   processed).
#' @param adjust_per_area_s logical vector. Whether to adjust p-values per area
#'   for each \code{pvalue_columns} entry.
#' @param adjust_globally_s logical vector. Whether to adjust p-values globally
#'   for each \code{pvalue_columns} entry.
#' @param pvalue_columns character vector. Column name(s) in the inference
#'   results to use as p-value filter.
#' @param adjustment_methods character vector. Multiple-testing correction
#'   method(s) applied (e.g. \code{"BH"}).
#' @param alphas numeric vector. Significance thresholds to iterate over.
#' @param study character. Study identifier used for labelling outputs.
#' @param significance logical. Whether to filter by significance.
#' @param statistic_parameter character. Statistic column used for ranking.
#' @param path_dbs character. Path database(s) for pathfindR.
#' @param phenolyzer_folder_bin character. Path to the Phenolyzer binary
#'   directory (required when \code{phenolyzer = TRUE} or related flags).
#' @param disease character. Disease keyword passed to Phenolyzer.
#' @param phenolyzer logical. Enable Phenolyzer backend (default \code{FALSE}).
#'   Requires Phenolyzer installed system-wide.
#' @param WebGestalt logical. Enable WebGestalt backend (default \code{FALSE}).
#' @param pathfindr logical. Enable pathfindR backend (default \code{FALSE}).
#' @param STRINGdb logical. Enable STRINGdb backend (default \code{FALSE}).
#' @param Phenolyzer_STRINGdb logical. Enable combined Phenolyzer+STRINGdb
#'   backend (default \code{FALSE}).
#' @param Phenolyzer_WebGestalt logical. Enable combined Phenolyzer+WebGestalt
#'   backend (default \code{FALSE}).
#' @param ctdR logical. Enable CTDbase backend via the \pkg{ctdR} package
#'   (default \code{FALSE}).
#' @param result_folder character. Path to the SEMseeker result folder.
#' @param maxResources numeric. Maximum percentage of CPU cores to use
#'   (default 90).
#' @param parallel_strategy character. Parallelisation backend (default
#'   \code{"multicore"}).
#' @param ... Additional named arguments passed to \code{core_init_env()}.
#' @return Invisibly \code{NULL}. Pathway enrichment results are written to the
#'   pathway sub-folder of \code{result_folder}.
#' @examples
#' result_dir <- tempdir()
#' \dontrun{
#' enrichment_analysis(
#'   inference_details  = inference_df,
#'   adjust_per_area_s  = TRUE,
#'   adjust_globally_s  = FALSE,
#'   pvalue_columns     = "PVALUE_ADJ_ALL_BH",
#'   adjustment_methods = "BH",
#'   alphas             = 0.05,
#'   result_folder      = "~/semseeker_results/",
#'   WebGestalt         = TRUE
#' )
#' }
#' @export
enrichment_analysis <- function(inference_details, adjust_per_area_s, adjust_globally_s, pvalue_columns, adjustment_methods,alphas,
  study, significance,statistic_parameter, path_dbs, phenolyzer_folder_bin,disease,
  phenolyzer=FALSE, WebGestalt=FALSE, pathfindr=FALSE,STRINGdb=FALSE,Phenolyzer_STRINGdb=FALSE,Phenolyzer_WebGestalt=FALSE,ctdR=FALSE,
  result_folder, maxResources = 90, parallel_strategy  = "multicore", ...)
{
  start_fresh <- FALSE
  ssEnv <- core_init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy,
    start_fresh = start_fresh, ...)

  arguments <- list(...)

  # The 4 parameter arrays must have the same length: each index i defines
  # one (pvalue_column, adjustment_method, adjust_per_area, adjust_globally)
  # tuple to iterate over. Diagnostic shows the actual length of each so the
  # user knows exactly which line of their setup needs fixing.
  lens <- c(
    pvalue_columns     = length(pvalue_columns),
    adjustment_methods = length(adjustment_methods),
    adjust_per_area_s  = length(adjust_per_area_s),
    adjust_globally_s  = length(adjust_globally_s)
  )
  if (length(unique(lens)) != 1) {
    msg <- sprintf(
      paste0("Parameter array lengths must match across:\n  %s\n",
             "Each entry of pvalue_columns needs a matching entry in ",
             "adjustment_methods, adjust_per_area_s, adjust_globally_s. ",
             "Fix this in the inference setup file (typically `02_setup_<study>.R`)."),
      paste(sprintf("%s (length=%d)", names(lens), lens), collapse = "\n  ")
    )
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " ", msg)
    stop(msg, call. = FALSE)
  }


  # AI-311: the invariant is declared in one place and refused here, at the
  # door, instead of being applied as a subset() five call sites downstream. A
  # folder with no per-instance gene row used to yield no enrichment at all,
  # silently, which reads as "nothing was significant" when it means "the input
  # was never computed". Since the two aggregation branches stopped being
  # additive (AI-308) that folder is easy to produce by accident: a run at
  # scope = "SAMPLE" is complete, legitimate, and feeds no pathway analysis.
  enrich_input_assert(inference_details)

  for (alpha in alphas)
  {
    for (id in seq_len(nrow(inference_details)))
    {
      inference_detail <- inference_details[id,]
      inference_detail_prettified <- t(inference_detail)
      # Generate a plain text table using kable
      inference_detail_prettified <- knitr::kable(inference_detail_prettified, format = "simple",
        align = "l",    # Left align for all columns
        digits = 2,     # Number of digits for numeric columns
        row.names = TRUE) # Suppress row names

      inference_detail_prettified <- paste(inference_detail_prettified, collapse = "\n")
      core_log_event("JOURNAL: ##############################################################################################################")
      core_log_event("JOURNAL: ", format(Sys.time(), "%a %b %d %X %Y"), " \nStarting pathway for inference detail:\n", inference_detail_prettified)
      core_log_event("JOURNAL: Gene are selected by the sql_area_selection in the inference details above, also filtered using ", pvalue_columns, "columns. \nWith an alpha limit having value: ", alpha)

      for (i in seq_along(pvalue_columns)){

        pvalue_column <- pvalue_columns[i]
        ssEnv$alpha <- alpha
        core_update_session_info(ssEnv)

        adjustment_method <- adjustment_methods[i]
        adjust_per_area <- adjust_per_area_s[i]
        adjust_globally <- adjust_globally_s[i]

        if(ctdR)
          enrich_ctdR(
            study =  study,
            statistic_parameter=statistic_parameter,
            adjust_per_area = adjust_per_area,
            adjust_globally = adjust_globally,
            adjustment_method = adjustment_method,
            pvalue_column = pvalue_column,
            inference_detail = inference_detail,
            significance = significance
          )

        if(phenolyzer)
          enrich_phenotype_phenolyzer(
            study =  study,
            disease = disease,
            phenolyzer_folder_bin = phenolyzer_folder_bin,
            minimum_score = 0.5,
            statistic_parameter = statistic_parameter,
            adjust_per_area = adjust_per_area,
            adjust_globally = adjust_globally,
            adjustment_method = adjustment_method,
            pvalue_column = pvalue_column,
            inference_detail = inference_detail,
            significance = significance
          )

        if(Phenolyzer_STRINGdb)
          enrich_Phenolyzer_STRINGdb(
            study =  study,
            statistic_parameter=statistic_parameter,
            adjust_per_area = adjust_per_area,
            adjust_globally = adjust_globally,
            adjustment_method = adjustment_method,
            pvalue_column = pvalue_column,
            inference_detail = inference_detail,
            significance = significance,
            disease = disease
          )


        if (Phenolyzer_WebGestalt)
          enrich_Phenolyzer_WebGestalt(
            study = study,
            types=c("BP","MF"),
            enrich_methods = c("ORA"),
            adjust_per_area = adjust_per_area,
            adjust_globally = adjust_globally,
            adjustment_method = adjustment_method,
            pvalue_column = pvalue_column,
            inference_detail = inference_detail,
            significance = significance,
            disease = disease
          )


        if(STRINGdb)
          enrich_STRINGdb(
            study =  study,
            statistic_parameter=statistic_parameter,
            adjust_per_area = adjust_per_area,
            adjust_globally = adjust_globally,
            adjustment_method = adjustment_method,
            pvalue_column = pvalue_column,
            inference_detail = inference_detail,
            significance = significance
          )

        if(pathfindr)
          enrich_pathfindR(
            study =  study,
            path_dbs = path_dbs,
            iterations = 20,
            statistic_parameter=statistic_parameter,
            adjust_per_area = adjust_per_area,
            adjust_globally = adjust_globally,
            adjustment_method = adjustment_method,
            pvalue_column = pvalue_column,
            inference_detail = inference_detail,
            significance = significance
          )


        if (WebGestalt)
          enrich_WebGestalt(
            study = study,
            types=c("BP","MF"),
            enrich_methods = c("ORA"),
            adjust_per_area = adjust_per_area,
            adjust_globally = adjust_globally,
            adjustment_method = adjustment_method,
            pvalue_column = pvalue_column,
            inference_detail = inference_detail,
            significance = significance
          )

      }
    }
  }

  core_close_env()
}
