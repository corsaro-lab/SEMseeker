#' Association analysis of SEMseeker results
#'
#' Run statistical association models between SEM metrics and a phenotype
#' variable. Supports group tests (Wilcoxon, t-test), GLM families (gaussian,
#' poisson, binomial), quantile regression, correlations (Pearson, Kendall,
#' Spearman), and multi-covariate formulas (e.g.
#' \code{MUTATIONS_* ~ covariate1 + covariate2}).
#'
#' @param inference_details data.frame. Each row defines one analysis run.
#'   Required columns:
#'   \describe{
#'     \item{independent_variable}{Sample sheet column used as grouping /
#'       covariate variable.}
#'     \item{family_test}{Statistical model: \code{"wilcoxon"},
#'       \code{"stats::t.test"}, \code{"gaussian"}, \code{"poisson"},
#'       \code{"binomial"}, \code{"pearson"}, \code{"kendall"},
#'       \code{"spearman"}, or quantile regression as
#'       \code{"quantreg_<tau>_<runs>"} (e.g. \code{"quantreg_0.25_2000"}).}
#'     \item{transformation_y}{Transformation applied to the dependent variable:
#'       \code{"none"}, \code{"scale"}, \code{"log"}, \code{"log2"},
#'       \code{"log10"}, \code{"exp"}, or
#'       \code{"quantile_<n>"} (e.g. \code{"quantile_3"}).}
#'     \item{scope}{Required. Which of the two aggregations to run, over the
#'       region classes of the call (\code{areas}, \code{subareas}).
#'       \code{"SAMPLE"} reduces the positions of each class to one number per
#'       sample: the burden, or its density, or a descriptor of the signal.
#'       \code{"INSTANCE"} reduces them to one number per instance of the class:
#'       one row per gene, per island, per cytoband, per probe.
#'
#'       The two are \strong{mutually exclusive}. They used to be produced
#'       together and written into the same file, which made every result the
#'       union of two different questions; a request that wants both now writes
#'       two rows, and each carries its own \code{SCOPE} in the output.
#'
#'       It has no default on purpose. Any default would answer half of what was
#'       asked and say nothing about the other half.}
#'     \item{aggregation}{Required. How the positions are reduced to the one
#'       number the model is fitted on: \code{"SUM"}, \code{"MEAN"},
#'       \code{"MEDIAN"}, \code{"VARIANCE"}, \code{"IQR"}, \code{"MODELOW"} or
#'       \code{"MODEHIGH"}. While every marker admitted exactly one operator
#'       this could stay implicit; a scope now carries several, so the request
#'       has to name the one it wants. Which are admissible depends on the
#'       marker: a count carries \code{SUM} (the burden) and \code{MEAN} (the
#'       density, the form comparable across regions of different size), while
#'       its median and IQR are degenerate; the two modes exist only for the
#'       signal on the beta scale. A request no marker of the run admits is
#'       dropped with a warning naming it, not answered with an empty result.}
#'   }
#'
#'   The markers are \strong{not} named in \code{inference_details} either.
#'   They are chosen by the \code{markers} argument of this call and every row
#'   of the request is tested on all of them, one result file per marker. A
#'   \code{marker} column was documented here for a long time and never existed
#'   in the vocabulary, so a request written from that description was rejected
#'   as carrying an unknown column.
#'
#'   The region classes are \strong{not} named in \code{inference_details}: they
#'   are the \code{(AREA, SUBAREA)} pairs of the run, declared with the
#'   \code{areas} and \code{subareas} arguments of this call and built at
#'   runtime. Both scopes range over the same pairs, so
#'   \code{scope = "SAMPLE"} with \code{areas = c("GENE", "PROBE")} gives the
#'   burden over the gene probes and the burden over the whole sample, and
#'   \code{scope = "INSTANCE"} gives one row per gene and one row per probe.
#'   The artefact is built on the way in if it does not exist yet, so a class no
#'   previous run foresaw costs one scan of the position pivot rather than a
#'   rerun.
#'
#'   Result rows carry the coordinates as columns: \code{MARKER},
#'   \code{FIGURE}, \code{SCOPE}, \code{AREA}, \code{SUBAREA},
#'   \code{AGGREGATION}, never a class name squashed into one of them.
#' @param result_folder character. Path to the SEMseeker result folder.
#' @param maxResources numeric. Maximum percentage of CPU cores to use
#'   (default 90).
#' @param parallel_strategy character. Parallelisation backend; possible
#'   values: \code{"none"}, \code{"multisession"}, \code{"sequential"},
#'   \code{"multicore"}, \code{"cluster"} (default \code{"multicore"}).
#' @param start_fresh logical. If \code{TRUE}, delete previous inference
#'   results before running (default \code{FALSE}).
#' @param ... Additional arguments passed to \code{core_init_env()}.
#'
#' @return Invisibly \code{NULL}. Inference result CSV files are written to
#'   the \code{Inference/} sub-folder of \code{result_folder}, one file per
#'   marker/area/family combination defined in \code{inference_details}.
#' @importFrom doRNG %dorng%
#' @examples
#' result_dir <- tempdir()
#' \dontrun{
#' association_analysis(
#'   inference_details = data.frame(
#'     independent_variable = "Sample_Group",
#'     family_test          = "wilcoxon",
#'     transformation_y     = "none",
#'     aggregation          = "MEAN",
#'     scope                = "INSTANCE"
#'   ),
#'   result_folder     = "~/semseeker_results/",
#'   markers           = "DELTARP",
#'   areas             = "GENE",
#'   subareas          = "WHOLE",
#'   multiple_test_adj = "BH"
#' )
#' }
#' @export
association_analysis <- function(inference_details, result_folder, maxResources = 90,
  parallel_strategy = "multicore", start_fresh = FALSE, ...) {

  arguments <- list(...)
  areas_selection <- c()
  if (!is.null(arguments[["areas_selection"]])) {
    areas_selection <- arguments$areas_selection
    arguments[["areas_selection"]] <- NULL
  }

  ssEnv <- core_init_env(result_folder = result_folder, maxResources = maxResources,
    parallel_strategy = parallel_strategy, start_fresh = FALSE, ...)

  core_log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"),
    " SEMseeker will perform the association analysys for project \n in ",
    ssEnv$result_folderData)

  if (start_fresh) unlink(ssEnv$result_folderInference, recursive = TRUE)
  io_dir_check_and_create(ssEnv$result_folderInference, c())

  localKeys <- ssEnv$keys_markers_figures

  sem_deltaX_get()
  anno_annotate_position_pivots()

  inference_details <- assoc_validate_inference_schema(unique(inference_details))
  # AI-248: shape first, then meaning. Refuse a request that cannot be honoured
  # before any result is written: checking it inside the per-marker loop would
  # surface the mistake after part of the output exists.
  #
  # AI-308: the scope goes first, because which aggregations are admissible
  # depends on it. The two peaks of a bimodal density need one big group, so
  # they exist at SCOPE = SAMPLE and nowhere else; asking for them per instance
  # used to travel all the way to io_pivot_build() and stop there, with the run
  # already under way.
  inference_details <- assoc_validate_scope(inference_details)
  inference_details <- assoc_validate_aggregation(inference_details)
  # AI-309: and the model. A family test that is absent or unknown used to make
  # the row vanish from the loop below, leaving a result file indistinguishable
  # from one where the test had run.
  inference_details <- assoc_validate_family(inference_details)

  for (z in seq_len(nrow(inference_details))) {
    start_time <- Sys.time()
    inference_detail <- inference_details[z, ]
    filter_p_value <- if (!is.null(inference_detail$filter_p_value))
      inference_detail$filter_p_value else TRUE

    core_log_inference_header(inference_detail)

    # AI-309: validated at the door, so there is nothing to check and nothing to
    # skip here. The `next` this replaces is the reason a malformed request
    # could produce a complete-looking file.
    family_test <- util_split_and_clean(inference_detail$family_test)

    # AI-255: the models read artefacts, not columns — assoc_run_marker() opens
    # the pivot for every key, collapsed or not. So what this needs from the
    # sample sheet is the phenotype and the covariates, and joining the
    # per-sample statistics onto it would build artefacts nobody then reads:
    # io_feature_colname() has exactly one caller left, the composer inside
    # sem_study_summary_get(), and nothing reads those names back.
    #
    # The join is still done when the request names a feature the plain sheet
    # does not have — adjusting for the global burden is a legitimate thing to
    # ask — but it is no longer paid for on every run by default.
    study_summary <- sem_study_summary_get(inference_detail$samples_sql_condition,
                                           with_sample_stats = FALSE)
    wanted_cols <- c(gsub(" ", "", as.character(inference_detail$independent_variable)),
                     util_split_and_clean(inference_detail$covariates))
    wanted_cols <- wanted_cols[nzchar(wanted_cols) & !is.na(wanted_cols)]
    if (!is.null(study_summary) && !all(wanted_cols %in% colnames(study_summary))) {
      core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " The request names ", paste(setdiff(wanted_cols, colnames(study_summary)),
                                             collapse = ", "),
                ", which the sample sheet does not carry: joining the per-sample ",
                "features as well.")
      # AI-308: the request no longer names region classes: they are the
      # (AREA, SUBAREA) pairs of the run. A covariate the sheet does not carry
      # can name any of them, plus "SAMPLE" for the unrestricted feature, so the
      # join offers the whole registry rather than a list the request no longer
      # has.
      study_summary <- sem_study_summary_get(
        inference_detail$samples_sql_condition,
        regions = unique(c("SAMPLE", as.character(ssEnv$keys_areas_subareas$COMBINED))))
    }
    prep <- sem_prepare_study_for_analysis(inference_detail, study_summary, family_test)
    if (is.null(prep)) next

    processed_items <- 0L
    last_results <- data.frame()
    last_filename <- NULL

    for (marker in unique(localKeys$MARKER)) {
      keys <- unique(localKeys[localKeys$MARKER == marker, ])
      fileNameResults <- io_inference_file_name(prep$inference_detail, marker,
        ssEnv$result_folderInference,
        prefix = ifelse(length(areas_selection) == 0, "",
          paste(areas_selection, "_", sep = "")))
      core_log_event("JOURNAL:", "Result saved into file:", fileNameResults, ".")

      # AI-255: one road. There used to be two calls here, chosen by
      # depth_analysis, because the collapsed artefact and the per-instance one
      # had different shapes — a table of columns against a pivot of rows. They
      # have the same shape now, so a model handed a row does not know, and has
      # no reason to ask, whether the key of that row is a gene symbol or
      # PROBE_WHOLE. It fits. The scope travels in the key; the batch-family
      # exclusion travels with it (see .assoc_marker_keys).
      dn <- assoc_run_marker(prep, marker, family_test, fileNameResults,
        filter_p_value, ssEnv, selected_areas = areas_selection,
        data.frame(), start_time, processed_items, ...)
      results <- dn$results
      processed_items <- dn$processed_items

      last_results  <- results
      last_filename <- fileNameResults

      # AI-061+ (2026-06-09): volcano plot for this marker right after the
      # CSV is finalised. One call per marker; assoc_volcano_plot_inference
      # splits internally by (AREA, SUBAREA) and writes one PNG per
      # combination under <result_folder>/Chart/VOLCANO/. Best-effort:
      # plot failure must not abort the analysis loop — log WARNING and
      # continue with the next marker.
      tryCatch(
        assoc_volcano_plot_inference(
          inference_detail = prep$inference_detail,
          result_folder    = ssEnv$result_folder,
          markers          = marker
        ),
        error = function(e) {
          core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                    " assoc_volcano_plot_inference failed for marker '", marker,
                    "': ", conditionMessage(e))
        }
      )
    }

    util_finalize_job_results(last_results, prep$inference_detail, family_test,
      filter_p_value, last_filename, start_time, processed_items)
  }

  core_close_env()
}
