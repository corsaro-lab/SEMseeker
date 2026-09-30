#' @title Intra-study association (focused replication across models)
#' @description Re-run association analysis focusing on genomic regions that are
#'   statistically significant in a reference study, then merge results across
#'   two datasets for comparative inference.
#' @param inference_details_origin data.frame. Inference parameters for the
#'   reference (origin) study whose significant regions define the search space.
#' @param inference_details data.frame. Inference parameters for the target
#'   study on which the focused association analysis is performed.
#' @param result_folder character. Path to the SEMseeker result folder.
#' @param maxResources numeric. Maximum percentage of CPU cores to use (default 90).
#' @param parallel_strategy character. Parallelisation backend passed to
#'   \code{future}; one of \code{"multisession"}, \code{"sequential"},
#'   \code{"cluster"} (default \code{"multisession"}).
#'   Asking for \code{"multicore"} is accepted and converted to
#'   \code{"multisession"}: it means fork(), which is unsafe with this
#'   package's native thread pool on every platform that offers it, and
#'   absent on Windows. The conversion is logged.
#' @param start_fresh logical. If \code{TRUE}, delete previous results before
#'   running (default \code{FALSE}).
#' @param ... Additional named arguments passed to \code{core_init_env()}.
#' @return Invisibly \code{NULL}. Results are written to the inference
#'   sub-folder of \code{result_folder}.
#' @examples
#' result_dir <- tempdir()
#' \dontrun{
#' meta_association_replication(
#'   inference_details_origin = inference_study1,
#'   inference_details        = inference_study2,
#'   result_folder            = "~/semseeker_comparison/"
#' )
#' }
#' @export
meta_association_replication <- function(inference_details_origin, inference_details,result_folder,
  maxResources = 90, parallel_strategy  = "multisession",start_fresh = FALSE, ...)
{

  ssEnv <- core_init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)
  localKeys <- ssEnv$keys_markers_figures
  markers <- unique(localKeys$MARKER)
  for (a in seq_along(markers))
  {

    marker <- markers[a]

    inference_filename <- io_inference_file_name(inference_detail = inference_details, folder = ssEnv$result_folderInference, marker = marker, file_extension = "csv")

    if(file.exists(inference_filename))
      next

    inference_source <- assoc_results_get(inference_details_origin, marker,
      pvalue_column="PVALUE_ADJ_ALL_FDR", adjustment_method = "BH",
      area = "GENE", scope = "INSTANCE",
      omit_na = TRUE, significance = TRUE)

    # C-06: stop if the origin results carry a GENOME_BUILD that differs from the
    # current session.  Missing column = legacy run without provenance (warn only).
    if (nrow(inference_source) > 0L && "GENOME_BUILD" %in% colnames(inference_source)) {
      origin_builds <- unique(as.character(inference_source$GENOME_BUILD))
      origin_builds <- origin_builds[!is.na(origin_builds) & nzchar(origin_builds)]
      current_build <- if (!is.null(ssEnv$genome_build) && nzchar(ssEnv$genome_build))
        ssEnv$genome_build else "hg19"
      if (length(origin_builds) > 0L && !all(origin_builds == current_build))
        stop(
          "[cross_study_association] genome_build mismatch: origin results use '",
          paste(origin_builds, collapse = "/"), "' but current session uses '",
          current_build, "'.\n",
          "  Run liftover (C-07) before combining sessions from different assemblies."
        )
    } else if (nrow(inference_source) > 0L) {
      core_log_event("WARNING: [cross_study_association] Origin inference results",
                " have no GENOME_BUILD column - legacy run without provenance.",
                " Proceeding but cross-build safety cannot be guaranteed.")
    }

    areas <- as.vector(unique(inference_source$AREA_OF_TEST))

    association_analysis( inference_details = inference_details, result_folder = result_folder, areas_selection=areas,
      maxResources = maxResources,  parallel_strategy  = parallel_strategy,start_fresh = start_fresh,
      areas = unique(ssEnv$keys_areas_subareas_markers_figures$AREA), subareas = unique(ssEnv$keys_areas_subareas_markers_figures$SUBAREA),
      markers = marker, verbosity= ssEnv$verbosity, figures=unique(ssEnv$keys_areas_subareas_markers_figures$FIGURE),
      showprogress = ssEnv$showprogress)



    inference_filename_origin <- io_inference_file_name(inference_detail = inference_details_origin, folder = ssEnv$result_folderInference, marker = marker, file_extension = "csv")
    data_origin <- utils::read.csv2(inference_filename_origin)
    data <- utils::read.csv2(inference_filename)

    # remove from data origin where figure, area, subare and are_of_test are not in data
    data_origin$KEY <- paste(data_origin$FIGURE, data_origin$AREA, data_origin$SUBAREA, data_origin$AREA_OF_TEST, sep = "_")
    data$KEY <- paste(data$FIGURE, data$AREA, data$SUBAREA, data$AREA_OF_TEST, sep = "_")
    data_origin <- data_origin[!(data_origin$KEY %in% data$KEY),]
    data_origin$FAMILY_TEST <- inference_details$family_test

    data <- rbind(data, data_origin)
    data$KEY <- NA
    assoc_analysis_save_results(data, inference_filename, inference_details$family_test, FALSE)

  }


}
