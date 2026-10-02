#' @title Probe coverage analysis report
#' @description Generate a coverage analysis report for Illumina methylation
#'   array data, summarising probe representation across genomic regions.
#'   WGBS data are not supported and will cause the function to stop with an
#'   informative error.
#' @param signal_data character. Path to the signal parquet file under
#'   \code{Data/Pivots/SIGNAL/}, or a data.frame already loaded into memory.
#'   The file is named \code{SIGNAL_<FIGURE>_PROBE_WHOLE_<GENOME_BUILD>.parquet},
#'   where \code{FIGURE} is the scale the run was written on - \code{BETA} for a
#'   proportion bounded in [0,1], \code{MVALUE} for the logit-transformed one.
#' @param result_folder character. Path to the SEMseeker result folder.
#' @param maxResources numeric. Maximum percentage of CPU cores to use
#'   (default 90).
#' @param parallel_strategy character. Parallelisation backend passed to
#'   \code{future} (default \code{"multisession"}).
#'   Asking for \code{"multicore"} is accepted and converted to
#'   \code{"multisession"}: it means fork(), which is unsafe with this
#'   package's native thread pool on every platform that offers it, and
#'   absent on Windows. The conversion is logged.
#' @param ... Additional named arguments passed to \code{core_init_env()}.
#' @return Invisibly \code{NULL}. Coverage tables and charts are written to
#'   the result folder.
#' @examples
#' result_dir <- tempdir()
#' \dontrun{
#' # BETA is the scale this particular run was written on; a run on the
#' # logit scale produces SIGNAL_MVALUE_PROBE_WHOLE_HG19.parquet instead.
#' sem_coverage_analysis_report(
#'   signal_data   = "~/semseeker_results/Data/Pivots/SIGNAL/SIGNAL_BETA_PROBE_WHOLE_HG19.parquet",
#'   result_folder = "~/semseeker_results/"
#' )
#' }
#' @export
sem_coverage_analysis_report <- function (signal_data, result_folder, maxResources = 90, parallel_strategy  = "multisession", ...)
{

  ssEnv <- core_init_env( result_folder =  result_folder, maxResources =  maxResources, ...)
  # signal_data_path <- io_pivot_file_name_parquet("SIGNAL", io_signal_figure(),"PROBE","WHOLE")
  # if (!file.exists(signal_data_path))
  # {
  #   core_log_event("ERROR:  ", format(Sys.time(), "%a %b %d %X %Y"), " Signal data is missing")
  #   core_close_env()
  #   stop()
  # }
  # signal_data <- polars::pl$read_parquet(signal_data_path)$to_data_frame()
  signal_data <- io_source_data_get(signal_data)
  ssEnv <- core_get_meth_tech(signal_data)
  if (ssEnv$tech=="WGBS")
  {
    core_log_event("ERROR:  ", format(Sys.time(), "%a %b %d %X %Y"), " WGBS data is not supported for coverage analysis")
    core_close_env()
    stop()
  }
  sem_coverage_analysis(signal_data$AREA)
  core_close_env()
}
