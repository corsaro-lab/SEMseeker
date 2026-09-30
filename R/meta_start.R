#' Open a meta analysis: its studies, its space, its selection and its files
#'
#' Every function in the \code{meta_} family starts the same way, and the six of
#' them started it six different ways, which is how they came to check different
#' things or nothing. This is that opening, once.
#'
#' The order is not incidental. Each step refuses what the next one could not
#' recover from:
#'
#' \enumerate{
#'   \item the studies become one table of label and folder, and a path that is
#'     not a study is refused here rather than read;
#'   \item their provenance is compared, because pooling results keyed by
#'     position across two genome builds pools positions that do not mean the
#'     same thing;
#'   \item the region classes each study analysed are read from its own session,
#'     with no selection applied, and intersected strictly: that is the space,
#'     and what strictness left out is reported;
#'   \item the caller's selection is checked against that space one coordinate
#'     at a time and then intersected with it;
#'   \item the artefacts that answer the request are named, per study, and their
#'     absence is refused before any of them is read.
#' }
#'
#' The meta analysis's own session is opened \strong{last}, and that is a
#' correctness requirement rather than tidiness: \code{core_init_env()} replaces
#' the current session, so reading the studies leaves the last study's session
#' current. A function that opened its own session first and read the studies
#' afterwards would write its results into the last study it read.
#'
#' @param studies Study paths, in any form [meta_studies_normalise()] accepts.
#' @param result_folder Where the meta analysis writes. This becomes a session
#'   of the same shape a study's is, which is what lets a meta analysis be the
#'   input of another one.
#' @param inference_detail One row of the inference specification: both the
#'   request and the name of the files that answer it.
#' @param artefact \code{"inference"} or \code{"enrichment"}.
#' @param markers,figures,areas,subareas Optional selection within the space.
#' @param enrichment_package Required for \code{artefact = "enrichment"}.
#' @param pvalue_column,alpha,significance Carried into the artefact names.
#' @param ... Passed to the session setup, for the studies and for the meta
#'   analysis alike.
#'
#' @return A list with \code{ssEnv} (the meta analysis's session, current on
#'   return), \code{studies}, \code{space}, \code{excluded}, \code{keys},
#'   \code{requested}, \code{dropped}, \code{files} and \code{artefact}.
#'
#' @keywords internal
#' @noRd
meta_start <- function(studies, result_folder, inference_detail,
                       artefact = c("inference", "enrichment"),
                       markers = NULL, figures = NULL, areas = NULL, subareas = NULL,
                       enrichment_package = NULL,
                       pvalue_column = "PVALUE_ADJ_ALL_BH",
                       alpha = NULL, significance = TRUE, ...)
{
  artefact <- match.arg(artefact)

  studies <- meta_studies_normalise(studies)
  core_check_session_compatibility(studies$STUDY_FOLDER)

  key_sets <- meta_key_sets_read(studies, ...)
  space <- meta_key_space(key_sets)
  selected <- meta_keys_select(space$space, markers = markers, figures = figures,
                               areas = areas, subareas = subareas)
  files <- meta_artefact_files(studies, selected$keys, inference_detail,
                               artefact = artefact,
                               enrichment_package = enrichment_package,
                               pvalue_column = pvalue_column,
                               alpha = alpha, significance = significance)

  # Last, so that it is the current session when this returns.
  ssEnv <- core_init_env(result_folder = result_folder, start_fresh = FALSE, ...)
  ssEnv$studies <- studies
  # Persisted, not merely held: it is what makes the chain of provenance
  # readable from a result folder alone, once this folder is itself an input.
  core_update_session_info(ssEnv, save_to_disk = TRUE)

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                 " meta analysis over ", nrow(studies), " studies (",
                 paste(studies$STUDY, collapse = ", "), "): ",
                 nrow(space$space), " region classes shared, ",
                 nrow(space$excluded), " declared by some study and not by all, ",
                 nrow(selected$keys), " selected of ", selected$requested,
                 " named, ", selected$dropped, " outside the space.")

  list(ssEnv = ssEnv, studies = studies,
       space = space$space, excluded = space$excluded,
       keys = selected$keys, requested = selected$requested,
       dropped = selected$dropped, files = files, artefact = artefact)
}
