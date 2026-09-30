#' Open the session of a meta analysis over a set of studies
#'
#' One door for the whole \code{meta_} family. Every function in it does the
#' same four things before it can start, and doing them in four places is how
#' they came to do them differently.
#'
#' The order matters. The studies are turned into one table and checked for
#' being studies at all; then their provenance is compared, because pooling
#' estimates keyed by position across two different genome builds pools
#' positions that do not mean the same thing; and only then is the meta
#' analysis's own session opened on \code{result_folder}, which is where it
#' writes and not where it reads.
#'
#' The studies table is kept in the session and written to disk with it, so the
#' functions downstream read it from there instead of each deriving it again
#' from whatever they were handed, and so the provenance survives the run.
#'
#' That last part is the point of the layer. A meta analysis writes its own
#' session, with the same shape a study's session has, which makes it an outer
#' layer rather than a report: its result folder satisfies the same check
#' [meta_studies_normalise()] applies to a study, so a meta analysis can be the
#' input of another one. What keeps that chain readable is \code{studies} on
#' disk: without it the second meta analysis would know its inputs and not know
#' theirs.
#'
#' @param studies Study paths in any of the forms [meta_studies_normalise()]
#'   accepts.
#' @param result_folder Directory the meta analysis writes to.
#' @param ... Passed to the session setup.
#'
#' @return The session environment, carrying \code{studies}, a table of
#'   \code{STUDY} and \code{STUDY_FOLDER}.
#'
#' @keywords internal
#' @noRd
meta_session_open <- function(studies, result_folder, ...)
{
  studies <- meta_studies_normalise(studies)

  # Stops on differing genome builds, warns on differing technology.
  core_check_session_compatibility(studies$STUDY_FOLDER)

  ssEnv <- core_init_env(result_folder = result_folder, start_fresh = FALSE, ...)
  ssEnv$studies <- studies
  # Persisted, not just held: this is what makes the chain of provenance
  # readable from a result folder alone.
  core_update_session_info(ssEnv, save_to_disk = TRUE)
  ssEnv
}
