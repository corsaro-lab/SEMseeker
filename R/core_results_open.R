#' Open the session of an existing result folder for a chart or a reader.
#'
#' The charts and readers take the session as they find it, so they used to
#' work only in the R session that had just run semseeker(): in a new one they
#' stopped with "called without result folder!", and results could not be
#' reopened the next day. Naming the folder opens its session from the file
#' the run wrote. core_init_env() rebuilds every sub-folder path from the
#' folder it is given, so a folder that was moved or copied works where it now
#' is, and an option set by the run is kept because an option is only set
#' while it is still empty.
#'
#' A folder without Log/session_info.rds is refused: opening it would create
#' an empty session there, and a mistyped path would end in an empty chart
#' rather than an error.
#'
#' @param result_folder character path of a folder written by semseeker(), or
#'   NULL to use the session already open. A folder whose session is already
#'   open is not opened again.
#' @return TRUE when this call opened the session, so the caller closes it on
#'   exit; FALSE when it used the one already open.
#' @noRd
core_results_open <- function(result_folder)
{
  if (is.null(result_folder))
    return(FALSE)

  open_session <- tryCatch(core_get_session_info(), error = function(e) NULL)
  if (length(open_session) > 0 && !is.null(open_session$result_folder) &&
      .core_same_folder(open_session$result_folder, result_folder))
    return(FALSE)

  if (!core_result_folder_is_own(result_folder))
    stop("'", result_folder, "' is not a SEMseeker result folder: it holds no ",
         "Log/session_info.rds. Name the folder given to semseeker().",
         call. = FALSE)

  core_init_env(result_folder = result_folder)
  TRUE
}
