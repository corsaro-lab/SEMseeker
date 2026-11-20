# core_folder_delete - the only way SEMseeker removes a folder of the user.
#
# start_fresh = TRUE used to call unlink(result_folder, recursive = TRUE) on
# whatever path the user passed: a typo, the home directory or a folder holding
# other work was deleted without a question. A package may delete only what it
# created, so every removal of a result folder, or of a folder inside one, goes
# through this guard.
#
# A folder is SEMseeker's own when it holds the session file that
# core_init_env() writes, Log/session_info.rds. Deleting is allowed when the
# folder is empty or missing (nothing can be lost), or when it lies inside a
# result folder that carries that marker. Everything else is refused before
# anything is touched, because a refusal that stops the run is recoverable and
# a deletion is not.

#' Whether a folder is a SEMseeker result folder.
#'
#' @param result_folder path of the folder to look at
#' @return TRUE when the folder holds Log/session_info.rds, the file
#'   core_init_env() writes into every result folder it creates
#' @noRd
core_result_folder_is_own <- function(result_folder) {
  file.exists(file.path(result_folder, "Log", "session_info.rds"))
}

#' Whether one path is the same as another or lies inside it.
#'
#' @param outer_path an absolute, normalised path
#' @param inner_path an absolute, normalised path
#' @return TRUE when inner_path equals outer_path or is below it
#' @noRd
core_path_contains <- function(outer_path, inner_path) {
  outer_prefix <- if (endsWith(outer_path, "/")) outer_path else paste0(outer_path, "/")
  identical(outer_path, inner_path) || startsWith(inner_path, outer_prefix)
}

#' Delete a folder of the user, but only one SEMseeker is entitled to delete.
#'
#' Refuses, without touching anything, a folder that is or contains the home
#' directory, the working directory or the root; and a non-empty folder that
#' does not lie inside a result folder carrying SEMseeker's session file.
#'
#' @param folder the folder to delete
#' @param result_folder the result folder that must contain `folder` and carry
#'   the session file; by default `folder` itself, which is the case of
#'   start_fresh on a whole result folder
#' @return TRUE, invisibly, when a folder was deleted; FALSE when there was
#'   nothing to delete
#' @noRd
core_folder_delete <- function(folder, result_folder = folder) {
  if (!dir.exists(folder))
    return(invisible(FALSE))

  folder_path <- normalizePath(folder, winslash = "/", mustWork = TRUE)
  protected_paths <- unique(normalizePath(c(path.expand("~"), getwd(), "/"), winslash = "/", mustWork = FALSE))
  for (protected_path in protected_paths) {
    if (core_path_contains(folder_path, protected_path))
      stop("SEMseeker: refusing to delete '", folder_path, "': it is or contains '", protected_path,
           "'. Choose a result folder of its own.", call. = FALSE)
  }

  if (length(list.files(folder_path, all.files = TRUE, no.. = TRUE)) > 0L) {
    result_path <- normalizePath(result_folder, winslash = "/", mustWork = FALSE)
    if (!core_path_contains(result_path, folder_path) || !core_result_folder_is_own(result_path))
      stop("SEMseeker: refusing to delete '", folder_path, "': it is not empty and it is not a SEMseeker result folder ",
           "(no Log/session_info.rds in '", result_path, "'). Empty it yourself or choose another result_folder.", call. = FALSE)
  }

  unlink(folder_path, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}
