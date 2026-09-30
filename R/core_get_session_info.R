core_get_session_info <- function(result_folder=NULL)
{

  # try from environment
  ssEnv <- .pkgglobalenv$ssEnv

  # When a folder is named, the folder identifies the session. An in-memory
  # session belonging to a DIFFERENT folder is not a recovery of this one: it is
  # the previous analysis still standing, and its options would be inherited by
  # an analysis that never asked for them. That is silent: no message, and the
  # second analysis runs with the first one's parameters.
  #
  # Resuming still works, because resuming means the same folder and the folders
  # then match. A session whose result_folder is not set yet is left alone: that
  # is the window inside core_init_env() before the paths are assigned.
  if (!is.null(result_folder) && !is.null(ssEnv) && length(ssEnv) > 0 &&
      !is.null(ssEnv$result_folder) &&
      !.core_same_folder(ssEnv$result_folder, result_folder))
    ssEnv <- NULL

  session_folder <- file.path(result_folder,"/Log/session_info.rds")
  # try from file
  if (is.null(ssEnv) | length(ssEnv)==0)
    if((!is.null(result_folder) & file.exists(session_folder)))
      ssEnv <- readRDS( file.path(session_folder))

  # try from doFuture
  # if (is.null(ssEnv) | length(ssEnv)==0)
  #   ssEnv <- getOption("ssEnv")

  if (is.null(ssEnv) | length(ssEnv)==0) {
    if (!is.null(result_folder))
      return(list())
    stop("ERROR: core_get_session_info called without result folder!")
  }

  return(ssEnv)

}


# Two paths name the same folder when they differ only by a trailing separator,
# by "." components, or by being one relative and one absolute.
.core_same_folder <- function(a, b)
{
  tidy <- function(p) sub("[/\\\\]+$", "",
                          normalizePath(as.character(p), winslash = "/", mustWork = FALSE))
  identical(tidy(a), tidy(b))
}
