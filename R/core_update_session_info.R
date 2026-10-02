#' Persist the session environment
#'
#' Stores ssEnv into `.pkgglobalenv` (in-memory cache, always) and optionally
#' onto disk as `session_info.rds` and `<YYYY-MM-DD>_session_info.rds` inside
#' `ssEnv$session_folder`.
#'
#' **Hot-path callers (inside foreach loops) MUST pass `save_to_disk = FALSE`**
#' to avoid hammering the disk with one 15 MB rds write per gene per worker.
#' The on-disk persistence is meant for end-of-job snapshots, not per-iteration
#' state syncs: the split exists because per-iteration writes were a
#' measured performance regression.
#'
#' @param ssEnv list. Session environment.
#' @param save_to_disk logical. If TRUE (default, backward-compatible) writes
#'   the session to disk as well as in-memory. If FALSE, only updates the
#'   in-memory cache - fast path for worker bodies inside `foreach \%dorng\%`.
#'
#' @return ssEnv, invisibly.
#' @keywords internal
core_update_session_info <- function(ssEnv, save_to_disk = TRUE)
{
  if (is.null(ssEnv) | length(ssEnv)==0)
    stop("DEBUG: I'm STOPPING HERE! You called update session info without ssEnv!")

  # save to environment (always; this is what core_get_session_info() reads first)
  assign("ssEnv", ssEnv, envir=.pkgglobalenv)

  if (save_to_disk) {
    saveRDS(ssEnv, file.path(ssEnv$session_folder,"session_info.rds"))
    today <- Sys.Date()
    saveRDS(ssEnv, file.path(ssEnv$session_folder,paste0(today, "_session_info.rds", sep="")))
  }

  return(ssEnv)
}
