core_close_env <- function()
{

  ssEnv <- tryCatch(core_get_session_info(), error = function(e) NULL)

  # Closing what is not open is not an error. Once a close empties the session
  # this becomes reachable, and it is reached often: these tests close in the
  # body and again from on.exit(), and a caller is entitled to close defensively.
  if (is.null(ssEnv) || length(ssEnv) == 0L)
    return(invisible(NULL))

  if (ssEnv$showprogress)
    progressr::handlers()

  # build all the folder tree in result_folder

  core_remove_empty_folders(ssEnv$result_folder)

  future::plan( future::sequential)
  unlink(ssEnv$temp_folder,recursive=TRUE)
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Job Completed !")
  core_log_event("DEBUG: --------------------------------------------------------------")

  # A close that leaves the session in memory closes nothing. The session is what
  # identifies an analysis, so the next read without a folder would be answered
  # with this one, and a caller that forgot to initialise would receive a stale
  # session instead of a failure. Emptied to the same value .onLoad() leaves.
  #
  # Last, and not earlier: the two log events above resolve their file through the
  # session, and core_log_event() swallows that lookup error, so emptying before
  # them would not fail. It would silently stop writing the log.
  assign("ssEnv", list(), envir = .pkgglobalenv)

  invisible(NULL)
}


core_remove_empty_folders <- function(path) {
  files <- list.files(path, full.names = TRUE)

  for (file in files) {
    if (file.info(file)$isdir) {
      core_remove_empty_folders(file)  # Recursively check subdirectories

      # After the check, remove the directory if empty
      if (length(list.files(file)) == 0) {
        unlink(file, recursive = TRUE)
        # cat("Removed empty directory:", file, "\n")
      }
    }
  }
}
