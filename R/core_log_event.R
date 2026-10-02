core_log_event <- function(...)
{
  # Robustness: core_log_event must be callable even when no session has been
  # initialised (e.g. anno_area_granges_build() invoked standalone, outside of a
  # full semseeker() run) AND when an earlier session left a stale
  # session_folder pointing at a tempdir that has since been unlink()-ed
  # (a common pattern in testthat runs). Silently no-op in both cases
  # instead of crashing on file().
  ssEnv <- tryCatch(core_get_session_info(), error = function(e) NULL)
  if (is.null(ssEnv) || length(ssEnv) == 0L ||
      is.null(ssEnv$session_folder) ||
      !dir.exists(ssEnv$session_folder))
    return(invisible(NULL))

  # append core_log_event to log file
  log_events <- list(...)
  log_event_to_save <- ""
  for (i in seq_along(log_events))
  {
    log_event_to_save <- paste0(log_event_to_save, log_events[i], sep =" ")
  }

  log_event_to_save <- gsub("  "," ", log_event_to_save)
  log_event_to_save <- gsub("  "," ", log_event_to_save)
  log_event_to_save <- gsub("  "," ", log_event_to_save)
  log_event_to_save <- gsub("  "," ", log_event_to_save)

  # Read session verbosity once: drives BOTH the augmentation policy
  # (below) and the final print filter (further down).
  if (is.null(ssEnv$verbosity))
    verbosity <- 4
  else
    verbosity <- as.numeric(ssEnv$verbosity)

  # Auto-augment messages with memory + CPU usage. AI-061+ (2026-06-09):
  # expanded augmentation policy on user request:
  #   - DEBUG: always augmented (legacy behaviour)
  #   - WARNING / ERROR: always augmented (mem/cpu snapshot at the
  #     moment a problem is logged is the most actionable diagnostic)
  #   - INFO: augmented only when the session is in DEBUG verbosity
  #     (verbosity == 4) - keeps INFO logs compact at normal levels
  #     but turns them into a full trace under DEBUG.
  augment_now <- grepl("^DEBUG",   log_event_to_save) ||
                 grepl("^WARNING", log_event_to_save) ||
                 grepl("^ERROR",   log_event_to_save) ||
                 (grepl("^INFO",   log_event_to_save) && verbosity == 4)
  if (augment_now) {
    mem_mb <- tryCatch(
      sum(gc(verbose = FALSE, full = FALSE)[, "(Mb)"]),
      error = function(e) NA_real_
    )
    cpu_pct <- if (.Platform$OS.type == "unix") {
      tryCatch(
        as.numeric(trimws(system2(
          "ps", args = c("-p", Sys.getpid(), "-o", "%cpu="),
          stdout = TRUE)[1])),
        error = function(e) NA_real_
      )
    } else NA_real_
    log_event_to_save <- paste0(log_event_to_save,
      " [mem=", round(mem_mb, 1), "MB cpu=", round(cpu_pct, 1), "%]")
  }

  file_name <- paste(as.character(Sys.info()["nodename"]),"_session_output.log", sep="")
  log_file <- file.path(ssEnv$session_folder,file_name)

  log_event_to_print <- log_event_to_save
  if (grepl("^BANNER", log_event_to_save))
  {
    log_event_to_print_dash <- "##############################################################################"
    log_event_to_print <- paste0( log_event_to_print_dash, "\n",
      gsub("BANNER: ", "", log_event_to_print), "\n",
      log_event_to_print_dash, '\n')
  }

  if (grepl("ERROR:", log_event_to_save))
  {
    log_event_to_print_dash <- "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    log_event_to_print <- paste0( log_event_to_print_dash, "\n",  log_event_to_print  , "\n",log_event_to_print_dash)
  }

  # `verbosity` already initialised at the top of the function for the
  # augmentation policy - reuse the same value here.

  cat(log_event_to_save, "\n", file = log_file, append = TRUE)

  file_name <- "console_session_output.log"
  log_file <- file.path(ssEnv$session_folder,file_name)
  cat(log_event_to_save, "\n", file = log_file, append = TRUE)

  if(grepl("^JOURNAL", log_event_to_save))
  {
    file_name <- "lab_journal.log"
    journal_file <- file.path(ssEnv$session_folder,file_name)
    log_event_to_save <- gsub("^JOURNAL: ", "", log_event_to_save)
    cat(log_event_to_save, "\n", file = journal_file, append = TRUE)
  }

  # AI-061+ (2026-06-09): rewritten with explicit parens and a verbosity
  # ladder. The previous form had operator-precedence inconsistencies
  # between the four branches (the verbosity == 1 line let BANNER through
  # by accident via && / || precedence, while the verbosity == 2/3 lines
  # gated BANNER inside the verbosity check) and printed every BANNER
  # twice on verbosity == 2/3.
  #
  # New rule:
  #   BANNER  → always printed (section separators the user wants to see)
  #   ERROR   → always printed (critical, regardless of verbosity)
  #   WARNING → verbosity >= 2
  #   INFO    → verbosity >= 3
  #   DEBUG   → verbosity >= 4
  # testthat::is_testing() suppresses all console output unconditionally.

  if (!testthat::is_testing()) {
    if (grepl("^BANNER",  log_event_to_save) ||
        grepl("^ERROR",   log_event_to_save)) {
      message(log_event_to_print)
      log_event_to_print <- ""
    } else if (verbosity >= 2 && grepl("^WARNING", log_event_to_save)) {
      message(log_event_to_print)
      log_event_to_print <- ""
    } else if (verbosity >= 3 && grepl("^INFO",    log_event_to_save)) {
      message(log_event_to_print)
      log_event_to_print <- ""
    } else if (verbosity >= 4 && grepl("^DEBUG",   log_event_to_save)) {
      message(log_event_to_print)
      log_event_to_print <- ""
    }
  }

}
