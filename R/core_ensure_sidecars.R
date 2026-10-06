#' Ensure every pivot parquet has its JSON sidecar
#'
#' Walks the \code{Pivots/} subtree under \code{data_root} (or any directory
#' tree containing \code{*.parquet} files) and writes a sidecar
#' \code{<base>_meta.json} for every parquet that is missing one. The sidecar
#' content is produced by \code{\link{core_pivot_sidecar_write}} from the current
#' \code{ssEnv} (genome_build, tech, semseeker_version, timestamp).
#'
#' Single point of responsibility for sidecar materialisation: callers that
#' write pivots (\code{anno_create_position_pivots}, \code{io_read_pivot},
#' \code{sem_core}, \code{core_recover}, ...) do not need
#' to write sidecars themselves - they simply call this function once at the
#' end of their work and let it scan for parquets without a sidecar.
#'
#' @param data_root character. Root folder to scan recursively. Defaults to
#'   \code{ssEnv$result_folderData} (i.e. the \code{Data/} folder of the
#'   current session).
#'
#' @return Invisibly, the character vector of sidecar paths that were written
#'   (possibly empty).
#'
#' @keywords internal
#' @noRd
core_ensure_sidecars <- function(data_root = NULL) {

  if (is.null(data_root)) {
    ssEnv <- core_get_session_info()
    data_root <- ssEnv$result_folderData
  }
  if (!dir.exists(data_root)) return(invisible(character(0)))

  parquets <- list.files(data_root, pattern = "\\.parquet$",
                         recursive = TRUE, full.names = TRUE)
  if (length(parquets) == 0L) return(invisible(character(0)))

  written <- character(0)
  for (p in parquets) {
    sidecar <- paste0(sub("\\.parquet$", "", p), "_meta.json")
    if (!file.exists(sidecar)) {
      core_pivot_sidecar_write(p)
      written <- c(written, sidecar)
    }
  }

  if (length(written) > 0L)
    core_log_event("INFO: ", Sys.time(),
              " core_ensure_sidecars wrote ", length(written),
              " missing sidecar(s) under ", data_root)

  invisible(written)
}
