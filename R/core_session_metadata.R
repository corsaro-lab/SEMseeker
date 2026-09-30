#' Write session provenance metadata
#'
#' Creates (or overwrites) \code{session_metadata.json} in
#' \code{result_folder}.  The JSON is human-readable and does not require R to
#' parse, making it suitable for pipeline auditing and automated compatibility
#' checks across studies.
#'
#' Fields written:
#' \describe{
#'   \item{genome_build}{Reference assembly used (\code{"hg19"} / \code{"hg38"}
#'     / \code{"mm10"}).}
#'   \item{tech}{Methylation technology detected or declared by the user
#'     (\code{"K850"}, \code{"K450"}, \code{"K27"}, \code{"WGBS"},
#'     \code{"LONGREAD"}, or \code{""} if not yet determined at write time).}
#'   \item{semseeker_version}{Package version string.}
#'   \item{created}{ISO-8601 timestamp of the run.}
#'   \item{sample_n}{Total number of samples in the run.}
#' }
#'
#' @param result_folder character. Path to the SEMseeker result folder.
#' @param sample_n integer. Total number of samples across all batches
#'   (default \code{0L}).
#'
#' @return Invisibly returns the metadata list that was serialised to JSON.
#' @keywords internal
core_session_metadata_write <- function(result_folder, sample_n = 0L) {
  ssEnv <- core_get_session_info()

  genome_build <- if (!is.null(ssEnv$genome_build) && nzchar(ssEnv$genome_build))
    as.character(ssEnv$genome_build) else "hg19"
  tech <- if (!is.null(ssEnv$tech) && nzchar(ssEnv$tech))
    as.character(ssEnv$tech) else ""

  meta <- list(
    genome_build      = genome_build,
    tech              = tech,
    semseeker_version = as.character(utils::packageVersion("SEMseeker")),
    created           = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    sample_n          = as.integer(sample_n)
  )

  out_path <- file.path(result_folder, "session_metadata.json")
  writeLines(jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE), out_path)
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " [session_metadata] Written: ", out_path)
  invisible(meta)
}


#' Write a pivot sidecar metadata file
#'
#' Writes a small JSON file alongside a pivot file (parquet or csv.gz) that
#' records the \code{genome_build} and \code{tech} of the session that created
#' it.  The sidecar path is constructed by appending \code{_meta.json} to the
#' pivot base name (before the extension).
#'
#' @param pivot_path character. Full path of the pivot file
#'   (e.g. \code{.../Pivots/MUTATIONS/MUTATIONS_HYPER_GENE_TSS1500_hg19.parquet}).
#'
#' @return Invisibly \code{NULL}.
#' @keywords internal
core_pivot_sidecar_write <- function(pivot_path) {
  ssEnv <- core_get_session_info()

  genome_build <- if (!is.null(ssEnv$genome_build) && nzchar(ssEnv$genome_build))
    as.character(ssEnv$genome_build) else "hg19"
  tech <- if (!is.null(ssEnv$tech) && nzchar(ssEnv$tech))
    as.character(ssEnv$tech) else ""

  # Build sidecar path: strip last extension, append _meta.json
  # e.g. foo_hg19.parquet  →  foo_hg19_meta.json
  base_no_ext <- sub("\\.[^.]+$", "", pivot_path)
  sidecar_path <- paste0(base_no_ext, "_meta.json")

  meta <- list(
    genome_build      = genome_build,
    tech              = tech,
    semseeker_version = as.character(utils::packageVersion("SEMseeker")),
    pivot_file        = basename(pivot_path),
    created           = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )
  writeLines(jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE), sidecar_path)
  invisible(NULL)
}


#' Check compatibility of SEMseeker sessions before meta-analysis
#'
#' Reads \code{session_metadata.json} from each path in \code{session_list}
#' and enforces provenance rules before combining results across studies:
#' \enumerate{
#'   \item \strong{Stops} if \code{genome_build} differs - coordinates from
#'     different assemblies are physically incomparable and would produce
#'     silently wrong intersection results.
#'   \item \strong{Warns} if \code{tech} differs - cross-array meta-analysis
#'     (e.g. K450 + K850) is statistically valid on the probe intersection but
#'     must be intentional.
#' }
#' A missing \code{session_metadata.json} raises a warning (legacy runs
#' without provenance data are tolerated but flagged).
#'
#' @param session_list character vector. Paths to SEMseeker result folders to
#'   compare (must be at least two for a meaningful comparison).
#'
#' @return Invisibly returns a \code{data.frame} with one row per session
#'   containing the \code{folder}, \code{genome_build}, and \code{tech}
#'   fields parsed from each metadata file.
#'
#' @examples
#' # Two fake session folders with provenance metadata
#' d1 <- file.path(tempdir(), "study_A"); dir.create(d1, showWarnings = FALSE)
#' d2 <- file.path(tempdir(), "study_B"); dir.create(d2, showWarnings = FALSE)
#' writeLines(jsonlite::toJSON(list(genome_build = "hg19", tech = "K450"),
#'                             auto_unbox = TRUE),
#'            file.path(d1, "session_metadata.json"))
#' writeLines(jsonlite::toJSON(list(genome_build = "hg19", tech = "K450"),
#'                             auto_unbox = TRUE),
#'            file.path(d2, "session_metadata.json"))
#' core_check_session_compatibility(c(d1, d2))
#' @export
core_check_session_compatibility <- function(session_list) {
  if (length(session_list) < 2L) {
    core_log_event("INFO: [core_check_session_compatibility] Only one session - nothing to compare.")
    return(invisible(NULL))
  }

  meta_list <- lapply(session_list, function(folder) {
    json_path <- file.path(folder, "session_metadata.json")
    if (!file.exists(json_path)) {
      core_log_event("WARNING: [core_check_session_compatibility]",
                " No session_metadata.json in: ", folder,
                " - provenance cannot be verified.")
      return(data.frame(folder        = folder,
                        genome_build  = NA_character_,
                        tech          = NA_character_,
                        stringsAsFactors = FALSE))
    }
    parsed <- jsonlite::fromJSON(json_path)
    gb <- parsed[["genome_build"]]
    tc <- parsed[["tech"]]
    data.frame(
      folder       = folder,
      genome_build = if (!is.null(gb)) as.character(gb) else NA_character_,
      tech         = if (!is.null(tc)) as.character(tc) else NA_character_,
      stringsAsFactors = FALSE
    )
  })

  meta_df <- do.call(rbind, meta_list)

  # ── genome_build: STOP if any two known builds differ ─────────────────────
  builds_known <- meta_df$genome_build[!is.na(meta_df$genome_build)]
  if (length(unique(builds_known)) > 1L) {
    stop(
      "[core_check_session_compatibility] INCOMPATIBLE sessions: genome_build differs.\n",
      "  Combining coordinates from different genome assemblies produces wrong results.\n",
      "  Found: ", paste(unique(builds_known), collapse = ", "), "\n",
      "  Sessions:\n",
      paste0("    ", meta_df$folder, "  [", meta_df$genome_build, "]",
             collapse = "\n")
    )
  }

  # ── tech: WARN if different ───────────────────────────────────────────────
  techs_known <- meta_df$tech[!is.na(meta_df$tech) & nzchar(meta_df$tech)]
  if (length(unique(techs_known)) > 1L) {
    warning(
      "[core_check_session_compatibility] Technologies differ across sessions: ",
      paste(unique(techs_known), collapse = ", "), ".\n",
      "  Cross-array meta-analysis is valid on the probe intersection ",
      "but must be intentional.\n",
      "  Proceed only if you have verified that probe sets overlap."
    )
  }

  if (length(builds_known) > 0L)
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " [core_check_session_compatibility] All sessions compatible.",
              " genome_build=", unique(builds_known)[1L],
              " tech=", paste(unique(techs_known), collapse = "+"))

  invisible(meta_df)
}
