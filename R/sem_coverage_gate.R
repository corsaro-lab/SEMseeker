#' Mandatory coverage check run before every SEM analysis (internal)
#'
#' Coverage used to be an opt-in report: any pipeline path could go
#' straight to the SEM rebuild on signal data that has little or no overlap
#' with the reference annotation (Nanopore input against an Illumina manifest,
#' a batch with a probe mismatch, the wrong genome build). The run then
#' produced thresholds computed on a handful of probes and every downstream
#' number was quietly wrong.
#'
#' This gate always runs, in this order:
#' \enumerate{
#'   \item the coverage charts (\code{\link{sem_coverage_analysis}}) are
#'     generated on EVERY run, whatever the verdict - they are the artefact the
#'     analyst reads to understand what was covered;
#'   \item a global coverage percentage (input positions found in the reference
#'     annotation) is logged as a BANNER and written to a JSON sidecar for
#'     later audit;
#'   \item the run STOPS when coverage is below \code{coverage_minimum}
#'     (default 80\%), naming the numbers involved.
#' }
#'
#' Coordinate-based technologies (WGBS / long reads) derive their annotation
#' from the positions themselves, so the overlap is 100\% by construction and
#' the threshold is not applied; the charts are still attempted and the JSON
#' sidecar still written, with \code{gate = "skipped"}.
#'
#' @param observed_probes character. Probe / position identifiers of the input.
#' @param tech character. Technology label; defaults to the session value.
#'
#' @return Invisibly the list of coverage metrics.
#' @keywords internal
#' @noRd
sem_coverage_gate <- function(observed_probes, tech = NULL) {

  ssEnv <- core_get_session_info()
  if (is.null(tech)) tech <- ssEnv$tech
  tech <- if (is.null(tech)) "" else as.character(tech)

  observed_probes <- unique(as.character(observed_probes))
  n_input <- length(observed_probes)

  minimum <- suppressWarnings(as.numeric(ssEnv$coverage_minimum))
  if (length(minimum) != 1L || is.na(minimum)) minimum <- 80

  # ---- (1) charts ALWAYS, whatever the verdict -----------------------------
  # A plotting failure must never take down a SEM run: log it and carry on to
  # the numeric gate, which is the part that protects the results.
  # Report on the FULL annotation set, not on ssEnv$keys_areas_subareas: the
  # latter is filtered down to the areas the run analyses (POSITION only by
  # default), which would leave the charts empty on a default run.
  keys_all <- ssEnv$default$keys_areas_subareas_default
  tryCatch(
    sem_coverage_analysis(observed_probes, keys = keys_all),
    error = function(e)
      core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
                " Coverage charts could not be generated: ",
                conditionMessage(e))
  )

  # ---- (2) global coverage -------------------------------------------------
  coord_tech <- tech %in% c("WGBS", "LONGREAD")
  n_covered <- NA_integer_
  pct <- NA_real_
  verdict <- "skipped"

  if (coord_tech) {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " Coverage gate not applicable to tech '", tech,
              "': positions carry their own coordinates, overlap is 100% by construction.")
    n_covered <- n_input
    pct <- 100
  } else {
    reference <- tryCatch(anno_probe_features_get("PROBE"),
                          error = function(e) NULL)
    if (is.null(reference) || nrow(reference) == 0) {
      core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                " Coverage gate could not load the probe annotation for tech '",
                tech, "' / build '", ssEnv$genome_build,
                "': the gate is not enforced for this run.")
    } else {
      n_covered <- sum(observed_probes %in% as.character(reference$PROBE))
      pct <- if (n_input > 0) 100 * n_covered / n_input else 0
      verdict <- if (pct >= minimum) "pass" else "fail"
    }
  }

  metrics <- list(
    tech             = tech,
    genome_build     = if (is.null(ssEnv$genome_build)) NA_character_ else ssEnv$genome_build,
    input_positions  = n_input,
    covered_positions = n_covered,
    coverage_percent = if (is.na(pct)) NA_real_ else round(pct, 2),
    coverage_minimum = minimum,
    gate             = verdict
  )

  core_log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"),
            " Coverage - input_positions=", n_input,
            " | covered_by_reference=", n_covered,
            " | coverage=", if (is.na(pct)) "NA" else paste0(round(pct, 2), "%"),
            " | minimum=", minimum, "%",
            " | gate=", verdict)

  # ---- (3) JSON sidecar for later audit ------------------------------------
  tryCatch({
    sidecar <- io_file_path_build(ssEnv$result_folderData, "coverage_gate", "json")
    writeLines(jsonlite::toJSON(metrics, auto_unbox = TRUE, pretty = TRUE), sidecar)
  }, error = function(e)
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " Coverage sidecar not written: ", conditionMessage(e)))

  # ---- (4) verdict ---------------------------------------------------------
  if (identical(verdict, "fail")) {
    msg <- paste0(
      "Coverage gate failed: only ", n_covered, " of ", n_input,
      " input positions (", round(pct, 2), "%) are present in the reference ",
      "annotation for tech '", tech, "' / genome_build '", ssEnv$genome_build,
      "', below the required ", minimum, "%. This usually means the input does ",
      "not match the declared technology or genome build. Fix the input, or ",
      "lower the threshold explicitly with coverage_minimum = <value> if a ",
      "partial overlap is expected (e.g. a deliberate cross-technology run)."
    )
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " ", msg)
    stop(msg)
  }

  invisible(metrics)
}
