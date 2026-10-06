#' Write the study sample sheet of the run (internal)
#'
#' this function used to append the per-sample burden
#' (`MUTATIONS_HYPER`, `DELTAS_HYPO`, ...) and `PROBES_COUNT` to the sample
#' sheet. Those columns were aggregated over the `AREA == "POSITION"` keys -
#' i.e. over every probe, with no genomic filter - which is exactly the
#' `SAMPLE` scope of the statistics sibling. They now live in
#' the `SCOPE = SAMPLE` artefacts, one row tall, and are joined back on read by
#' [sem_study_summary_get()], so the sample sheet stays the description of the
#' study.
#'
#' The taxonomy unification brought **one** of them back here: `N_PROBES`. It is a property of the
#' imputation - how many positions of that sample survived the treatment of
#' missing values - not an aggregation of a marker over a region class, so it
#' belongs with the descriptive properties of the sample rather than in the
#' taxonomy.
#'
#' @return Invisibly the path of the file written.
#' @keywords internal
#' @noRd
sem_study_summary_total <- function()
{
  ssEnv <- core_get_session_info()
  study_summary <- sem_study_summary_get(with_sample_stats = FALSE)

  if (is.null(study_summary) || nrow(study_summary) == 0) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " No sample sheet available; SAMPLE_SHEET_RESULT not written.")
    return(invisible(NULL))
  }

  # N_PROBES is an imputation property of the sample, so it travels with
  # the sample sheet. Absent when the position pivot has not been produced yet,
  # which is not an error: the sheet is usable without it.
  n_probes <- .sem_sample_n_probes()
  if (!is.null(n_probes) && "Sample_ID" %in% colnames(study_summary)) {
    study_summary <- study_summary[, colnames(study_summary) != "N_PROBES", drop = FALSE]
    study_summary <- merge(study_summary, n_probes, by = "Sample_ID", all.x = TRUE)
  }

  # io_file_path_build() uppercases via core_name_cleaning() - on-disk name is
  # SAMPLE_SHEET_RESULT.csv. Linux ext4 is case-sensitive; readers that
  # hard-code the path must use the uppercase form. See io_file_path_build().
  summary_file <- io_file_path_build( ssEnv$result_folderData, "sample_sheet_result","csv")
  utils::write.csv2(study_summary, summary_file, row.names = FALSE)

  invisible(summary_file)
}
