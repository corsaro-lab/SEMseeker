#' Bayesian posterior probability analysis of SEMseeker mutations and lesions
#'
#' Computes P(case | epimutated) via empirical Bayes for each marker and
#' genomic area, providing a probabilistic complement to the frequentist
#' association analysis.
#'
#' @param result_folder character. Path to the SEMseeker result folder.
#' @param independent_variable character. Sample sheet column defining the
#'   case/control grouping variable (default \code{"Sample_Group"}).
#' @param maxResources numeric. Maximum percentage of CPU cores to use
#'   (default 90).
#' @param parallel_strategy character. Parallelisation backend; one of
#'   \code{"multisession"}, \code{"sequential"}, \code{"cluster"}
#'   (default \code{"multisession"}).
#'   Asking for \code{"multicore"} is accepted and converted to
#'   \code{"multisession"}: it means fork(), which is unsafe with this
#'   package's native thread pool on every platform that offers it, and
#'   absent on Windows. The conversion is logged.
#' @param bayes_case_threshold numeric. Minimum P(case | epimutated) required
#'   to report a hit (default 0.9).
#' @param bayes_control_threshold numeric. Maximum P(control | epimutated)
#'   allowed to report a hit (default 0.1).
#' @param ... Additional arguments passed to \code{core_init_env()}.
#'
#' @return Invisibly \code{NULL}. Bayesian posterior probability tables
#'   (\code{bayes_analysis_*.csv}) are written to the \code{Euristic/}
#'   sub-folder of \code{result_folder}.
#' @importFrom doRNG %dorng%
#' @examples
#' result_dir <- tempdir()
#' \dontrun{
#' assoc_bayes_analysis(
#'   result_folder        = "~/semseeker_results/",
#'   independent_variable = "Sample_Group"
#' )
#' }
#' @export
assoc_bayes_analysis <- function(
    result_folder,
    independent_variable    = "Sample_Group",
    maxResources            = 90,
    parallel_strategy       = "multisession",
    bayes_case_threshold    = 0.9,   # A-09 fix 9: exposed as parameter
    bayes_control_threshold = 0.1,   # A-09 fix 9: exposed as parameter
    ...)
{
  markers <- c("MUTATIONS", "LESIONS")
  ssEnv <- core_init_env(
    result_folder     = result_folder,
    maxResources      = maxResources,
    parallel_strategy = parallel_strategy,
    start_fresh       = FALSE,
    markers           = markers,
    ...
  )

  arguments       <- list(...)
  areas_selection <- if (!is.null(arguments[["areas_selection"]])) arguments$areas_selection else c()

  study_summary <- sem_study_summary_get()
  if (independent_variable == "Sample_Group")
    study_summary <- study_summary[, c("Sample_Group", "Sample_ID")]
  else
    study_summary <- study_summary[, c("Sample_Group", "Sample_ID", independent_variable)]

  # ── Outer loop: one pass per marker (MUTATIONS, LESIONS) ──────────────────
  for (a in seq_along(markers)) {   # A-09 fix 1: seq_along, not length()

    fileNameResults <- io_file_path_build(
      baseFolder      = ssEnv$result_folderEuristic,
      detailsFilename = c(markers[a], "assoc_bayes_analysis"),   # A-09 fix 7: typo
      extension       = "csv"
    )

    localKeys_1 <- ssEnv$keys_areas_subareas_markers_figures
    keys <- localKeys_1[localKeys_1$MARKER == markers[a], ]

    # ── Resume: skip combos already written to disk ────────────────────────
    results <- data.frame()   # initialize; possibly populated from existing file below
    if (file.exists(fileNameResults)) {
      results    <- utils::read.csv2(fileNameResults, header = TRUE)
      done_keys  <- unlist(apply(
        unique(results[, c("MARKER", "FIGURE", "AREA", "SUBAREA")]),
        1, function(x) paste(x, collapse = "_")
      ))
      todo_keys  <- unlist(apply(
        keys[, c("MARKER", "FIGURE", "AREA", "SUBAREA")],
        1, function(x) paste(x, collapse = "_")
      ))
      keys <- keys[!(todo_keys %in% done_keys), ]
    }

    if (nrow(keys) == 0) next

    # ── Inner loop: one pass per area/subarea/figure combo ─────────────────
    for (k in seq_len(nrow(keys))) {

      key           <- keys[k, ]
      pivot_filename <- io_pivot_file_name(key$MARKER, key$FIGURE, key$AREA, key$SUBAREA)

      if (!file.exists(pivot_filename)) next

      core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
        " [assoc_bayes_analysis] Reading pivot: ", pivot_filename)

      # A-09 fix 2: use pivot_filename (variable), not io_pivot_file_name (function)
      pivot <- readr::read_delim(
        pivot_filename,
        col_types      = readr::cols(.default = readr::col_double(),
                                     AREA     = readr::col_character()),
        show_col_types = FALSE,
        progress       = FALSE
      )

      # A-09 fix 3: assign pivot → tempDataFrame before using it
      tempDataFrame          <- pivot
      row.names(tempDataFrame) <- tempDataFrame$AREA

      if (length(areas_selection) > 0)
        tempDataFrame <- tempDataFrame[tempDataFrame[, 1] %in% areas_selection, ]

      tempDataFrame <- tempDataFrame[, -1]          # drop AREA name column
      if (is.null(dim(tempDataFrame)))   next
      if (plyr::empty(tempDataFrame) || nrow(tempDataFrame) == 0) next

      tempDataFrame           <- t(tempDataFrame)   # samples × areas
      tempDataFrame           <- as.data.frame(tempDataFrame)
      tempDataFrame$Sample_ID <- rownames(tempDataFrame)
      tempDataFrame           <- merge(study_summary, tempDataFrame,
                                       by = "Sample_ID", all.x = TRUE)

      # A-09 fix 3 (cont.): remove only Sample_ID (merge key), keep everything else.
      # The original c(-1,-3) accidentally dropped the first data column or the
      # independent_variable column depending on the session configuration.
      tempDataFrame <- tempDataFrame[, colnames(tempDataFrame) != "Sample_ID",
                                     drop = FALSE]

      # A-09 fix 4: column reference, not string literal - "x" != "y" is always TRUE
      tempDataFrame <- subset(tempDataFrame, Sample_Group != "Reference")
      tempDataFrame <- subset(tempDataFrame, Sample_Group != 0)
      tempDataFrame <- as.data.frame(tempDataFrame)

      if (nrow(tempDataFrame) == 0) next

      tempDataFrame$Sample_Group <- tempDataFrame$Sample_Group == "Case"

      tempDataFrame[is.na(tempDataFrame)] <- 0

      # Sanitise column names (spaces/dashes/etc. → underscore)
      colnames(tempDataFrame) <- gsub("[[:space:]\\-:/']", "_", colnames(tempDataFrame))

      phenotype <- as.logical(tempDataFrame[, independent_variable])
      n_case    <- sum(phenotype,  na.rm = TRUE)
      n_control <- sum(!phenotype, na.rm = TRUE)

      if (n_case == 0 || n_control == 0) {
        core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
          " [assoc_bayes_analysis] Skipping ", key$MARKER, "/", key$FIGURE,
          "/", key$AREA, "/", key$SUBAREA,
          " - n_case=", n_case, " n_control=", n_control)
        next
      }

      if (ssEnv$showprogress)
        progress_bar <- progressr::progressor(along = seq_len(ncol(tempDataFrame)))
      else
        progress_bar <- ""

      # Identify data columns (skip grouping / phenotype meta-columns)
      meta_cols      <- unique(c("Sample_Group", independent_variable))
      data_col_idx   <- which(!colnames(tempDataFrame) %in% meta_cols)

      var_to_export  <- c("tempDataFrame", "ssEnv", "progress_bar",
                          "phenotype", "keys", "k", "n_case", "n_control",
                          "bayes_case_threshold", "bayes_control_threshold",
                          "independent_variable", "meta_cols")

      # A-09 fix 8: loop variable renamed col_idx (was 'c', which shadows base c())
      results_temp <- foreach::foreach(
        col_idx  = data_col_idx,
        .combine = rbind,
        .export  = var_to_export
      ) %dorng% {
        # AI-056: workers must NOT saveRDS on every iteration (the worker-copy
        # of ssEnv can't propagate to master anyway, and the per-iter disk
        # write is the AI-041 performance trap).
        core_update_session_info(ssEnv, save_to_disk = FALSE)

        area <- names(tempDataFrame)[col_idx]
        if (ssEnv$showprogress)
          progress_bar(sprintf("genomic area: %s",
            stringr::str_pad(area, 20, side = "left", pad = " ")))

        epimutated <- as.numeric(tempDataFrame[, col_idx]) != 0

        # Bayes theorem:  P(A|B) = P(B|A) * P(A) / P(B)
        # A = being a Case;  B = being epimutated in this area
        P_B <- sum(epimutated) / length(epimutated)
        if (P_B == 0) return(NULL)   # no epimutations → skip (avoid 0/0)

        P_A      <- n_case    / (n_case + n_control)
        P_B_A    <- sum(epimutated &  phenotype) / sum( phenotype)
        P_A_B    <- (P_B_A * P_A) / P_B    # P(Case    | Epimutated)

        P_notA   <- n_control / (n_case + n_control)
        P_B_notA <- sum(epimutated & !phenotype) / sum(!phenotype)
        P_notA_B <- (P_B_notA * P_notA) / P_B  # P(Control | Epimutated)

        # A-09 fix 9: configurable thresholds (were hardcoded 0.9 / 0.1)
        if (P_A_B >= bayes_case_threshold && P_notA_B < bayes_control_threshold)
          data.frame(
            MARKER                               = as.character(keys[k, "MARKER"]),
            FIGURE                               = as.character(keys[k, "FIGURE"]),
            AREA                                 = as.character(keys[k, "AREA"]),
            SUBAREA                              = as.character(keys[k, "SUBAREA"]),
            AREA_OF_TEST                         = area,
            P_to_be_Case_cond_to_be_Epimutated   = P_A_B,
            P_to_be_Control_cond_to_be_Epimutated = P_notA_B
          )
      }

      # AI-056: post-foreach end-of-batch snapshot (matches AI-041 pattern).
      core_update_session_info(ssEnv, save_to_disk = TRUE)

      if (is.null(dim(results_temp))) next
      results_temp <- as.data.frame(results_temp)
      colnames(results_temp) <- c("MARKER", "FIGURE", "AREA", "SUBAREA",
                                  "AREA_OF_TEST",
                                  "P_to_be_Case_cond_to_be_Epimutated",
                                  "P_to_be_Control_cond_to_be_Epimutated")

      results <- if (nrow(results) > 0) plyr::rbind.fill(results, results_temp) else results_temp
      utils::write.csv2(x = results, file = fileNameResults, row.names = FALSE)
    }

    # ── End-of-marker cleanup: final file + filtered file ──────────────────
    # A-09 fix 5: exists() scoped to local env (same pattern as E-01)
    if (!exists("results", envir = environment(), inherits = FALSE)) next
    if (nrow(results) == 0) next

    results <- subset(results, results$AREA != "CHR")
    results <- results[, colSums(is.na(results)) != nrow(results), drop = FALSE]

    results$P_to_be_Case_cond_to_be_Epimutated    <-
      as.numeric(results$P_to_be_Case_cond_to_be_Epimutated)
    results$P_to_be_Control_cond_to_be_Epimutated <-
      as.numeric(results$P_to_be_Control_cond_to_be_Epimutated)

    utils::write.csv2(x = results, file = fileNameResults, row.names = FALSE)

    # A-09 fix 6: two DISTINCT max variables (was duplicate of Case)
    max_P_case    <- max(results$P_to_be_Case_cond_to_be_Epimutated,    na.rm = TRUE)
    max_P_control <- max(results$P_to_be_Control_cond_to_be_Epimutated, na.rm = TRUE)

    results_filtered <- subset(results,
      P_to_be_Case_cond_to_be_Epimutated    != 0 &
      P_to_be_Control_cond_to_be_Epimutated != 0 &
      P_to_be_Case_cond_to_be_Epimutated    == max_P_case
    )

    fileNameFiltered <- io_file_path_build(
      baseFolder      = ssEnv$result_folderEuristic,
      detailsFilename = c(markers[a], "filtered_bayes_analysis"),  # A-09 fix 7: typo
      extension       = "csv"
    )
    utils::write.csv2(x = results_filtered, file = fileNameFiltered, row.names = FALSE)

    rm("results", envir = environment())   # A-09 fix 5: scoped rm
  }

  core_close_env()
  invisible(NULL)
}
