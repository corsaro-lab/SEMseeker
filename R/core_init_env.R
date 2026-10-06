# ---------------------------------------------------------------------------
# Module-private helpers and tables used by core_init_env()
#
# Refactor history:
# - 2026-05-31 phase 1: collapsed 5 repeated patterns (defaults table,
#   folders table, .core_log_info, Filter, .core_pop_arg). 223 -> 159 LOC.
# - 2026-05-31 phase 2: extracted 7 setup phases as .init_env_*() helpers.
#   159 -> ~40 LOC orchestrator; every helper is under the 50-LOC threshold.
#
# Behaviour is preserved end-to-end; only the shape changes. Existing
# tests in test-0-init_env.R, test-0-set-env_variable.R, test-0-log-event.R
# pass without modification.
# ---------------------------------------------------------------------------

# Per-option defaults applied by .core_apply_defaults().
.SS_DEFAULTS <- list(
  verbosity              = list(value = 1),
  q_b_param              = list(value = data.frame("DELTAP_B"=4,"DELTARP_B"=4,"DELTAQ_Q"=4,"DELTARQ_Q"=4)),
  DELTAP_B               = list(value = 4),
  DELTARP_B              = list(value = 4),
  DELTAQ_Q               = list(value = 4),
  DELTARQ_Q              = list(value = 4),
  inpute                 = list(value = "none"),
  plot_format            = list(value = "png"),
  plot_resolution        = list(value = "print"),
  plot_resolution_ppi    = list(value = 600),
  alpha                  = list(value = 0.05),
  sex_chromosome_remove  = list(value = FALSE),
  opencl                 = list(value = FALSE),
  bonferroni_threshold   = list(value = 0.05),
  iqrTimes               = list(value = 3),
  LESIONS_BP             = list(value = 5000L),  # bp-based window radius. Default 5000 bp = 5 kbp (literature-aligned; review pending - Bock 2012, Jaffe 2012 bumphunter, Aryee 2014 minfi DMR).
  tech                   = list(value = ""),
  genome_build           = list(value = "hg19", choices = c("hg19","hg38","mm10","legacy")),
  showprogress           = list(value = FALSE),
  openai_api_key         = list(value = ""),
  multiple_test_adj      = list(value = "q", choices = c("BY","fdr","BH","bonferroni","q")),
  coverage_minimum       = list(value = 80),      # minimum % of input positions that must be present in the reference annotation; below it sem_coverage_gate() stops the run.
  bulk_population        = list(value = TRUE)    # vectorized population is the default (no per-sample bed dump). Set FALSE only to core_recover the legacy per-sample loop.
)

.SS_FOLDERS <- c(
  result_folderData       = "Data",
  result_folderChart      = "Chart",
  result_folderInference  = "Inference",
  result_folderEnrichment = "Enrichment",
  result_folderPhenotype  = "Phenotype",
  result_folderEuristic   = "Euristic",
  session_folder          = "Log"
)

.SS_FN_EXPORT <- c("sem_analyze_single_sample","sem_deltar_single_sample",
                   "io_dump_sample_as_bed_file", "sem_delta_single_sample",
                   "io_dir_check_and_create", "io_file_path_build",
                   "sem_analyze_single_sample_both", "anno_sort_by_chr_and_start",
                   "util_test_match_order", "sem_lesions_get", "sem_mutations_get")

.core_apply_defaults <- function(arguments, defaults) {
  for (key in names(defaults)) {
    d <- defaults[[key]]
    if (is.null(d$choices))
      arguments <- core_set_env_variable(arguments, key, d$value)
    else
      arguments <- core_set_env_variable(arguments, key, d$value, d$choices)
  }
  arguments
}

.core_log_info <- function(...) {
  core_log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y"), ...)
}

.core_pop_arg <- function(args, name, default) {
  val <- if (!is.null(args[[name]])) args[[name]] else default
  args[[name]] <- NULL
  list(value = val, args = args)
}

.core_init_env_silence_warnings <- function() {
  PKGs <- c("future","doRNG","doParallel","progressr","data.table","ggplot2","dplyr",
            "readr","readxl","stringr","tidyr","tibble","purrr","ggpubr","ggrepel","ggsci",
            "foreach","VennDiagram")
  invisible(lapply(PKGs, suppressWarnings(suppressMessages)))
  invisible(lapply(PKGs, suppressPackageStartupMessages))
}

.core_init_env_clean_args <- function(arguments) {
  if (length(arguments) > 0) {
    arguments <- lapply(arguments, function(x) if (is.character(x)) gsub(" ", "", x) else x)
    arguments <- lapply(arguments, function(x) if (is.character(x)) x[x != ""] else x)
    arguments <- arguments[vapply(arguments, function(x) length(x) > 0, logical(1))]
    arguments <- arguments[vapply(arguments, function(x) !is.null(x), logical(1))]
  }
  arguments[["areas_selection"]] <- NULL
  arguments
}

.core_init_env_bootstrap_session <- function(result_folder, start_fresh) {
  if (start_fresh) {
    unlink(result_folder, recursive = TRUE, force = TRUE)
    ssEnv <- list()
  } else if (dir.exists(result_folder)) {
    ssEnv <- core_get_session_info(result_folder)
  } else {
    ssEnv <- list()
  }
  ssEnv$session_id <- if (is.null(ssEnv$session_id)) 0 else ssEnv$session_id + 1
  ssEnv$session_folder <- io_dir_check_and_create(result_folder, c("Log"))
  ssEnv$seed <- 7658776
  core_update_session_info(ssEnv)
  ssEnv
}

.core_init_env_apply_computed_defaults <- function(arguments) {
  original_colors <- c('#b9e192', '#b3c7f7', '#f8b8d0','#f194b8', '#ffefb6', '#cfebb6','#b9ef92')
  original_colors <- rep(original_colors, 2)
  arguments <- core_set_env_variable(arguments, "color_palette", original_colors)
  darker_colors <- grDevices::adjustcolor(original_colors, alpha.f = 0.5)
  darker_colors <- c("blue","red","purple","green","yellow","orange","brown")
  arguments <- core_set_env_variable(arguments, "color_palette_darker", darker_colors)
  arguments <- core_set_env_variable(arguments, "cluster_workers", NULL)
  model_metrics <- toupper(as.vector(SEMseeker::metrics_properties$Metric))
  arguments <- core_set_env_variable(arguments, "model_metrics", model_metrics)
  arguments
}

.core_init_env_setup_paths <- function(ssEnv, result_folder) {
  tmp <- tempdir()
  ssEnv$temp_folder   <- paste(tmp, "/semseeker/",
                               stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),
                               sep = "")
  ssEnv$result_folder <- result_folder
  for (key in names(.SS_FOLDERS))
    ssEnv[[key]] <- io_dir_check_and_create(result_folder, .SS_FOLDERS[[key]])
  ssEnv
}

.core_init_env_setup_log_sink <- function(session_folder) {
  if (identical(Sys.getenv("SEMSEEKER_CHILD"), "1")) return(invisible())
  if (sink.number() != 0) sink(NULL)
  file_name <- paste(as.character(Sys.info()["nodename"]), "_session_output.log", sep = "")
  sink(file.path(session_folder, file_name), split = TRUE, append = TRUE)
  invisible()
}

.core_init_env_setup_progress <- function(showprogress) {
  if (!showprogress) return(invisible())
  if (exists("cli", mode = "function", inherits = TRUE)) return(invisible())
  if (testthat::is_testing()) return(invisible())
  handler_settings <- progressr::handlers()
  if (!("cli" %in% handler_settings$handler)) {
    progressr::handlers(global = TRUE)
    progressr::handlers("cli")
  }
  invisible()
}

.core_init_env_validate_args <- function(arguments) {
  arguments <- Filter(function(x) !is.null(x) && !identical(x, character(0)), arguments)
  if (length(arguments) == 0) {
    return(invisible())
  }

  # Build "name = value" pairs so the diagnostic shows WHICH argument is
  # unrecognised, not just its value. Previously the error printed only
  # the values (e.g. "ERROR: This options are not recognized: FALSE"),
  # which made typos like `phenolyser = FALSE` undebuggable without
  # diving into the source.
  unknown_names <- names(arguments)
  pairs <- vapply(seq_along(arguments), function(i) {
    n <- if (!is.null(unknown_names) && nzchar(unknown_names[i]))
           unknown_names[i] else "?"
    v <- paste(format(arguments[[i]]), collapse = ", ")
    if (nchar(v) > 60) v <- paste0(substr(v, 1, 57), "...")
    sprintf("%s = %s", n, v)
  }, character(1))

  # Identify the user-facing caller (semseeker / association_analysis /
  # enrichment_analysis / ...) by walking up the call stack until we
  # leave the .init_env_* internals. This makes the error point at the
  # SEMseeker entry point the user actually called.
  caller <- "SEMseeker"
  for (depth in seq_len(20)) {
    parent <- tryCatch(sys.call(-depth), error = function(e) NULL)
    if (is.null(parent)) break
    fn <- tryCatch(deparse(parent[[1]])[1], error = function(e) "")
    if (!grepl("^(\\.init_env|core_init_env|eval|tryCatch|do\\.call|withCallingHandlers|sys\\.call|local|source|withVisible)", fn)) {
      caller <- fn
      break
    }
  }

  msg <- sprintf(
    paste0("Unrecognised argument(s) in call to %s():\n  %s\n",
           "Check spelling and case. Common causes: typo ",
           "(e.g. 'phenolyser' vs 'phenolyzer'), case mismatch ",
           "(e.g. 'stringdb' vs 'STRINGdb'), or a parameter renamed in a ",
           "newer SEMseeker version."),
    caller, paste(pairs, collapse = "\n  ")
  )
  .core_log_info(msg)
  stop(msg, call. = FALSE)
}

.core_init_env_handle_dry_run <- function(ssEnv) {
  knitr::kable(as.data.frame(ssEnv$keys_areas_subareas_markers_figures),
               format = "pipe", caption = "Selection:")
  message(ssEnv$keys_areas_subareas_markers_figures)
  stop("INFO: Dry run is requested. Exiting now.")
}

.core_init_env_check_kwargs <- function(args) {
  tryCatch(
    { test_it <- args },
    error = function(cond) {
      core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
                " Function's arguments must be passed explicitily !")
      core_log_event(cond)
      stop("Function's arguments must be passed explicitily !")
    }
  )
}

.core_init_env_log_focus <- function(ssEnv) {
  .core_log_info(" I will focus on:",
            paste(unique(ssEnv$keys_markers_figures$MARKER), collapse = " ", sep = " "),
            " due to ",
            paste(unique(ssEnv$keys_markers_figures$FIGURE), collapse = " ", sep = " "),
            " of ",
            paste(unique(ssEnv$keys_areas_subareas_markers_figures$AREA),
                  collapse = " ", sep = " "))
}

#' init ssEnvonment
#'
#' @param result_folder where result of semseeker will be stored
#' @param maxResources percentage of how many available cores will be used
#'   (default 90 percent, rounded to lowest integer)
#' @param ... additional session options, including:
#'   \describe{
#'     \item{\code{parallel_strategy}}{parallelisation strategy for \pkg{future}:
#'       \code{"sequential"} (default), \code{"multisession"}, \code{"multicore"},
#'       \code{"cluster"}}
#'     \item{\code{genome_build}}{reference genome assembly: \code{"hg19"} (default,
#'       matches Illumina array annotation), \code{"hg38"} (GRCh38, typical for
#'       long-read / Nanopore data), \code{"mm10"} (mouse - requires C-05).
#'       Stored in \code{ssEnv$genome_build} and written to session provenance
#'       metadata (C-06).}
#'     \item{\code{tech}}{override technology auto-detection: \code{"K850"},
#'       \code{"K450"}, \code{"K27"}, \code{"WGBS"}, \code{"LONGREAD"}.
#'       Required for long-read data because LONGREAD cannot be distinguished from
#'       WGBS by probe-ID pattern alone.  Example:
#'       \code{core_init_env(folder, tech = "LONGREAD", genome_build = "hg38")}}
#'   }
#'
#' @return the working ssEnvonment
core_init_env <- function(result_folder, maxResources = 90, ...) {
 
  gc()
  .core_init_env_check_kwargs(list(...))
  withr::local_options(list(digits = 22))
  .core_init_env_silence_warnings()

  arguments <- .core_init_env_clean_args(list(...))
  popped       <- .core_pop_arg(arguments, "start_fresh", FALSE)
  start_fresh  <- popped$value
  arguments    <- popped$args
  ssEnv        <- .core_init_env_bootstrap_session(result_folder, start_fresh)

  arguments <- .core_apply_defaults(arguments, .SS_DEFAULTS)
  if (!is.null(ssEnv$openai_api_key) && nzchar(ssEnv$openai_api_key))
    message("SEMseeker: set OPENAI_API_KEY in your environment to enable OpenAI features.")
  arguments <- .core_init_env_apply_computed_defaults(arguments)

  ssEnv     <- core_get_session_info()
  popped    <- .core_pop_arg(arguments, "dry_run", FALSE)
  dry_run   <- popped$value
  arguments <- popped$args
  if (dry_run) ssEnv$verbosity <- 4

  .core_log_info(" data will saved in this folder:", result_folder)
  ssEnv <- .core_init_env_setup_paths(ssEnv, result_folder)
  .core_init_env_setup_log_sink(ssEnv$session_folder)

  arguments <- util_keys_create(ssEnv, arguments)
  ssEnv     <- core_get_session_info()
  ssEnv$functionToExport <- .SS_FN_EXPORT
  .core_init_env_setup_progress(ssEnv$showprogress)

  arguments <- core_set_env_variable(arguments, "maxResources", maxResources)
  arguments <- core_set_env_variable(arguments, "parallel_strategy", "sequential")

  # Flush BEFORE the workers exist, not after. core_set_env_variable() no longer
  # writes to disk on every option (it used to, and that was most of the cost of
  # this function), so the session file is stale at this point; multisession
  # workers rebuild ssEnv by reading it, and core_parallel_session() is where
  # they are spawned. The final flush at the end of this function is too late
  # for them.
  core_update_session_info(core_get_session_info())

  core_parallel_session()
  ssEnv <- core_get_session_info()

  .core_init_env_log_focus(ssEnv)
  .core_init_env_validate_args(arguments)
  # one-line WARNING when R is linked against a single-thread BLAS.
  # Hot for the batch families (limma_/voom_) - solve()/crossprod()
  # inside lmFit scale ~linearly with cores on Accelerate/OpenBLAS/MKL.
  .core_warn_blas_single_thread()
  if (dry_run) .core_init_env_handle_dry_run(ssEnv)

  core_update_session_info(ssEnv)
  return(ssEnv)
}
