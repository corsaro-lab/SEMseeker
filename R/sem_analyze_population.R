#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param signal_data whole matrix of data to analyze.
#' @param sample_sheet name of samplesheet's column to use as control population
#' selector followed by selection value,
#' @param probe_features probe_features detail from 27 to EPIC illumina dataset
#' @param signal_thresholds thresholds defined to calculate epimutations
#' @return files into the result folder with pivot table and bedgraph.
#'   A BANNER is logged once per batch before the per-sample loop showing:
#'   input_positions, beta_range_positions, covered_by_inner_join - allows
#'   immediate audit of cross-run coverage (e.g. Nanopore sample vs Illumina reference).
#' @importFrom doRNG %dorng%
#'
sem_analyze_population <- function(signal_data, sample_sheet,signal_thresholds, probe_features) {

  ssEnv <- core_get_session_info()
  # #
  start_time <- Sys.time()
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " AnalyzePopulation warmingUP ")

  nrow_before <- nrow(signal_data)
  signal_data <- stats::na.omit(signal_data)
  nrow_after <- nrow(signal_data)
  if (nrow_before != nrow_after) {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Removed ", nrow_before - nrow_after, " rows with NA values")
    stop()
  }

  ### get signal_values ########################################################
  # idempotent normalisation of both sides - the per-sample subset
  # below is name-based and must not depend on the caller having cleaned them.
  .normalized <- core_normalize_sample_ids(sample_sheet, signal_data)
  sample_sheet <- .normalized$sample_sheet
  signal_data  <- .normalized$signal_data
  rm(.normalized)

  sample_sheet <- sample_sheet[order(sample_sheet[, "Sample_ID"], decreasing = FALSE), ]
  existent_samples <- colnames(signal_data)
  sample_names <- sample_sheet$Sample_ID
  missed_samples <- setdiff(setdiff(sample_names, existent_samples), "PROBE")

  if (length(missed_samples) != 0) {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " These samples data are missed: ", paste0(missed_samples, sep = " "))
  }

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " WarmedUP AnalyzePopulation")
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Start population analysis")

  # progress_bar <- progress::progress_bar$new(
  #   format = paste("INFO: Performing population analysis [:bar] :percent eta: :eta"),
  #   total = nrow(sample_sheet),
  #   clear = FALSE,
  #   width= 60)

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = seq_len(nrow(sample_sheet)))

  # precompute the set of bed/bedgraph files already on disk for each
  # (marker, figure) combo USED in this population, in ONE list.files() scan
  # per combo. Replaces the per-sample file.exists() (4-5 stat syscalls per
  # iteration x N samples = up to ~20k stat calls per population on big SCs)
  # with O(1) in-memory %in% lookups. Idempotent when folder is empty: set
  # stays character(0) -> %in% is always FALSE -> all samples get processed.
  .existing_bed_set <- function(marker, figure) {
    pattern <- paste0(marker, "_", figure)
    candidates <- list.files(
      ssEnv$result_folderData,
      pattern    = "\\.(bed|bedgraph)(\\.gz)?$",
      recursive  = TRUE,
      full.names = TRUE
    )
    candidates[grepl(paste0("/", pattern, "/"), candidates, fixed = TRUE)]
  }
  existing_signal_mean <- .existing_bed_set("SIGNAL", io_signal_figure())
  existing_mut_hyper   <- .existing_bed_set("MUTATIONS", "HYPER")
  existing_mut_hypo    <- .existing_bed_set("MUTATIONS", "HYPO")
  existing_deltas_hypo <- .existing_bed_set("DELTAS", "HYPO")
  existing_deltar_hypo <- .existing_bed_set("DELTAR", "HYPO")

  # pass skip_dir_create=TRUE to io_bed_file_name ONLY when the
  # existing_set for that (marker, figure) is non-empty (= prova certa che
  # la dir esiste). Empty set => skip_dir_create=FALSE => io_dir_check_and_create
  # corre normalmente per creare la dir. No assunzioni cieche.
  dir_known_signal_mean <- length(existing_signal_mean) > 0L
  dir_known_mut_hyper   <- length(existing_mut_hyper)   > 0L
  dir_known_mut_hypo    <- length(existing_mut_hypo)    > 0L
  dir_known_deltas_hypo <- length(existing_deltas_hypo) > 0L
  dir_known_deltar_hypo <- length(existing_deltar_hypo) > 0L

  variables_to_export <- c("sample_sheet", "signal_data", "sem_analyze_single_sample", "ssEnv",
    "signal_superior_thresholds","sem_deltar_single_sample","signal_inferior_thresholds","iqr","signal_median_values",
    "bt","bonferroni_threshold", "probe_features", "sem_analyze_single_sample_both", "sem_delta_single_sample", "progress_bar",
    "progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace","sem_signal_single_sample",
    "core_get_session_info","io_bed_file_name","signal_thresholds","core_update_session_info","anno_normalize_chr","io_signal_figure",
    "existing_signal_mean","existing_mut_hyper","existing_mut_hypo","existing_deltas_hypo","existing_deltar_hypo",
    "dir_known_signal_mean","dir_known_mut_hyper","dir_known_mut_hypo","dir_known_deltas_hypo","dir_known_deltar_hypo")
  i <- 1

  for(i in seq_len(nrow(sample_sheet))) {
  # foreach::foreach(i =1:nrow(sample_sheet), .export = variables_to_export) %dorng% {
    local_sample_detail <- sample_sheet[i,]
    ssEnv <- core_get_session_info()
    bed_filename <- io_bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "SIGNAL", io_signal_figure(), skip_dir_create = dir_known_signal_mean)
    if(!(bed_filename %in% existing_signal_mean)) {
      signal_values <- signal_data[,local_sample_detail$Sample_ID]
      sem_signal_single_sample( signal_values,local_sample_detail,probe_features)
    }
    if(ssEnv$showprogress)
      progress_bar(sprintf("Saving signal of sample: %s",local_sample_detail$Sample_ID))
  }
  gc()

  # ── Coverage banner - emitted ONCE before per-sample analysis ───────────────
  # Shows how many positions in the current run are covered by signal_thresholds.
  # Especially important for cross-run analysis (e.g. Nanopore sample vs an
  # Illumina reference batch passed via populationControlRangeBetaValues):
  # the user can immediately see if there is a partial or zero overlap.
  # Uses the first sample's signal bed file (just written above) for positions.
  {
    coverage_first_sample <- sample_sheet[1L, ]
    coverage_bed <- io_bed_file_name(coverage_first_sample$Sample_ID,
                                  coverage_first_sample$Sample_Group,
                                  "SIGNAL", io_signal_figure())
    if (file.exists(coverage_bed)) {
      coverage_sig <- utils::read.delim(coverage_bed, header = FALSE, sep = "\t")
      colnames(coverage_sig) <- c("CHR", "START", "END", "VALUE")
      coverage_sig$CHR <- anno_normalize_chr(coverage_sig$CHR, "internal")
      coverage_n_input  <- nrow(coverage_sig)
      coverage_n_ranges <- nrow(signal_thresholds)
      coverage_sig_lf <- polars::as_polars_df(
        coverage_sig[, c("CHR", "START", "END")])$
        with_columns(polars::pl$col("CHR")$cast(polars::pl$String))$lazy()
      coverage_thr_lf <- polars::as_polars_df(
        signal_thresholds[, c("CHR", "START", "END")])$
        with_columns(polars::pl$col("CHR")$cast(polars::pl$String))$lazy()
      coverage_n_covered <- coverage_sig_lf$join(
        coverage_thr_lf, on = c("CHR", "START", "END"), how = "inner"
      )$collect()$height
      core_log_event(
        "BANNER: ", format(Sys.time(), "%a %b %d %X %Y"),
        " [sem_analyze_population] Coverage -",
        " input_positions=", coverage_n_input,
        " | beta_range_positions=", coverage_n_ranges,
        " | covered_by_inner_join=", coverage_n_covered,
        " | analysis will run on ", coverage_n_covered, "/", coverage_n_input, " positions"
      )
      rm(coverage_sig, coverage_sig_lf, coverage_thr_lf,
         coverage_first_sample, coverage_bed,
         coverage_n_input, coverage_n_ranges, coverage_n_covered)
    }
  }

  progress_bar <- NULL
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = seq_len(nrow(sample_sheet)))

  rm(signal_data)
  # for(i in 1:nrow(sample_sheet)) {
  # How the internal helpers reach a worker. .packages attaches SEMseeker there,
  # and an attach exposes the EXPORTED functions only, so it is not what carries
  # them. .export is: it resolves each name in the frame this call is made from,
  # and that frame's enclosure is the namespace, so a name living there is found
  # and the function object travels with the task.
  #
  # Measured, because an earlier version of this comment asserted the opposite and
  # the calls below were qualified with the namespace on that basis - while the
  # same names were already in the list above. The prefix was doing nothing, and
  # R CMD check was reporting it.
  foreach::foreach(
    i = seq_len(nrow(sample_sheet)),
    .export = variables_to_export,
    .packages = "SEMseeker"
  ) %dorng% {
    # CRITICAL (E-14): multisession workers are fresh R processes where
    # .pkgglobalenv$ssEnv is empty. All internal helpers (io_bed_file_name,
    # sem_analyze_single_sample, etc.) call core_get_session_info() which reads from
    # .pkgglobalenv - NOT from the exported `ssEnv` variable. Without this
    # call, multisession workers fail with "core_get_session_info called without
    # result folder". See engineering-decisions.md §1.3.
    # in-memory only; saveRDS would happen N_samples × N_workers
    # times per SEM step otherwise (15 MB per write → catastrophic I/O).
    core_update_session_info(ssEnv, save_to_disk = FALSE)

    local_sample_detail <- sample_sheet[i,]

    # precomputed existing-bed sets => O(1) %in% lookups instead of
    # 4 stat syscalls per sample (~16k saved on a 4000-sample population).
    # skip_dir_create=TRUE because all destination dirs were ensured ONCE at
    # the top of sem_analyze_population - no per-sample dir_check_and_create.
    bed_mut_hyper   <- io_bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "MUTATIONS","HYPER", skip_dir_create = dir_known_mut_hyper)
    bed_mut_hypo    <- io_bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "MUTATIONS","HYPO", skip_dir_create = dir_known_mut_hypo)
    bed_deltas_hypo <- io_bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "DELTAS","HYPO", skip_dir_create = dir_known_deltas_hypo)
    bed_deltar_hypo <- io_bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "DELTAR","HYPO", skip_dir_create = dir_known_deltar_hypo)
    need_mut_hyper   <- !(bed_mut_hyper   %in% existing_mut_hyper)
    need_mut_hypo    <- !(bed_mut_hypo    %in% existing_mut_hypo)
    need_deltas_hypo <- !(bed_deltas_hypo %in% existing_deltas_hypo)
    need_deltar_hypo <- !(bed_deltar_hypo %in% existing_deltar_hypo)

    # Only pay the read.delim cost when at least one figure needs computing.
    if (need_mut_hyper || need_mut_hypo || need_deltas_hypo || need_deltar_hypo) {
      # Every call inside the worker has to be namespace-qualified: the foreach
      # expression is evaluated in an environment that does not have the package
      # namespace on its search path. io_signal_figure() was left unqualified
      # when it was introduced with the explicit aggregation, which resolves under an attached package but
      # not in a parallel worker.
      bed_filename <- io_bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "SIGNAL", io_signal_figure())
      signal_values <- utils::read.delim(bed_filename, header = FALSE, sep = "\t")
      colnames(signal_values) <- c("CHR", "START", "END", "VALUE")
      signal_values$CHR <- anno_normalize_chr(signal_values$CHR, "internal")

      if (need_mut_hyper)
        sem_analyze_single_sample( values = signal_values,thresholds = signal_thresholds, figure="HYPER", sample_detail = local_sample_detail)
      if (need_mut_hypo)
        sem_analyze_single_sample( values = signal_values,thresholds = signal_thresholds, figure="HYPO", sample_detail = local_sample_detail)
      if (need_deltas_hypo)
        sem_delta_single_sample( values = signal_values,thresholds = signal_thresholds , sample_detail = local_sample_detail)
      if (need_deltar_hypo)
        sem_deltar_single_sample ( values = signal_values, thresholds = signal_thresholds,sample_detail = local_sample_detail)
    }

    if(ssEnv$showprogress)
      progress_bar(sprintf("Performed sample: %s",local_sample_detail$Sample_ID))
  }

  gc()
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Row count result:", nrow(sample_sheet))
  if (exists("signal_data", envir = environment(), inherits = FALSE))
    rm("signal_data", envir = environment())

  # end-of-batch disk snapshot (workers used save_to_disk=FALSE
  # inside the foreach; here we persist the session exactly once).
  core_update_session_info(ssEnv, save_to_disk = TRUE)

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Completed population analysis ")
  end_time <- Sys.time()
  time_taken <- difftime(end_time,start_time, units = "mins")
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Completed population with summary - Time taken: ", time_taken, " minutes.")

}


