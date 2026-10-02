# AI-096 (2026-06-09): lazy passthrough for the resume path. The legacy
# behaviour materialised the SIGNAL pivot into a ~12 GB R data.frame even
# in resume mode (where io_signal_save would just early-return), then doubled
# that with an R-side match+sort to align signal_data with probe_features.
# Peak ~24 GB R-side + Polars residuo on 367k × 4013 → macOS silent jetsam
# kill on 64 GB Macs (v18, v21 confirmed). This refactor:
#
#   - In RESUME mode (POSITION pivot in cache): keep signal_data as a
#     Polars LazyFrame end-to-end. Extract schema (colnames, nrow) lazily.
#     Reference population thresholds use a lazy column subset materialised
#     ONLY for the ~10% reference samples → 1-2 GB peak instead of 12.
#     io_signal_save call is skipped entirely (would have early-returned).
#     sem_analyze_population_bulk reads the pivot lazily from disk and ignores
#     the signal_data arg.
#
#   - In FRESH mode (pivot must be built): keep R-side materialisation
#     for sem_inpute_missing_values (median needs row-wise access) and for the
#     io_signal_save fresh-path write. BUT skip the R-side probe_features
#     match+sort - io_signal_save's per-chr chunked sort is the canonical
#     sort gate (see `single-sort-gate-at-pivot-save` memory) so the input
#     row order doesn't matter.
#
# Sort gate policy: there is exactly ONE canonical sort in the pipeline,
# inside `io_signal_save()`, chunked per-chr by START → sink_parquet. Every
# downstream consumer reads from disk and trusts that order. No re-sort.

sem_analyze_batch <- function(signal_data, sample_sheet)
{

 
  ssEnv <- core_get_session_info()
  batch_id <- ssEnv$running_batch_id
  # AI-224: idempotent - sem_core() already normalised the sheet, but
  # sem_analyze_batch() is also called directly (tests, resume tooling).
  sample_sheet <- core_normalize_sample_ids(sample_sheet)$sample_sheet

  # AI-027: read via unified dispatcher. CASE 2 (streaming merge) lets
  # the SEM step pick up raw bed/bedgraph files when the SIGNAL
  # pivot has not been materialised yet.
  signal_pivot <- io_read_pivot("SIGNAL", io_signal_figure(), "POSITION", "WHOLE")
  resume_mode  <- !is.null(signal_pivot)

  if (resume_mode) {
    # ----------------------------------------------------------------
    # RESUME PATH - lazy passthrough
    # ----------------------------------------------------------------
    # AI-061+ (2026-06-09): extract schema + row count from the RAW
    # POSITION pivot (signal_pivot) BEFORE the anno_position_pivot_to_probe
    # join. Reason: signal_pivot is a direct scan_parquet LazyFrame, so
    # collect_schema() and select(pl$len())$collect() resolve from the
    # parquet footer (O(1) metadata read). After anno_position_pivot_to_probe
    # builds the join LazyFrame, those same calls FORCE the join to
    # execute - Polars cannot infer the post-join schema/count without
    # running the join, which on ewas-scale (367k × 4014 cols) allocates
    # ~12 GB of Rust heap and triggers jetsam silently (v18, v21, v25-v30).
    # Sample columns are the SAME pre/post join (only CHR/START/END are
    # dropped and PROBE is added), so taking them from the raw pivot is
    # equivalent and cheap.
    schema_cols_position <- names(signal_pivot$collect_schema())
    # polars 1.x R bindings: $to_data_frame() does NOT exist; coerce via
    # as.data.frame() instead. (Pattern was inherited from older polars 0.x
    # samples and silently broke when this code path was finally exercised.)
    n_probes    <- as.integer(
      as.data.frame(signal_pivot$select(polars::pl$len())$collect())$len[1]
    )
    sample_cols <- setdiff(schema_cols_position,
                           c("CHR", "START", "END", "PROBE"))

    # anno_position_pivot_to_probe returns a LazyFrame post AI-096; the PROBE
    # column is the probe identifier, sample columns follow.
    if ("CHR" %in% schema_cols_position) {
      signal_lazy <- anno_position_pivot_to_probe(signal_pivot)
    } else {
      signal_lazy <- signal_pivot
    }
    rm(signal_pivot)

    # If anno_position_pivot_to_probe (legacy or another caller) returned a
    # DataFrame instead of LazyFrame, coerce.
    if (inherits(signal_lazy, "polars_data_frame")) {
      signal_lazy <- signal_lazy$lazy()
    }

    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " working on batch:", batch_id, " of ", n_probes,
              " rows and ", length(sample_cols), " samples (resume mode, lazy).")

    # Tech detection: in resume mode we MUST trust ssEnv$tech (set at
    # core_init_env time or pre-declared). Auto-detection requires reading
    # probe IDs + a sample of values which would materialise rows we
    # don't otherwise need. If ssEnv$tech is missing, we can fall back
    # to a tiny lazy collect of first 10k rows + 1 column.
    if (is.null(ssEnv$tech) || !nzchar(ssEnv$tech)) {
      probe_ids <- as.character(
        as.data.frame(signal_lazy$select("PROBE")$head(10000L)$collect())$PROBE
      )
      tech_signal <- data.frame(
        PROBE = probe_ids,
        stringsAsFactors = FALSE
      )
      rownames(tech_signal) <- probe_ids
      ssEnv <- core_get_meth_tech(tech_signal)
      ssEnv <- core_get_session_info()
    } else {
      core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " technology pre-declared as '", ssEnv$tech,
                "'; using ssEnv$tech in resume mode.")
    }

    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " I will work on:", n_probes, " PROBES.")

    # Build probe_features by reading ONLY the PROBE column (tiny).
    core_log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"),
              " pre-probe_ids_collect mem_MB=", round(sum(gc()[, "(Mb)"]), 1))
    probe_ids_vec <- as.character(
      as.data.frame(signal_lazy$select("PROBE")$collect())$PROBE
    )
    core_log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"),
              " post-probe_ids_collect mem_MB=", round(sum(gc()[, "(Mb)"]), 1),
              " n_probe_ids=", length(probe_ids_vec))

    # AI-074: the gate applies to the resume path too - resuming from a SIGNAL
    # pivot still runs the SEM detection downstream.
    sem_coverage_gate(probe_ids_vec)

    if (ssEnv$tech %in% c("WGBS", "LONGREAD")) {
      probe_features <- io_coord_probe_features(probe_ids_vec)
    } else {
      probe_features <- anno_probe_features_get("PROBE")
      core_log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"),
                " post-anno_probe_features_get mem_MB=", round(sum(gc()[, "(Mb)"]), 1),
                " n_rows=", nrow(probe_features))
      probe_features <- probe_features[probe_features$PROBE %in% probe_ids_vec, ]
    }
    core_log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"),
              " post-probe_features_filter mem_MB=", round(sum(gc()[, "(Mb)"]), 1),
              " n_rows=", nrow(probe_features))
    # NO anno_sort_by_chr_and_start - sort gate is io_signal_save (already written).
    # NO signal_data row reorder - pivot rows are already canonical.

    # sem_sample_group_check expects something with colnames(signal_data) →
    # pass a zero-row placeholder with the right sample-column names.
    signal_data_check <- as.data.frame(matrix(numeric(0),
                                              nrow = 0, ncol = length(sample_cols)))
    colnames(signal_data_check) <- sample_cols
    sample_group_checkResult <- sem_sample_group_check(sample_sheet, signal_data_check)
    if (!is.null(sample_group_checkResult)) stop(sample_group_checkResult)
    rm(signal_data_check)

    # io_signal_save would early-return (POSITION pivot exists). Skip the
    # call entirely - saves a function frame and the misleading log line.
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " Signal data already saved (resume mode); skipping signal_save.")

    # Reference population thresholds: lazy select reference sample
    # columns only, then collect into a small R matrix (~10% of full
    # width). For ewas (4013 sample) reference ~ 400 sample → 1.5 GB
    # instead of 12 GB peak.
    referencePopulationSampleSheet <- sample_sheet[sample_sheet$Sample_Group == "Reference", ]
    ref_sample_ids <- intersect(sample_cols, referencePopulationSampleSheet$Sample_ID)
    if (length(ref_sample_ids) < 1L) {
      core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
                " No Reference sample columns found in SIGNAL pivot.")
      stop()
    }
    core_log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"),
              " pre-ref-collect mem_MB=", round(sum(gc()[, "(Mb)"]), 1),
              " n_ref=", length(ref_sample_ids))
    referencePopulationMatrix <- as.data.frame(
      signal_lazy$select(c("PROBE", ref_sample_ids))$collect()
    )
    rownames(referencePopulationMatrix) <- referencePopulationMatrix$PROBE
    referencePopulationMatrix$PROBE     <- NULL
    core_log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"),
              " post-ref-collect mem_MB=", round(sum(gc()[, "(Mb)"]), 1),
              " dim=", nrow(referencePopulationMatrix), "x", ncol(referencePopulationMatrix))

    populationControlRangeBetaValues <- as.data.frame(
      sem_signal_range_values(referencePopulationMatrix, batch_id, probe_features)
    )
    rm(referencePopulationMatrix); gc()
    core_log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"),
              " post-thresholds mem_MB=", round(sum(gc()[, "(Mb)"]), 1))

    # Sample-sheet cleanup (drop duplicates between Reference and other).
    referenceSamples <- sample_sheet[sample_sheet$Sample_Group == "Reference", ]
    otherSamples     <- sample_sheet[sample_sheet$Sample_Group != "Reference", ]
    referenceSamples <- referenceSamples[!(referenceSamples$Sample_ID %in% otherSamples$Sample_ID), ]
    sample_sheet     <- rbind(otherSamples, referenceSamples)
    core_log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"),
              " post-sample_sheet_rbind mem_MB=", round(sum(gc()[, "(Mb)"]), 1),
              " n_sample_sheet_rows=", nrow(sample_sheet))

    # bulk_population path: sem_analyze_population_bulk reads the pivot
    # lazily from disk and ignores the signal_data argument. Pass NULL
    # explicitly to avoid the caller assuming an R data.frame.
    if (isTRUE(ssEnv$bulk_population)) {
      core_log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"),
                " pre-apb_call mem_MB=", round(sum(gc()[, "(Mb)"]), 1))
      sem_analyze_population_bulk(
        signal_data       = NULL,
        sample_sheet      = sample_sheet,
        signal_thresholds = populationControlRangeBetaValues,
        probe_features    = probe_features
      )
      # AI-061+ (2026-06-09): release the thresholds R data.frame after
      # the bulk pass - sem_analyze_population_bulk copied the data into
      # polars (Rust heap) and rm()'d its own local binding, but R
      # would otherwise keep this parent-frame reference alive for the
      # rest of sem_analyze_batch's body (which we are about to return from
      # anyway, but make the intent explicit so future edits don't
      # accidentally rely on this binding).
      rm(populationControlRangeBetaValues)
      invisible(gc(verbose = FALSE))
      core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " Batch completed (bulk mode, lazy resume):", batch_id)
      return(invisible(NULL))
    }

    # Non-bulk path: per-population subset is computed via lazy select +
    # collect on demand. Each population sees its own materialisation,
    # bounded by population size.
    for (i in seq_along(ssEnv$keys_sample_groups[, 1])) {
      sample_group <- ssEnv$keys_sample_groups[i, 1]
      populationSampleSheet <- sample_sheet[sample_sheet$Sample_Group == sample_group, ]
      populationMatrixColumns <- intersect(sample_cols, populationSampleSheet$Sample_ID)

      if (length(populationMatrixColumns) == 0L) {
        core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                  " Population ", sample_group,
                  " is empty in resume mode; samples may be in another group.")
        next
      }
      core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " Working on population ", sample_group, " with ",
                n_probes, " probes (lazy resume).")
      population_signal <- as.data.frame(
        signal_lazy$select(c("PROBE", populationMatrixColumns))$collect()
      )
      rownames(population_signal) <- population_signal$PROBE
      population_signal$PROBE     <- NULL
      sem_analyze_population(
        signal_data       = population_signal,
        sample_sheet      = populationSampleSheet,
        signal_thresholds = populationControlRangeBetaValues,
        probe_features    = probe_features
      )
      rm(population_signal); gc()
    }

    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " Batch completed (lazy resume):", batch_id)
    return(invisible(NULL))
  }

  # ------------------------------------------------------------------
  # FRESH PATH - R-side materialisation required for inpute + io_signal_save
  # ------------------------------------------------------------------
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " working on batch:", batch_id, " of ", nrow(signal_data),
            " rows and ", ncol(signal_data), " samples (fresh mode).")
  # AI-224: normalise the signal columns with the SAME function used on the
  # sheet identifiers, so the name-based subsetting below cannot silently miss.
  signal_data <- core_normalize_sample_ids(sample_sheet = NULL,
                                           signal_data = signal_data)$signal_data

  # Transparent conversion: WGBS/LONGREAD coordinate input → synthetic probe IDs
  signal_data <- io_normalize_signal_input(signal_data)
  signal_data <- util_substitute_infinite(signal_data)
  signal_data <- sem_inpute_missing_values(signal_data)
  signal_data <- as.data.frame(signal_data)
  ssEnv <- core_get_meth_tech(signal_data)
  ssEnv <- core_get_session_info()

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " I will work on:", nrow(signal_data), " PROBES.")

  # AI-074: mandatory coverage gate - charts are produced on every run and the
  # run stops when the input barely overlaps the reference annotation.
  sem_coverage_gate(rownames(signal_data))

  # AI-106+ (2026-06-09): single source of truth for input → annotation
  # alignment. sem_prepare_batch_signal() centralises:
  #   - tech-specific probe_features build (manifest for Illumina,
  #     io_coord_probe_features for WGBS / LONGREAD)
  #   - dmr_annotation duplicate-PROBE collapse
  #   - intersection with input rownames
  #   - uniform sex-chromosome removal across all techs
  #   - signal_data ⇔ probe_features alignment with strict invariant
  #     nrow(signal_data) == nrow(probe_features) and matching row order
  # Replaces the ~50 lines of scattered annotation/filter/align logic
  # that previously lived here and was the source of the v35-v43 silent
  # drift between signal_data and probe_features.
  signal_data <- sem_prepare_batch_signal(
    signal_data           = signal_data,
    tech                  = ssEnv$tech,
    sex_chromosome_remove = isTRUE(ssEnv$sex_chromosome_remove)
  )
  probe_features <- attr(signal_data, "probe_features")
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " prepared batch signal: ", nrow(signal_data),
            " probes after sex-chr + manifest alignment (tech=",
            attr(signal_data, "tech"), ").")

  sample_group_checkResult <- sem_sample_group_check(sample_sheet, signal_data)
  if (!is.null(sample_group_checkResult)) {
    stop(sample_group_checkResult)
  }

  # AI-224: last gate before every name-based column subset below
  # (io_signal_save, reference matrix, per-population matrices). Fails with the
  # offending identifiers instead of R's opaque "undefined columns selected".
  core_normalize_sample_ids(sample_sheet, signal_data, require_all_ids = TRUE)

  io_signal_save(signal_data, sample_sheet, batch_id, probe_features = probe_features)
  core_log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"),
            " post-io_signal_save  mem_MB=", round(sum(gc()[, "(Mb)"]), 1))

  # Reference population subset and thresholds. probe IDs are kept as
  # rownames (preserved by data.frame column subsetting). NO extra
  # PROBE column wrapper - sem_signal_range_values reads probe IDs from
  # rownames.
  referencePopulationSampleSheet <- sample_sheet[sample_sheet$Sample_Group == "Reference", ]
  core_log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"),
            " pre-subset       mem_MB=", round(sum(gc()[, "(Mb)"]), 1),
            " n_ref=", nrow(referencePopulationSampleSheet))
  referencePopulationMatrix <- signal_data[, referencePopulationSampleSheet$Sample_ID, drop = FALSE]
  core_log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"),
            " post-subset      mem_MB=", round(sum(gc()[, "(Mb)"]), 1),
            " dim=", nrow(referencePopulationMatrix), "x", ncol(referencePopulationMatrix))

  if (plyr::empty(referencePopulationMatrix) || ncol(referencePopulationMatrix) < 1) {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
              " Empty signal_data ",
              format(Sys.time(), "%a %b %d %X %Y"))
    stop()
  }

  core_log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"),
            " pre-thresholds   mem_MB=", round(sum(gc()[, "(Mb)"]), 1))
  populationControlRangeBetaValues <- as.data.frame(
    sem_signal_range_values(referencePopulationMatrix, batch_id, probe_features))
  core_log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"),
            " post-thresholds  mem_MB=", round(sum(gc()[, "(Mb)"]), 1))
  rm(referencePopulationMatrix); gc()

  # Sample-sheet cleanup (drop duplicates between Reference and other).
  referenceSamples <- sample_sheet[sample_sheet$Sample_Group == "Reference", ]
  otherSamples     <- sample_sheet[sample_sheet$Sample_Group != "Reference", ]
  referenceSamples <- referenceSamples[!(referenceSamples$Sample_ID %in% otherSamples$Sample_ID), ]
  sample_sheet     <- rbind(otherSamples, referenceSamples)

  # AI-042: bulk_population path skips the per-sample loop entirely.
  if (isTRUE(ssEnv$bulk_population)) {
    sem_analyze_population_bulk(
      signal_data       = signal_data,
      sample_sheet      = sample_sheet,
      signal_thresholds = populationControlRangeBetaValues,
      probe_features    = probe_features
    )
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " Batch completed (bulk mode):", batch_id)
    return(invisible(NULL))
  }

  i <- 0
  for (i in seq_along(ssEnv$keys_sample_groups[, 1])) {
    sample_group <- ssEnv$keys_sample_groups[i, 1]
    populationSampleSheet <- sample_sheet[sample_sheet$Sample_Group == sample_group, ]
    populationMatrixColumns <- intersect(colnames(signal_data), populationSampleSheet$Sample_ID)

    if (length(populationMatrixColumns) == 0) {
      core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                " Population ", sample_group, " is empty.")
    } else {
      core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " Working on population ", sample_group, " with ",
                nrow(signal_data), " probes.")
      population_signal <- signal_data[, populationMatrixColumns, drop = FALSE]
      sem_analyze_population(
        signal_data       = population_signal,
        sample_sheet      = populationSampleSheet,
        signal_thresholds = populationControlRangeBetaValues,
        probe_features    = probe_features
      )
      rm(population_signal); gc()
    }
  }
  if (exists("signal_data", envir = environment(), inherits = FALSE)) {
    rm("signal_data", envir = environment())
  }
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " Batch completed:", batch_id)
}
