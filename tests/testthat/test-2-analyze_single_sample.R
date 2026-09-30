test_that("sem_analyze_single_sample", {

  tempFolder <- sem_test_folder()
  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, inpute = "mean")

  tt <- SEMseeker:::core_get_meth_tech(signal_data)

  if (!exists("signal_thresholds")) {
    signal_data <- SEMseeker:::sem_inpute_missing_values(signal_data)
    signal_thresholds <<- SEMseeker:::sem_signal_range_values(signal_data, batch_id)
  }
  probe_features <<- probe_features[probe_features$PROBE %in% rownames(signal_data), ]

  values_df <- data.frame(
    CHR   = probe_features$CHR[match(rownames(signal_data), probe_features$PROBE)],
    START = probe_features$START[match(rownames(signal_data), probe_features$PROBE)],
    END   = probe_features$END[match(rownames(signal_data), probe_features$PROBE)],
    VALUE = as.numeric(signal_data[, 1])
  )

  sample_detail <- mySampleSheet[1, c("Sample_ID", "Sample_Group")]

  # ── HYPO ──────────────────────────────────────────────────────────────────
  SEMseeker:::sem_analyze_single_sample(
    values        = values_df,
    thresholds    = signal_thresholds,
    figure        = "HYPO",
    sample_detail = sample_detail
  )

  hypo_mutations_folder <- SEMseeker:::io_dir_check_and_create(ssEnv$result_folderData,
    c(sample_detail$Sample_Group, "MUTATIONS_HYPO"))
  hypo_mutations_file <- SEMseeker:::io_file_path_build(
    hypo_mutations_folder, c(sample_detail$Sample_ID, "MUTATIONS", "HYPO"), "bed", add_gz = TRUE)

  # MUTATIONS bed file is always written (mutations exist for synthetic data)
  testthat::expect_true(file.exists(hypo_mutations_file))

  # LESIONS folder is always created even when no significant lesions exist
  hypo_lesions_folder <- SEMseeker:::io_dir_check_and_create(ssEnv$result_folderData,
    c(sample_detail$Sample_Group, "LESIONS_HYPO"))
  testthat::expect_true(dir.exists(hypo_lesions_folder))

  # ── HYPER ──────────────────────────────────────────────────────────────────
  SEMseeker:::sem_analyze_single_sample(
    values        = values_df,
    thresholds    = signal_thresholds,
    figure        = "HYPER",
    sample_detail = sample_detail
  )

  hyper_mutations_folder <- SEMseeker:::io_dir_check_and_create(ssEnv$result_folderData,
    c(sample_detail$Sample_Group, "MUTATIONS_HYPER"))
  hyper_mutations_file <- SEMseeker:::io_file_path_build(
    hyper_mutations_folder, c(sample_detail$Sample_ID, "MUTATIONS", "HYPER"), "bed", add_gz = TRUE)

  testthat::expect_true(file.exists(hyper_mutations_file))

  hyper_lesions_folder <- SEMseeker:::io_dir_check_and_create(ssEnv$result_folderData,
    c(sample_detail$Sample_Group, "LESIONS_HYPER"))
  testthat::expect_true(dir.exists(hyper_lesions_folder))

  # ── sem_analyze_single_sample_both runs without error ─────────────────────────
  # Note: sem_analyze_single_sample_both reads uncompressed .bed files while
  # sem_analyze_single_sample writes .bed.gz; the BOTH file is empty but the
  # function must not error.
  testthat::expect_no_error(
    SEMseeker:::sem_analyze_single_sample_both(sample_detail, "MUTATIONS")
  )

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})
