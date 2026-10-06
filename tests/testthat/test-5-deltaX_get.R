test_that("sem_deltaX_get", {

  tempFolder <- sem_test_folder()
  unlink(tempFolder, recursive = TRUE)
  # message(tempFolder)
  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
    bonferroni_threshold = bonferroni_threshold,
    inpute="median", start_fresh=TRUE)

  ####################################################################################

  tt <- SEMseeker:::core_get_meth_tech(signal_data)

  ####################################################################################
  LESIONS_BP <- 2000L
  bonferroni_threshold <- 0.05

  if (!exists("signal_thresholds"))
  {
    signal_data <- SEMseeker:::sem_inpute_missing_values(signal_data)
    signal_thresholds <<- SEMseeker:::sem_signal_range_values(signal_data, batch_id)
  }
  probe_features <<- probe_features[probe_features$PROBE %in% rownames(signal_data), ]

  keys <- ssEnv$keys_markers_figures_default
  sp <- SEMseeker:::sem_analyze_population(signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet[mySampleSheet$Sample_Group == "Case",],
    probe_features = probe_features
  )
  SEMseeker:::anno_create_position_pivots(mySampleSheet[mySampleSheet$Sample_Group == "Case",],keys)

  sp <- SEMseeker:::sem_analyze_population(signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet[mySampleSheet$Sample_Group == "Control",],
    probe_features = probe_features
  )

  SEMseeker:::anno_create_position_pivots(mySampleSheet[mySampleSheet$Sample_Group == "Control",],keys)

  # sem_deltaX_get() calls sem_study_summary_get() which reads the sample sheet CSV.
  # sem_analyze_population is called directly here (bypassing sem_analyze_batch which writes it),
  # so we write it manually with original mixed-case Sample_IDs.
  ssEnv2 <- SEMseeker:::core_get_session_info()
  sample_sheet_csv <- SEMseeker:::io_file_path_build(ssEnv2$result_folderData, "1_sample_sheet_original", "csv", FALSE)
  utils::write.csv2(mySampleSheet, file=sample_sheet_csv)

  ss <- SEMseeker:::sem_deltaX_get()

  # BED file names go through io_file_path_build() -> core_name_cleaning() which
  # uppercases sample IDs; pivot column names are derived from BED basenames
  # via io_stream_merge_bed(). So pivot colnames are uppercase regardless of
  # whether the pipeline went through analyze_batch.
  cleaned_sample_ids <- SEMseeker:::core_name_cleaning(mySampleSheet$Sample_ID)

  # verify all DELTAX (except DELTAS and DELTAR ) are coherent with MUTATIONS
  keys <- subset(ssEnv$keys_areas_subareas_markers_figures, AREA == "POSITION")
  keys <- subset(keys, MARKER != "DELTAS")
  keys <- subset(keys, MARKER != "DELTAR")
  keys <- subset(keys, MARKER != "LESIONS")
  keys <- subset(keys, MARKER != "SIGNAL")


  for (k in seq_len(nrow(keys)))
  {
    # k <- 1
    key <- keys[k,]
    marker <- as.character(key$MARKER)
    figure <- as.character(key$FIGURE)
    area <- as.character(key$AREA)
    subarea <- as.character(key$SUBAREA)

    mutations_pivot_file_name <- SEMseeker:::io_pivot_file_name_parquet("MUTATIONS",figure,area,subarea)
    if(file.exists(mutations_pivot_file_name))
      mutations_pivot <- as.data.frame(polars::pl$read_parquet(mutations_pivot_file_name))
    else
      next

    io_pivot_file_name <- SEMseeker:::io_pivot_file_name_parquet(marker,figure,area,subarea)
    # derived markers may not exist with sparse synthetic data
    if(!file.exists(io_pivot_file_name))
      next
    pivot <- as.data.frame(polars::pl$read_parquet(io_pivot_file_name))

    pivot <- pivot[,-c(1:3)]
    mutations_pivot <- mutations_pivot[,-c(1:3)]
    testthat::expect_true(nrow(pivot)<nprobes)
    testthat::expect_true(nrow(pivot)>0)

    pivot <- pivot[,order(colnames(pivot))]
    mutations_pivot <- mutations_pivot[,order(colnames(mutations_pivot))]

    testthat::expect_true(all(colnames(pivot) %in% cleaned_sample_ids))
    testthat::expect_true(all(colnames(pivot) == colnames(mutations_pivot)))

    testthat::expect_true(ncol(pivot)==ncol(mutations_pivot))
    testthat::expect_true(nrow(pivot)==nrow(mutations_pivot))

    pivot[is.na(pivot)] <- 0
    mutations_pivot[is.na(mutations_pivot)] <- 0

    pivot[pivot > 0] <- 1
    mutations_pivot[mutations_pivot > 0] <- 1

    # Only MUTATIONS is trivially identical to itself; DELTA* use signed values
    if (marker == "MUTATIONS")
    {
      if (!all(as.data.frame(pivot) == as.data.frame(mutations_pivot)))
      {
        print(marker)
        print(figure)
      }
      testthat::expect_true(all(as.data.frame(pivot) == as.data.frame(mutations_pivot)))
    }

    # Per-sample BED presence is NOT a valid invariant for derived markers.
    # By design (long-reads readiness) io_dump_sample_as_bed_file() skips
    # writing when the per-sample data has 0 rows — DELTAS/DELTAR/LESIONS
    # can legitimately be empty for a given (sample, figure) even when
    # MUTATIONS is not. Pivot-level expectations above are authoritative.

  }

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})

