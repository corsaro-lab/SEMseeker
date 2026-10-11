test_that("annotations", {

  tempFolder <- sem_test_folder()
  # message(tempFolder)
  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
    bonferroni_threshold = bonferroni_threshold,
    inpute="median", start_fresh=TRUE)

  tt <- SEMseeker:::core_get_meth_tech(signal_data)

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

  # sem_deltaX_get needs sample sheet CSV (normally written by sem_analyze_batch)
  ssEnv2 <- SEMseeker:::core_get_session_info()
  sample_sheet_csv <- SEMseeker:::io_file_path_build(ssEnv2$result_folderData, "1_sample_sheet_original", "csv", FALSE)
  utils::write.csv2(mySampleSheet, file=sample_sheet_csv)

  SEMseeker:::sem_deltaX_get()

  SEMseeker:::anno_annotate_position_pivots()

  # reload keys some subareas may have been removed because not presente as positions
  keys <- subset(ssEnv$keys_areas_subareas_markers_figures, AREA != "POSITION")

  for (k in seq_len(nrow(keys)))
  {
    # k <- 1
    key <- keys[k,]
    marker <- key$MARKER
    figure <- key$FIGURE
    area <- key$AREA
    subarea <- key$SUBAREA

    mutations_pivot_file_name_parquet <- SEMseeker:::io_pivot_file_name_parquet("MUTATIONS",figure,area,subarea)
    if(file.exists(mutations_pivot_file_name_parquet))
      mutations_pivot <- polars::pl$read_parquet(mutations_pivot_file_name_parquet)|> as.data.frame()
    else
      next

    mutations_position_pivot_file_name_parquet <- SEMseeker:::io_pivot_file_name_parquet("MUTATIONS",figure,"POSITION","")
    if(file.exists(mutations_position_pivot_file_name_parquet))
      mutations_position_pivot <- polars::pl$read_parquet(mutations_position_pivot_file_name_parquet)|> as.data.frame()
    else
      next

    io_pivot_file_name_parquet <- SEMseeker:::io_pivot_file_name_parquet(marker,figure,area,subarea)
    testthat::expect_true(file.exists(io_pivot_file_name_parquet))
    if(file.exists(io_pivot_file_name_parquet))
      pivot <- polars::pl$read_parquet(io_pivot_file_name_parquet)|> as.data.frame()
    else
      pivot <- NULL


    pivot <- pivot[,-c(1:3)]
    mutations_pivot <- mutations_pivot[,-c(1:3)]
    testthat::expect_true(nrow(pivot)<nprobes)
    testthat::expect_true(nrow(pivot)>0)

    pivot <- pivot[,order(colnames(pivot))]
    mutations_pivot <- mutations_pivot[,order(colnames(mutations_pivot))]
    mutations_position_pivot <- mutations_position_pivot[,order(colnames(mutations_position_pivot))]

    testthat::expect_true(all(colnames(pivot) == colnames(mutations_position_pivot)))
    testthat::expect_true(all(colnames(mutations_position_pivot) == colnames(mutations_pivot)))

    testthat::expect_true(all(colnames(pivot) %in% (mySampleSheet$Sample_ID)))
    testthat::expect_true(all(colnames(pivot) == colnames(mutations_pivot)))

    testthat::expect_true(ncol(pivot)==ncol(mutations_pivot))
    testthat::expect_true(nrow(pivot)==nrow(mutations_pivot))

    pivot[is.na(pivot)] <- 0
    mutations_pivot[is.na(mutations_pivot)] <- 0

    pivot[pivot > 0] <- 1
    mutations_pivot[mutations_pivot > 0] <- 1

    testthat::expect_true(all(as.data.frame(pivot) == as.data.frame(mutations_pivot)))
  }

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})
