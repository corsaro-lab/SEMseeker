test_that("sem_signal_range_values", {

  tempFolder <- sem_test_folder()


  # test range calculation with missed values: sem_signal_range_values must error when data has NAs
  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, iqrTimes = iqrTimes)
  # anno_probe_features_get requires tech to be set; build equivalent manually
  probe_features_local <- probe_features[probe_features$PROBE %in% rownames(signal_data), c("CHR","START","END","PROBE")]
  probe_features_local <- probe_features_local[!is.na(probe_features_local$CHR),]
  # PROBE_WHOLE = canonical (AREA, SUBAREA) alias of PROBE — kept on probe_features
  # so sem_signal_range_values can drop it pre-write (and downstream association can
  # do `probe_features[[area_subarea]]` lookup uniformly).
  probe_features_local$PROBE_WHOLE <- probe_features_local$PROBE
  signal_data_with_na <- signal_data
  signal_data_with_na[1, 1] <- NA  # inject one NA to trigger the error
  testthat::expect_error(SEMseeker:::sem_signal_range_values(signal_data_with_na, batch_id, probe_features_local), "^ERROR:")

  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, iqrTimes = iqrTimes, inpute = "median")
  signal_data <- SEMseeker:::sem_inpute_missing_values(signal_data)

  ####################################################################################
  signal_thresholds <<- SEMseeker:::sem_signal_range_values(signal_data, batch_id, probe_features_local)
  testthat::expect_true(sum(colnames(signal_thresholds) %in% c("signal_inferior_thresholds","signal_superior_thresholds","signal_median_values","iqr","q1","q3"))==6)

  ####################################################################################
  # check thresholds file exists
  testthat::expect_true(file.exists(SEMseeker:::io_file_path_build(ssEnv$result_folderData ,c(batch_id, "signal_thresholds"),"parquet")))

  ####################################################################################
  # test no probe are lost
  testthat::expect_true(nrow(signal_thresholds)==nrow(signal_data))
  ####################################################################################
  # test lower limits are different from upper limits
  testthat::expect_true(signal_thresholds$signal_inferior_thresholds[1]!=signal_thresholds$signal_superior_thresholds[1])
  ####################################################################################
  # test median is different from lower limits
  testthat::expect_true(signal_thresholds$signal_inferior_thresholds[1]!=signal_thresholds$signal_median_values[1])
  ####################################################################################
  # test median is different from IQR
  testthat::expect_true(signal_thresholds$iqr[1]!=signal_thresholds$signal_median_values[1])
  ####################################################################################
  # test median is different from Q1
  testthat::expect_true(signal_thresholds$q1[1]!=signal_thresholds$signal_median_values[1])
  ####################################################################################
  # test median is different from Q3
  testthat::expect_true(signal_thresholds$q3[1]!=signal_thresholds$signal_median_values[1])
  ####################################################################################
  # test median is different from upper limits
  testthat::expect_true(sum(signal_thresholds$signal_inferior_thresholds<signal_thresholds$signal_superior_thresholds)/nrow(signal_thresholds)>0.9)

  ####################################################################################
  # test thresholds row are the same of probes of input matrix
  testthat::expect_true( nrow(signal_thresholds)==nprobes)

  ####################################################################################
  # test median is lower than Q3
  testthat::expect_true(sum(signal_thresholds$signal_median_values<signal_thresholds$q3)/nrow(signal_thresholds)>0.9)

  ####################################################################################
  # test median is higher than Q1
  testthat::expect_true(sum(signal_thresholds$signal_median_values>signal_thresholds$q1)/nrow(signal_thresholds)>0.9)

  SEMseeker:::core_close_env()
  unlink(tempFolder,recursive = TRUE)
})

