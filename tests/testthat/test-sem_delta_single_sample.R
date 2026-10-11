test_that("sem_delta_single_sample",{

  tempFolder <- sem_test_folder()
  ssEnv <- SEMseeker:::core_init_env(tempFolder, inpute="median")

  ####################################################################################

  tt <- SEMseeker:::core_get_meth_tech(signal_data)

  probe_features <<- probe_features[probe_features$PROBE %in% rownames(signal_data), ]
  if (!exists("signal_thresholds"))
  {
    signal_data <- SEMseeker:::sem_inpute_missing_values(signal_data)
    signal_thresholds <- SEMseeker:::sem_signal_range_values(signal_data, batch_id,probe_features)
  }

  signal_data$PROBE <- rownames(signal_data)
  signal_data <- merge(signal_data, probe_features, by = "PROBE", all.x = TRUE)
  values <- cbind(signal_data[, c("PROBE","CHR", "START", "END")],signal_data[, colnames(signal_data)[2]])
  colnames(values) <- c("PROBE","CHR", "START", "END", "VALUE")
  values <- values[,c("CHR", "START", "END", "VALUE")]

  # Pick the first Control sample; mySampleSheet may have 16 rows (Reference reuse
  # pattern) so we cannot use row 1 directly — it may be a Reference sample.
  sample_detail <- mySampleSheet[mySampleSheet$Sample_Group == "Control", c("Sample_ID","Sample_Group")][1, , drop = FALSE]

  dss <- SEMseeker:::sem_delta_single_sample(
    values = values,
    thresholds = signal_thresholds,
    sample_detail = sample_detail
  )

  result_folderData  <-  SEMseeker:::io_dir_check_and_create(tempFolder, "Data")
  outputFolder <- SEMseeker:::io_dir_check_and_create(result_folderData, c(sample_detail$Sample_Group, "DELTAS_HYPER"))
  fileName <- SEMseeker:::io_file_path_build(outputFolder, c(sample_detail$Sample_ID, "DELTAS", "HYPER"), "bedgraph", add_gz = TRUE)
  testthat::expect_true(file.exists(fileName))

  # message("fileName: ", fileName)
  # test I can open it
  res <- read.table(gzfile(fileName), header = FALSE)
  # message("res: ", res)
  testthat::expect_true(nrow(res)>0)

  ####################################################################################
  # result_folderData  <-  SEMseeker:::io_dir_check_and_create(tempFolder, "Data")
  # outputFolder <- SEMseeker:::io_dir_check_and_create(result_folderData,c("Control","DELTAS_HYPERS"))
  # fileName <- SEMseeker:::io_file_path_build(outputFolder,c(mySampleSheet[1,c("Sample_ID")],"DELTAS","HYPERS"), "bedgraph", add_gz = TRUE)
  # testthat::expect_true(file.exists(fileName))
  #
  # # message("fileName: ", fileName)
  # # test I can open it
  # res <- read.table(gzfile(fileName), header = FALSE)
  # # message("res: ", res)
  # testthat::expect_true(nrow(res)> 0)

  ####################################################################################

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)

})
