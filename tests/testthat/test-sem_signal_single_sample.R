test_that("signal",{

  tempFolder <- sem_test_folder()
  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, inpute="median")

  sample_detail <- mySampleSheet[1,]

  res <- SEMseeker:::sem_signal_single_sample(values = signal_data[,sample_detail$Sample_ID], sample_detail = sample_detail, probe_features = probe_features)

  # The FIGURE of SIGNAL is the scale of the session, resolved by the same
  # function the writer uses: the test follows BETA or MVALUE, it does not
  # decide it.
  figure <- SEMseeker:::io_signal_figure()

  folder_to_save <- SEMseeker:::io_dir_check_and_create(ssEnv$result_folderData,c(sample_detail$Sample_Group,paste0("SIGNAL","_", figure, sep = "")))
  signal_file <- SEMseeker:::io_file_path_build(folder_to_save,c(sample_detail$Sample_ID,"SIGNAL",figure),"bedgraph", add_gz = TRUE)
  testthat::expect_true(file.exists(signal_file))

  signal_file <- read.table(signal_file, sep="\t")
  testthat::expect_true(nrow(signal_file)==nrow(signal_data))


})
