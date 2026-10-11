test_that("signal-save",{

  tempFolder <- sem_test_folder()
  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, inpute="median")

  # In the normal pipeline sem_analyze_batch() calls core_get_meth_tech() before io_signal_save()
  # so that ssEnv$tech is already set when anno_probe_features_get() runs inside io_signal_save().
  SEMseeker:::core_get_meth_tech(signal_data)
  SEMseeker:::io_signal_save(signal_data,mySampleSheet,batch_id )

  # The FIGURE of SIGNAL is the scale of the session, resolved by the same
  # function the writer uses: the test follows BETA or MVALUE, it does not
  # decide it.
  figure <- SEMseeker:::io_signal_figure()

  # io_signal_save writes the probe-level parquet with subarea "WHOLE"
  signal_file <- SEMseeker:::io_pivot_file_name_parquet("SIGNAL", figure, "PROBE", "WHOLE")
  testthat::expect_true(file.exists(signal_file))
  # it also writes the position-level file
  position_file <- SEMseeker:::io_pivot_file_name_parquet("SIGNAL", figure, "POSITION", "WHOLE")
  testthat::expect_true(file.exists(position_file))


})
