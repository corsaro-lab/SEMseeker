test_that(" SEMseeker:::sem_sample_group_check", {


  tempFolder <- sem_test_folder()
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, iqrTimes = iqrTimes, inpute="median")

  ####################################################################################

  SEMseeker:::core_get_meth_tech(signal_data)

  ####################################################################################

  # all fine
  testthat::expect_true(is.null( SEMseeker:::sem_sample_group_check(mySampleSheet, signal_data)))

  ####################################################################################

  #Sample_ID has NA
  mySampleSheet$Sample_ID[1] <- NA
  testthat::expect_true(!is.null( SEMseeker:::sem_sample_group_check(mySampleSheet, signal_data)))

  ####################################################################################

  #Sample_Group has NA
  mySampleSheet$Sample_Group[1] <- NA
  testthat::expect_true(!is.null( SEMseeker:::sem_sample_group_check(mySampleSheet, signal_data)))

  ####################################################################################

  #Lost Sample_Group Values
  mySampleSheet$Sample_Group <- NA
  testthat::expect_true(!is.null( SEMseeker:::sem_sample_group_check(mySampleSheet, signal_data)))

  ####################################################################################

  #Duplicated Samnple Rows
  # mySampleSheet <- rbind(mySampleSheet,mySampleSheet)
  # expect_error( SEMseeker:::sem_sample_group_check(mySampleSheet, signal_data), "I'm STOPPING HERE!")
  # testthat::expect_true(!is.null( SEMseeker:::sem_sample_group_check(mySampleSheet, signal_data)))

  ####################################################################################

  #Duplicated Samnple Rows
  # expect_error( SEMseeker:::sem_sample_group_check(mySampleSheet, signal_data), "I'm STOPPING HERE!")

  ####################################################################################

  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(computationCluster)
  SEMseeker:::core_close_env()
  unlink(tempFolder,recursive = TRUE)
})

