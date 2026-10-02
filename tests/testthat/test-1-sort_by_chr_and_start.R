test_that("anno_sort_by_chr_and_start", {


  tempFolder <- sem_test_folder()
  probe_features$ABSOLUTE <- paste(probe_features$CHR, probe_features$START, sep="_")

  #order not matching
  second <- SEMseeker:::anno_sort_by_chr_and_start( probe_features[order(probe_features$START),])

  testthat::expect_true( SEMseeker:::util_test_match_order( SEMseeker:::anno_sort_by_chr_and_start(probe_features)$ABSOLUTE,second$ABSOLUTE))

  ####################################################################################

  # SEMseeker:::core_close_env()
}
)
