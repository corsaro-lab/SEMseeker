test_that("SEMseeker:::util_test_match_order", {

  ####################################################################################
  #same order
  testthat::expect_true( SEMseeker:::util_test_match_order( probe_features$ABSOLUTE,probe_features$ABSOLUTE  ) )

  ####################################################################################
  #ordre not matching
  testthat::expect_true( !SEMseeker:::util_test_match_order( probe_features$ABSOLUTE,sort(probe_features$ABSOLUTE, decreasing = TRUE)))

  ####################################################################################
  #values not matching
  testthat::expect_true( !SEMseeker:::util_test_match_order( probe_features[-nrow(probe_features),"ABSOLUTE"],probe_features[-1, "ABSOLUTE"] ))

  ####################################################################################
  #one of two is NULL
  testthat::expect_false( SEMseeker:::util_test_match_order( probe_features$ABSOLUTE,NULL) )

}
)
