test_that("set-env-variable", {

  # run this test as first other wise the reuse of session is not testable
  tempFolder <- sem_test_folder()

  # sex_chromosome_remove is FALSE by default
  ssEnv <- SEMseeker:::core_init_env(result_folder= tempFolder, parallel_strategy = parallel_strategy, areas = c("GENE"), subareas= c("BODY"), start_fresh=TRUE)
  testthat::expect_false(ssEnv$sex_chromosome_remove)

  # create list arguments
  arguments <- list(areas = c("GENE"), subareas = c("BODY"), sex_chromosome_remove = TRUE)
  arguments  <- SEMseeker:::core_set_env_variable(arguments =  arguments, var_name =  "sex_chromosome_remove",var_default_value =  TRUE)
  ssEnv  <- SEMseeker:::core_get_session_info(tempFolder)
  testthat::expect_true(ssEnv$sex_chromosome_remove)
  # check sex_chromosome_remove is not in the list
  testthat::expect_true(is.null(arguments[["sex_chromosome_remove"]]))

  SEMseeker:::core_close_env()

  # check reopening not fresh is FALSE
  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, areas = c("GENE"), subareas= c("BODY"), start_fresh=FALSE)
  # FALSE, which is what the comment above says and what the disk holds. The TRUE
  # set through core_set_env_variable() a few lines up lives in memory only: that
  # function defers the write and the caller is the one that flushes, and here
  # nobody did. This read used to be answered by the session the close had left
  # standing, which is what made it look like a round-trip through the file.
  #
  # Provisional, deliberately: if the session manager stops taking save_to_disk
  # from its caller and flushes on close instead, the TRUE does reach disk and the
  # assertion here goes back to expect_true.
  testthat::expect_false(ssEnv$sex_chromosome_remove)
  SEMseeker:::core_close_env()

  # change sex_chromosome_remove to TRUE
  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, areas = c("GENE"), subareas= c("BODY"), sex_chromosome_remove=TRUE, start_fresh=FALSE)
  ssEnv  <- SEMseeker:::core_get_session_info(tempFolder)
  arguments  <- SEMseeker:::core_set_env_variable(arguments =  arguments, var_name =  "sex_chromosome_remove",var_default_value =  TRUE)
  testthat::expect_true(is.null(arguments[["sex_chromosome_remove"]]))
  ssEnv  <- SEMseeker:::core_get_session_info(tempFolder)
  testthat::expect_true(ssEnv$sex_chromosome_remove=="TRUE")

  SEMseeker:::core_close_env()

})
