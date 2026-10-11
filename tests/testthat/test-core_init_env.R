test_that("test_init_env", {

  # run this test as first other wise the reuse of session is not testable
  tempFolder <- sem_test_folder()

  # test minimum number of parameters
  ssEnv <- SEMseeker:::core_init_env(tempFolder)
  testthat::expect_true(ssEnv$result_folderData == normalizePath(file.path(tempFolder, "Data"), mustWork = FALSE))
  SEMseeker:::core_close_env()

  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, areas = c("PROBE"), subareas= c(""), showprogress=TRUE, start_fresh=TRUE, maxResources = 10)
  # POSITION is always injected internally by util_keys_create() (needed for sem_study_summary_total).
  # So requesting areas=c("PROBE") yields keys containing PROBE + POSITION only.
  testthat::expect_true(nrow(subset(ssEnv$keys_areas_subareas_markers_figures, !AREA %in% c("PROBE", "POSITION")))==0)


  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, areas = c("GENE"), subareas= c("BODY"), showprogress=TRUE, start_fresh=TRUE, maxResources = 10)
  testthat::expect_true(ssEnv$maxResources == "10")
  SEMseeker:::core_close_env()

  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,areas = c("GENE","DMR"), showprogress=TRUE, start_fresh=TRUE)
  testthat::expect_true(length(unique(ssEnv$keys_areas_subareas[,"AREA"]))!=1)

  unlink(tempFolder,recursive = TRUE)
  assign("ssEnv", NULL, envir=.pkgglobalenv)

  ####################################################################################

  testthat::expect_error( SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, figures="HYPPO"), "Unrecognised argument")
  ####################################################################################

  testthat::expect_error( SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, markers="HYPPO"), "Unrecognised argument")

  ####################################################################################

  testthat::expect_error( SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, areas="HYPPO"), "Unrecognised argument")

  ####################################################################################

  # expect_error( SEMseeker:::core_init_env(tempFolder, parallel_strategy ="cluster"), "Unrecognised argument")

  ####################################################################################

  testthat::expect_error( SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, subareas="HYPPO"), "Unrecognised argument")


  testthat::expect_error( SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, test_it="HYPPO"), "Unrecognised argument")

  ####################################################################################

  unlink(tempFolder,recursive = TRUE)
  # SEMseeker:::core_close_env()
}
)
