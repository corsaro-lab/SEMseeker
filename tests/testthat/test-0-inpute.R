test_that("inpute", {

  tempFolder <- sem_test_folder()

  # With the real GSE133774 fixture (ncol = 10), the removal threshold is
  # 0.1 * 10 = 1.  Any row with > 1 NA triggers the removal block, which then
  # also removes rows with exactly 1 NA (< 1 is FALSE).  We therefore split
  # the two behaviours into independent sub-tests so they don't interfere.

  # ── Sub-test A: imputation of partial NAs (no row removal) ────────────────
  # Inject exactly 1 NA per partial row (= 10% of ncol, NOT strictly > 10%,
  # so the removal block is never entered).
  local_imp <- signal_data[1:20, ]
  local_imp[1,  1] <- NA   # 1 NA  → imputed (1 > 1.0 = FALSE, no removal trigger)
  local_imp[5,  1] <- NA   # 1 NA  → imputed
  local_imp[10, 1] <- NA   # 1 NA  → imputed
  na_ante <- sum(is.na(local_imp))  # = 3

  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, iqrTimes = iqrTimes, inpute = "median")
  imp_median <- SEMseeker:::sem_inpute_missing_values(local_imp)
  testthat::expect_true(na_ante != 0)
  testthat::expect_true(sum(is.na(imp_median)) == 0)
  testthat::expect_true(nrow(imp_median) == nrow(local_imp))   # no rows removed

  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, iqrTimes = iqrTimes, inpute = "mean")
  imp_mean <- SEMseeker:::sem_inpute_missing_values(local_imp)
  testthat::expect_true(sum(is.na(imp_mean)) == 0)
  testthat::expect_true(nrow(imp_mean) == nrow(local_imp))

  # ── Sub-test B: removal of all-NA rows ────────────────────────────────────
  # Inject only all-NA rows (no partial NAs) so that the removal block removes
  # exactly the expected number of rows and does not also catch partial-NA rows.
  local_rm <- signal_data[1:20, ]
  local_rm[15, ] <- NA    # all-NA → removed
  local_rm[18, ] <- NA    # all-NA → removed
  nrow_missed <- sum(apply(local_rm, 1, function(x) all(is.na(x))))  # = 2

  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, iqrTimes = iqrTimes, inpute = "median")
  rm_median <- SEMseeker:::sem_inpute_missing_values(local_rm)
  testthat::expect_true(sum(is.na(rm_median)) == 0)
  testthat::expect_true(nrow(local_rm) - nrow(rm_median) == nrow_missed)

  ssEnv <- SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy, iqrTimes = iqrTimes, inpute = "mean")
  rm_mean <- SEMseeker:::sem_inpute_missing_values(local_rm)
  testthat::expect_true(sum(is.na(rm_mean)) == 0)
  testthat::expect_true(nrow(local_rm) - nrow(rm_mean) == nrow_missed)

  ####################################################################################

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})
