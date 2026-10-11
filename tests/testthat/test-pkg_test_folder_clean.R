# Meta-test: no run writes inside tests/testthat/.
#
# A test that passes NA, "" or a relative path as result folder makes the
# package write its output into the working directory, which during a test run
# is this folder. Nothing fails when that happens: the run completes and the
# output stays behind. A folder named NA/ grew here that way to 42 files, from
# tests that had been handed a result folder that did not exist.
#
# The check lists what this folder may contain and fails on anything else,
# naming it. Test files run in parallel, so a write made later in the same run
# is caught on the next local run, where the folder is reused; under
# R CMD check it is caught when the writer happens to run first. Neither is
# complete, both are better than the silence there was.

test_that("tests/testthat contains only test sources", {
  allowed <- c(
    "^test-.*\\.R$",
    "^helper-.*\\.R$",
    "^setup.*\\.R$",
    "^teardown.*\\.R$",
    "^test_template\\.R$",
    "^_snaps$",
    "^_disabled_.*\\.txt$",
    "^testthat-problems\\.rds$"
  )
  entries <- list.files(testthat::test_path(), all.files = TRUE,
                        no.. = TRUE)
  ok <- Reduce(`|`, lapply(allowed, grepl, x = entries))
  expect_true(all(ok),
              info = paste("written into tests/testthat by a run:",
                           paste(entries[!ok], collapse = ", ")))
})
