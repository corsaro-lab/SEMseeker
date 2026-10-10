# =============================================================================
# The whole output of a fixed run is the same as the stored reference.
#
# The run, the fingerprint and the comparison are in helper-golden.R. A failure
# lists the artefacts and the cells that moved. If the change is intended,
# rewrite the reference with dev/golden-update.R and commit it with the change.
# =============================================================================

test_that("the output of the reference run matches the golden master", {
  skip_on_cran()
  reference_file <- file.path(.golden_dir(), "fingerprint.rds")
  expect_true(file.exists(reference_file))

  result_folder <- sem_test_folder()
  on.exit({
    try(SEMseeker:::core_close_env(), silent = TRUE)
    unlink(result_folder, recursive = TRUE)
  }, add = TRUE)

  .golden_run(result_folder)
  moved <- .golden_compare(.golden_read_reference(), .golden_fingerprint(result_folder))

  expect(length(moved) == 0L,
         paste0(length(moved), " difference(s) from the golden master:\n  ",
                paste(utils::head(moved, 40), collapse = "\n  "),
                if (length(moved) > 40) "\n  ..." else "",
                "\nIf the change is intended: dev/golden-update.R, and commit ",
                "tests/testthat/golden/ with it."))
})

test_that("the comparison catches a single moved cell, a lost row and a lost file", {
  ref <- list(
    "Inference/a.csv" = list(nrow = 2, columns = c("AREA", "PVALUE"),
                             values = data.frame(AREA = c("G1", "G2"), PVALUE = c(0.01, 0.2))),
    "Data/p.parquet"  = list(nrow = 10, columns = c("S1"),
                             summary = data.frame(column = "S1", n = 10, non_zero = 4, sum = 7,
                                                  sum_sq = 9, min = 0, max = 3)))
  expect_length(.golden_compare(ref, ref), 0L)

  cell <- ref; cell[["Inference/a.csv"]]$values$PVALUE[2] <- 0.3
  expect_match(.golden_compare(ref, cell), "column PVALUE row 2")

  row <- ref; row[["Inference/a.csv"]]$nrow <- 1
  expect_match(.golden_compare(ref, row), "1 rows, reference 2")

  big <- ref; big[["Data/p.parquet"]]$summary$sum <- 8
  expect_match(.golden_compare(ref, big), "column S1 summary moved")

  lost <- ref; lost[["Data/p.parquet"]] <- NULL
  expect_match(.golden_compare(ref, lost), "no longer written: Data/p.parquet")
})
