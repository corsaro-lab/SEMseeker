## Coverage for pure, session-free IO/enrichment reader helpers.
##
## Covered:
##   enrich_results_get()          read + filter a pathway result CSV
##   io_source_data_get()          resolve a data.frame or a file path to a df
##   util_finalize_job_results()   final per-job save (empty-results guard path)

# ---------------------------------------------------------------------------
# enrich_results_get
# ---------------------------------------------------------------------------

test_that("enrich_results_get returns an empty frame for a missing file", {
  expect_equal(nrow(SEMseeker:::enrich_results_get(tempfile(fileext = ".csv"))), 0L)
})

test_that("enrich_results_get reads, significance-filters and top-N-truncates", {
  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  utils::write.csv2(
    data.frame(PATHWAY = c("p1", "p2", "p3"), FDR = c(0.01, 0.20, 0.001)),
    csv, row.names = FALSE)

  expect_equal(nrow(SEMseeker:::enrich_results_get(csv)), 3L)

  sig <- SEMseeker:::enrich_results_get(csv, significance_column = "FDR", alpha = 0.05)
  expect_equal(nrow(sig), 2L)                       # 0.01 and 0.001 pass

  top1 <- SEMseeker:::enrich_results_get(csv, significance_column = "FDR",
                                         alpha = 0.05, top = 1)
  expect_equal(nrow(top1), 1L)
  expect_equal(top1$PATHWAY, "p3")                  # smallest FDR first

  # required column missing -> empty
  expect_equal(nrow(SEMseeker:::enrich_results_get(csv, required_columns = "MISSING")), 0L)
})

# ---------------------------------------------------------------------------
# io_source_data_get
# ---------------------------------------------------------------------------

test_that("io_source_data_get returns a data.frame unchanged from an in-memory object", {
  res <- SEMseeker:::io_source_data_get(data.frame(a = 1:3, b = 4:6))
  expect_s3_class(res, "data.frame")
  expect_equal(res$a, 1:3)
})

test_that("io_source_data_get returns NULL for a non-existent path", {
  expect_message(
    expect_null(SEMseeker:::io_source_data_get(file.path(tempdir(), "no_such_file.csv"))),
    "not found")
})

test_that("io_source_data_get reads a CSV and promotes PROBE to rownames", {
  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  utils::write.csv2(
    data.frame(PROBE = c("cg1", "cg2"), S1 = c(0.1, 0.2), S2 = c(0.3, 0.4)),
    csv, row.names = FALSE)

  res <- SEMseeker:::io_source_data_get(csv)
  expect_s3_class(res, "data.frame")
  expect_equal(rownames(res), c("cg1", "cg2"))
  expect_false("PROBE" %in% colnames(res))
})

test_that("io_source_data_get with check_is_numeric stops on non-numeric columns", {
  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  utils::write.csv2(
    data.frame(PROBE = c("cg1", "cg2"), LABEL = c("x", "y")),
    csv, row.names = FALSE)

  expect_error(SEMseeker:::io_source_data_get(csv, check_is_numeric = TRUE))
})

# ---------------------------------------------------------------------------
# util_finalize_job_results  (empty-results guard path)
# ---------------------------------------------------------------------------

test_that("util_finalize_job_results is a no-op (NULL) when there are no results", {
  expect_null(SEMseeker:::util_finalize_job_results(
    results         = data.frame(),
    inference_detail = list(transformation_x = ""),
    family_test     = "gaussian",
    filter_p_value  = FALSE,
    fileNameResults = tempfile(fileext = ".csv"),
    start_time      = Sys.time(),
    processed_items = 0L))
})
