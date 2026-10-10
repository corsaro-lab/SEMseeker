## Coverage for association/enrichment ORCHESTRATORS via real sessions.
##
## The happy path of these functions requires full-pipeline artifacts
## (pivots + inference CSVs) whose generation is (a) tied to parallel_strategy
## in a way that is flaky under devtools::load_all and (b) affected by a known
## depth=3 regression that empties results (see test-7). These tests therefore
## exercise the deterministic ENTRY + guard/degenerate paths against a live
## session: argument handling, inference-path construction, missing-file and
## empty-input early returns. Full happy-path behaviour stays covered by the
## end-to-end tests (test-6/7/8-*).
##

# ---------------------------------------------------------------------------
# assoc_results_get — missing inference file -> empty data.frame
# ---------------------------------------------------------------------------

test_that("assoc_results_get returns an empty frame when the inference file is absent", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  detail <- list(
    covariates            = "",
    covariates_dummy      = "",
    covariates_pca        = FALSE,
    family_test           = "gaussian",
    transformation_y      = "",
    independent_variable  = "AGE",
    transformation_x      = "",
    samples_sql_condition = NULL,
    areas_sql_condition   = ""
  )

  res <- SEMseeker:::assoc_results_get(inference_detail = detail, marker = "MUTATIONS",
                                       area = "GENE", scope = "INSTANCE")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0L)
})

test_that("assoc_results_get refuses to guess which region class and which extent", {
  # `area` used to default to "GENE" and `scope` to NULL, meaning no
  # filter at all: a caller that said neither read the genes of every extent,
  # collapsed rows included, and had nothing to tell it so. Both are required,
  # and the refusal names what it wants.
  detail <- list(
    covariates            = "",
    covariates_dummy      = "",
    covariates_pca        = FALSE,
    family_test           = "gaussian",
    transformation_y      = "",
    independent_variable  = "AGE",
    transformation_x      = "",
    samples_sql_condition = NULL,
    areas_sql_condition   = ""
  )

  expect_error(
    SEMseeker:::assoc_results_get(inference_detail = detail, marker = "MUTATIONS",
                                  scope = "INSTANCE"),
    "'area' is required")
  expect_error(
    SEMseeker:::assoc_results_get(inference_detail = detail, marker = "MUTATIONS",
                                  area = "GENE"),
    "'scope' is required")
})

# ---------------------------------------------------------------------------
# assoc_data_extractor — no inference rows -> empty result
# ---------------------------------------------------------------------------

test_that("assoc_data_extractor returns an empty frame for zero inference rows", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  res <- SEMseeker:::assoc_data_extractor(
    inference_details = data.frame(),
    destination_folder = "",
    result_folder = "")             # "" -> reuse the active session

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0L)
})

# ---------------------------------------------------------------------------
# sem_available_metrics — no inference rows -> NULL
# ---------------------------------------------------------------------------

test_that("sem_available_metrics is a no-op (NULL) when there are no inference rows", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  res <- SEMseeker:::sem_available_metrics(
    inference_details = data.frame(),
    result_folder     = tf)

  expect_null(res)
})
