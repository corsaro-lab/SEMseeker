# Tests for assoc_validate_inference_schema()
#
# AI-035: strict schema validation of inference_details with fuzzy-match
# suggestions on typos. These tests pin the user-facing diagnostic so
# regressions on the error message are caught.

minimal_inference <- function() {
  data.frame(
    independent_variable = "Tumour_Stage_N",
    family_test          = "polynomial_2_1",
    covariates           = "Horvath",
    covariates_dummy     = "Tissue_Locus",
    transformation_y     = "none",
    filter_p_value       = FALSE,
    samples_sql_condition = "",
    collinearity_check   = TRUE,
    covariates_pca       = TRUE,
    stringsAsFactors     = FALSE
  )
}

test_that("assoc_validate_inference_schema accepts a minimal valid inference_details", {
  df <- minimal_inference()
  res <- expect_no_error(SEMseeker:::assoc_validate_inference_schema(df))
  # all 13 expected columns must be present after validation
  expected <- c("independent_variable", "family_test", "covariates",
                "covariates_dummy", "covariates_pca", "collinearity_check",
                "transformation_y", "transformation_x",
                "filter_p_value", "samples_sql_condition",
                "areas_sql_condition", "association_results_sql_condition")
  expect_true(all(expected %in% colnames(res)))
})

test_that("assoc_validate_inference_schema fills missing optional columns with NA", {
  df <- minimal_inference()
  # df has no transformation_x or areas_sql_condition
  expect_false("transformation_x" %in% colnames(df))
  expect_false("areas_sql_condition" %in% colnames(df))
  res <- SEMseeker:::assoc_validate_inference_schema(df)
  expect_true("transformation_x" %in% colnames(res))
  expect_true("areas_sql_condition" %in% colnames(res))
  expect_true(all(is.na(res$transformation_x)))
  expect_true(all(is.na(res$areas_sql_condition)))
})

test_that("assoc_validate_inference_schema errors on unknown columns in strict mode", {
  df <- minimal_inference()
  df$totally_unknown <- "x"
  expect_error(
    SEMseeker:::assoc_validate_inference_schema(df),
    "unknown column"
  )
})

test_that("assoc_validate_inference_schema fuzzy-suggests close-match typos", {
  df <- minimal_inference()
  df$samples_sql_condtion <- ""  # missing 'i' in 'condition'
  # error should mention the typo AND suggest the correct name
  expect_error(
    SEMseeker:::assoc_validate_inference_schema(df),
    "samples_sql_condtion.*samples_sql_condition"
  )
})

test_that("assoc_validate_inference_schema suggests covariates_dummy for covariate_dummy typo", {
  df <- minimal_inference()
  df$covariate_dummy <- "Tissue"  # missing 's' (-> close to covariates_dummy)
  err <- tryCatch(SEMseeker:::assoc_validate_inference_schema(df), error = identity)
  expect_s3_class(err, "error")
  # the error message should suggest the close match
  expect_match(conditionMessage(err), "covariate_dummy")
  expect_match(conditionMessage(err), "did you mean")
})

test_that("assoc_validate_inference_schema in lenient mode only warns, doesn't stop", {
  df <- minimal_inference()
  df$totally_unknown <- "x"
  expect_warning(
    res <- SEMseeker:::assoc_validate_inference_schema(df, strict = FALSE),
    "unknown column"
  )
  # in lenient mode the validation still fills missing required cols
  expect_true("transformation_x" %in% colnames(res))
})

test_that("assoc_validate_inference_schema whitelists engine-stamped runtime columns", {
  df <- minimal_inference()
  # These are stamped by the engine on a previously-run inference. They
  # should NOT trigger the unknown-column error when validating again
  # (e.g. for a resume scenario where validation runs on a CSV-roundtripped
  # inference_details).
  df$node_name        <- "host.local"
  df$session_id       <- 1L
  df$start_time       <- Sys.time()
  df$end_time         <- Sys.time()
  df$processed_items  <- 0L
  df$processed_time   <- 0
  expect_no_error(SEMseeker:::assoc_validate_inference_schema(df))
})
