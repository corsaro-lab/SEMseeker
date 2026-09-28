# Tests for the parameter array length-consistency check at the start of
# enrichment_analysis() (AI-035).
#
# enrichment_analysis iterates over 4 parallel arrays:
#   pvalue_columns, adjustment_methods, adjust_per_area_s, adjust_globally_s
# all must have the same length (each index i defines one tuple).
# Pre-AI-035 the error message said only "must have the same length";
# post-AI-035 it lists each array with its actual length and tells the
# user where to fix it.

# Minimal inference_details: enough to get past the first subset() filter of
# enrichment_analysis() and reach the length check. It used to say
# "with depth_analysis == 3"; that column is retired, and the row never carried
# it in the first place.
minimal_inf_details <- function() {
  data.frame(
    independent_variable = "Tumour_Stage_N",
    family_test          = "polynomial_2_1",
    covariates           = "Horvath",
    covariates_dummy     = "Tissue_Locus",
    transformation_y     = "none",
    transformation_x     = "none",
    filter_p_value       = FALSE,
    samples_sql_condition       = "",
    areas_sql_condition         = "",
    association_results_sql_condition = "",
    collinearity_check          = TRUE,
    covariates_pca              = TRUE,
    stringsAsFactors            = FALSE
  )
}

test_that("enrichment_analysis rejects mismatched array lengths", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]

  expect_error(
    SEMseeker::enrichment_analysis(
      inference_details   = minimal_inf_details(),
      adjust_per_area_s   = c(FALSE, FALSE, FALSE),
      adjust_globally_s   = c(FALSE, FALSE, FALSE),
      pvalue_columns      = c("PVALUE_A", "PVALUE_B", "PVALUE_C"),
      adjustment_methods  = c("NONE", "NONE"),  # length 2, mismatched
      alphas              = 0.05,
      study               = "test",
      significance        = TRUE,
      statistic_parameter = "",
      path_dbs            = c(),
      phenolyzer_folder_bin = "",
      disease             = "",
      result_folder       = tempFolder,
      maxResources        = 10,
      showprogress        = FALSE,
      verbosity           = 2
    ),
    "lengths must match"
  )
})

test_that("enrichment_analysis error lists each array's actual length", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]

  err <- tryCatch(
    SEMseeker::enrichment_analysis(
      inference_details   = minimal_inf_details(),
      adjust_per_area_s   = c(FALSE, FALSE, FALSE),
      adjust_globally_s   = c(FALSE, FALSE, FALSE),
      pvalue_columns      = c("PVALUE_A", "PVALUE_B", "PVALUE_C"),
      adjustment_methods  = c("NONE", "NONE"),
      alphas              = 0.05,
      study               = "test",
      significance        = TRUE,
      statistic_parameter = "",
      path_dbs            = c(),
      phenolyzer_folder_bin = "",
      disease             = "",
      result_folder       = tempFolder,
      maxResources        = 10,
      showprogress        = FALSE,
      verbosity           = 2
    ),
    error = identity
  )

  expect_s3_class(err, "error")
  msg <- conditionMessage(err)
  # Each of the 4 array names must appear with their length annotation.
  expect_match(msg, "pvalue_columns \\(length=3\\)")
  expect_match(msg, "adjustment_methods \\(length=2\\)")
  expect_match(msg, "adjust_per_area_s \\(length=3\\)")
  expect_match(msg, "adjust_globally_s \\(length=3\\)")
  # Hint pointing the user at the setup file.
  expect_match(msg, "setup file")
})

test_that("enrichment_analysis accepts matched array lengths (does not throw on length check)", {
  # We pass matched lengths and EXPECT the length check to pass. The function
  # will of course fail later (no real Inference CSV on disk) so we tolerate
  # any non-"length"-related error.
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]

  err <- tryCatch(
    SEMseeker::enrichment_analysis(
      inference_details   = minimal_inf_details(),
      adjust_per_area_s   = c(FALSE),
      adjust_globally_s   = c(FALSE),
      pvalue_columns      = c("PVALUE_A"),
      adjustment_methods  = c("NONE"),
      alphas              = 0.05,
      study               = "test",
      significance        = TRUE,
      statistic_parameter = "",
      path_dbs            = c(),
      phenolyzer_folder_bin = "",
      disease             = "",
      result_folder       = tempFolder,
      maxResources        = 10,
      showprogress        = FALSE,
      verbosity           = 2
    ),
    error = identity
  )

  if (inherits(err, "error")) {
    msg <- conditionMessage(err)
    expect_false(grepl("lengths must match", msg, fixed = TRUE))
  }
  # if no error, that's also fine
  succeed()
})
