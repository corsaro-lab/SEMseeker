# =============================================================================
# At scope SAMPLE every sample of the request reaches the model.
#
# A SAMPLE pivot carries one row per aggregation and one column per sample,
# with AREA as its LAST column. The association loop used to take the FIRST
# column as the area name and drop it: the first sample vanished, came back
# from the merge with the sample sheet as NA, and the NA was turned into 0.
# Nothing failed, the statistic was simply computed on a sample whose burden
# had been set to zero.
#
# The expected statistic is computed here by hand, from the pivot read by
# column NAME and from the sample sheet the run declared, with cor.test().
# It is not obtained by running the package's own model code again.
# =============================================================================

test_that("scope SAMPLE: the statistic uses every sample, the first one included", {
  tempFolder <- sem_test_folder()
  on.exit({
    try(SEMseeker:::core_close_env(), silent = TRUE)
    unlink(tempFolder, recursive = TRUE)
  }, add = TRUE)

  samples <- mySampleSheet
  samples$Phenotest <- seq_len(nrow(samples)) * 1.5 + (seq_len(nrow(samples)) %% 3)

  SEMseeker::semseeker(
    input             = signal_data,
    sample_sheet      = samples,
    result_folder     = tempFolder,
    parallel_strategy = "sequential",
    areas             = "POSITION",
    markers           = "LESIONS",
    start_fresh       = TRUE,
    inpute            = "median",
    showprogress      = showprogress,
    verbosity         = verbosity
  )

  inference_details <- data.frame(
    independent_variable = "Phenotest",
    family_test          = "spearman",
    transformation_y     = "",
    transformation_x     = "",
    aggregation          = "SUM",
    scope                = "SAMPLE",
    filter_p_value       = FALSE,
    stringsAsFactors     = FALSE
  )
  SEMseeker::association_analysis(
    inference_details = inference_details,
    result_folder     = tempFolder,
    parallel_strategy = "sequential",
    markers           = "LESIONS",
    figures           = c("HYPER", "HYPO"),
    areas             = "POSITION",
    multiple_test_adj = "BH",
    showprogress      = showprogress,
    verbosity         = verbosity
  )

  data_folder <- file.path(tempFolder, "Data")
  sheet <- utils::read.csv2(file.path(data_folder, "SAMPLE_SHEET_RESULT.csv"))
  result_file <- list.files(file.path(tempFolder, "Inference"),
                            pattern = "^LESIONS_PHENOTEST_SPEARMAN\\.csv$",
                            recursive = TRUE, full.names = TRUE)
  expect_length(result_file, 1L)
  results <- utils::read.csv2(result_file)

  first_sample_counts <- FALSE
  for (figure in c("HYPER", "HYPO")) {
    pivot_file <- list.files(file.path(data_folder, "Pivots", "LESIONS"),
                             pattern = paste0("^LESIONS_", figure,
                                              "_SAMPLE_PROBE_WHOLE_SUM_.*\\.parquet$"),
                             full.names = TRUE)
    expect_length(pivot_file, 1L)
    pivot <- as.data.frame(polars::pl$read_parquet(pivot_file))

    # The case the defect needs: a first sample column that is not zero. If no
    # figure has it, this test could not tell the defect from the fix.
    first_sample <- setdiff(names(pivot), "AREA")[1]
    first_sample_counts <- first_sample_counts || pivot[[first_sample]][1] != 0

    expect_true(all(sheet$Sample_ID %in% names(pivot)))
    burden <- vapply(sheet$Sample_ID, function(id) as.numeric(pivot[[id]][1]),
                     numeric(1))
    expected <- suppressWarnings(
      stats::cor.test(sheet$Phenotest, burden, method = "spearman"))

    row <- results[results$FIGURE == figure, ]
    expect_equal(nrow(row), 1L, info = figure)
    expect_equal(row$RHO,    unname(expected$estimate), tolerance = 1e-9, info = figure)
    expect_equal(row$PVALUE, expected$p.value,          tolerance = 1e-9, info = figure)
  }
  expect_true(first_sample_counts)
})
