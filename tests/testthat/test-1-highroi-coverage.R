## High-line-ROI coverage for previously 0%-covered files.
##
## Covered:
##   .core_bulk_model_memory_gate()   pure memory-budget decision (no deps)
##   io_plot_file_name()              plot-path assembly (session)
##   assoc_results_save()             write inference CSV to canonical path
##   util_exploratory_analysis()      end-to-end exploratory report (215 lines)
##

# ---------------------------------------------------------------------------
# .core_bulk_model_memory_gate  (pure)
# ---------------------------------------------------------------------------

test_that(".core_bulk_model_memory_gate returns a coherent decision + breakdown", {
  g <- SEMseeker:::.core_bulk_model_memory_gate(n_probes = 1000, n_samples = 10, n_coef = 2)

  expect_type(g, "list")
  expect_true(all(c("decision", "y_mat_GB", "lmfit_GB", "mono_peak_GB",
                    "chunk_y_GB", "chunk_peak_GB", "fit_object_GB",
                    "total_RAM_GB", "mem_frac", "available_GB") %in% names(g)))
  expect_true(g$decision %in% c("monolithic", "chunked", "abort"))
  # exact numeric identities
  expect_equal(g$y_mat_GB, 1000 * 10 * 8 / 1024^3)
  expect_equal(g$lmfit_GB, 1.5 * g$y_mat_GB)
  expect_equal(g$mono_peak_GB, g$y_mat_GB + g$lmfit_GB)
  expect_true(g$mem_frac > 0 && g$mem_frac <= 1)
  # a tiny matrix always fits (or falls back to monolithic when RAM unknown)
  expect_equal(g$decision, "monolithic")
})

test_that(".core_bulk_model_memory_gate escalates away from monolithic when huge", {
  h <- SEMseeker:::.core_bulk_model_memory_gate(n_probes = 1e8, n_samples = 1000, n_coef = 2)
  # only assert the escalation when system RAM was actually detected
  if (!is.na(h$total_RAM_GB))
    expect_true(h$decision %in% c("chunked", "abort"))
  expect_gt(h$mono_peak_GB, h$chunk_peak_GB)   # monolithic peak > single-chunk peak
})

# ---------------------------------------------------------------------------
# io_plot_file_name  (session)
# ---------------------------------------------------------------------------

test_that("io_plot_file_name builds a marker/figure plot path", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  detail <- list(
    family_test          = "gaussian",
    transformation_y     = "",
    independent_variable = "AGE",
    samples_sql_condition = NULL
  )
  key <- list(AREA = "GENE", SUBAREA = "TSS200", MARKER = "MUTATIONS", FIGURE = "HYPER")

  path <- SEMseeker:::io_plot_file_name(detail, folder = tf, key = key)

  expect_type(path, "character")
  # AI-255: DEPTH left the path with the concept it named.
  expect_false(grepl("DEPTH", path, ignore.case = TRUE))
  expect_match(path, "MUTATIONS", ignore.case = TRUE)
  expect_match(path, "\\.png$", ignore.case = TRUE)
})

# ---------------------------------------------------------------------------
# assoc_results_save  (session, writes CSV)
# ---------------------------------------------------------------------------

test_that("assoc_results_save writes a readable inference CSV at the canonical path", {
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
    samples_sql_condition = NULL
  )
  df <- data.frame(AREA = "GENE", SUBAREA = "TSS200", VALUE = 0.42)

  path <- SEMseeker:::assoc_results_save(df, detail, marker = "MUTATIONS")

  expect_type(path, "character")
  expect_true(file.exists(path))
  back <- utils::read.csv2(path, stringsAsFactors = FALSE)
  expect_equal(nrow(back), 1L)
  expect_true("AREA" %in% colnames(back))
})

# ---------------------------------------------------------------------------
# util_exploratory_analysis  (end-to-end on the bundled GEO fixture)
# ---------------------------------------------------------------------------

test_that("util_exploratory_analysis runs end-to-end and writes an exploratory report", {
  tf <- sem_test_folder()
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE); unlink(tf, recursive = TRUE) },
          add = TRUE)

  expect_no_error(
    SEMseeker:::util_exploratory_analysis(
      categorical_variables = c("Sample_Group"),
      numerical_variables   = c("Phenotest", "Covariates1"),
      sample_sheet          = mySampleSheet,
      signal_data           = signal_data,
      result_folder         = tf,
      parallel_strategy     = "sequential",
      start_fresh           = TRUE,
      showprogress          = FALSE,
      verbosity             = 1
    )
  )

  expect_true(dir.exists(file.path(tf, "Data", "Exploratory_0")))
})
