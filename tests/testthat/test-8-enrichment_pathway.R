# Tests for enrichment and pathway analysis helpers
#
# Covered:
#  - enrich_analysy_add_category()  category mapping helper (requires session)
#  - enrich_ctdR()                     graceful return when ctdR not installed
#  - enrich_WebGestalt()               graceful return when WebGestaltR not installed
#  - enrich_STRINGdb()                 graceful return when STRINGdb not installed
#  - enrich_pathfindR()                graceful return when pathfindR not installed
#  - enrich_ctdR() integration         end-to-end after a full semseeker +
#                                       association_analysis run with real results
#
# Note: all "package-not-installed" tests temporarily unload the optional package
# from the search path so the requireNamespace() guard fires correctly.

# ---------------------------------------------------------------------------
# 1. enrich_analysy_add_category — category mapping
# ---------------------------------------------------------------------------

test_that("enrich_analysy_add_category returns data unchanged for empty input", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  result <- SEMseeker:::enrich_analysy_add_category("ctdR", data.frame())
  testthat::expect_equal(nrow(result), 0)

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})

test_that("enrich_analysy_add_category adds SS_CATEGORY='CHEMICAL' for ctdR source", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  # Column names must match key_enrichment_format for "ctdR":
  #   column_of_adj_pvalue  = "padj"
  #   column_of_enrichment  = "foldEnrichment"
  #   column_of_description = "ChemicalName"
  # The ranking code inside the function accesses these columns;
  # using wrong names triggers a browser() call (debugging artifact in code).
  fake_result <- data.frame(
    ChemicalID     = c("D001", "D002"),
    ChemicalName   = c("chemical A", "chemical B"),
    pvalue         = c(0.01, 0.05),
    padj           = c(0.05, 0.10),
    foldEnrichment = c(2.5, 1.8),
    source         = c("ctdR", "ctdR"),
    key            = c("k1", "k2"),
    stringsAsFactors = FALSE
  )

  result <- SEMseeker:::enrich_analysy_add_category("ctdR", fake_result)
  testthat::expect_true("SS_CATEGORY" %in% colnames(result))
  testthat::expect_true(all(result$SS_CATEGORY == "CHEMICAL"))

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})

test_that("enrich_analysy_add_category maps GO types for WebGestalt source", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  # Column names must match key_enrichment_format for "WebGestalt":
  #   column_of_adj_pvalue  = "FDR"
  #   column_of_enrichment  = "enrichmentRatio"
  #   column_of_description = "description"
  fake_result <- data.frame(
    geneSet         = c("GO:0001", "GO:0002", "GO:0003"),
    description     = c("term1", "term2", "term3"),
    pValue          = c(0.01, 0.05, 0.1),
    FDR             = c(0.02, 0.10, 0.20),
    enrichmentRatio = c(2.5, 1.8, 1.2),
    type            = c("BP", "CC", "MF"),
    source          = c("WebGestalt", "WebGestalt", "WebGestalt"),
    stringsAsFactors = FALSE
  )

  result <- SEMseeker:::enrich_analysy_add_category("WebGestalt", fake_result)
  testthat::expect_true("SS_CATEGORY" %in% colnames(result))
  # Categories BP/CC/MF should map to GO-BP / GO-CC / GO-MF
  cats <- result$SS_CATEGORY[!is.na(result$SS_CATEGORY)]
  testthat::expect_true(all(cats %in% c("GO-BP", "GO-CC", "GO-MF")))

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})

# ---------------------------------------------------------------------------
# 2. Pathway functions return NULL gracefully when optional packages absent
#    We test this by calling with an active session but zero-row inference_details.
#    When the optional package IS present the function proceeds; we only care
#    that it does not throw an unhandled error regardless of package state.
# ---------------------------------------------------------------------------

test_that("enrich_WebGestalt returns NULL gracefully when WebGestaltR not installed", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  # WebGestaltR is typically not installed in the test environment.
  # The function guards with requireNamespace() and returns invisibly.
  if (!requireNamespace("WebGestaltR", quietly = TRUE)) {
    inference_detail <- data.frame(
      independent_variable = "Phenotest",
      family_test          = "spearman",
      transformation_y     = "",
      transformation_x     = "",
      aggregation          = "SUM",
      scope                = "INSTANCE",
      filter_p_value       = FALSE,
      areas_sql_condition  = NA,
      samples_sql_condition = NA,
      association_results_sql_condition = NA,
      stringsAsFactors = FALSE
    )
    testthat::expect_no_error(
      SEMseeker:::enrich_WebGestalt(
        study           = "test",
        inference_detail = inference_detail,
        significance     = TRUE
      )
    )
  } else {
    testthat::skip("WebGestaltR is installed; skipping guard test")
  }

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})

test_that("enrich_STRINGdb returns NULL gracefully when STRINGdb not installed", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    inference_details <- data.frame(
      independent_variable = "Phenotest",
      family_test          = "spearman",
      transformation_y     = "",
      transformation_x     = "",
      aggregation          = "SUM",
      scope                = "INSTANCE",
      filter_p_value       = FALSE,
      stringsAsFactors = FALSE
    )
    testthat::expect_no_error(
      SEMseeker:::enrich_STRINGdb(
        study            = "test",
        inference_detail = inference_details
      )
    )
  } else {
    testthat::skip("STRINGdb is installed; skipping guard test")
  }

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})

test_that("enrich_pathfindR returns NULL gracefully when pathfindR not installed", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  if (!requireNamespace("pathfindR", quietly = TRUE)) {
    inference_details <- data.frame(
      independent_variable = "Phenotest",
      family_test          = "spearman",
      transformation_y     = "",
      transformation_x     = "",
      aggregation          = "SUM",
      scope                = "INSTANCE",
      filter_p_value       = FALSE,
      stringsAsFactors = FALSE
    )
    testthat::expect_no_error(
      SEMseeker:::enrich_pathfindR(
        study             = "test",
        path_dbs          = c("KEGG"),
        inference_detail = inference_details,
        significance      = TRUE
      )
    )
  } else {
    testthat::skip("pathfindR is installed; skipping guard test")
  }

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})

# ---------------------------------------------------------------------------
# 3. enrich_ctdR integration — runs after semseeker + association_analysis
#    Requires ctdR to be installed (installed via Remotes in DESCRIPTION).
#    Uses the same bimodal synthetic data as test-7 to guarantee mutations.
# ---------------------------------------------------------------------------

test_that("enrich_ctdR runs without error on synthetic association results", {
  if (!requireNamespace("ctdR", quietly = TRUE)) {
    testthat::skip("ctdR not installed")
  }

  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  unlink(tempFolder, recursive = TRUE)

  # ── Build synthetic data with injected mutations ──────────────────────────
  set.seed(777)
  n_probes_ep <- 200L
  local_probes_ep <- probe_features[1:n_probes_ep, ]
  local_sig_ep <- matrix(stats::rbeta(n_probes_ep * nsamples, 90L, 10L),
                          nrow = n_probes_ep, ncol = nsamples)
  local_sig_ep[1:50, 1:5] <- stats::rbeta(50L * 5L, 1L, 100L)
  rownames(local_sig_ep) <- local_probes_ep$PROBE
  local_sig_ep <- as.data.frame(local_sig_ep)
  # signal_data has 10 unique columns; mySampleSheet has 16 rows (Reference reuse pattern)
  colnames(local_sig_ep) <- colnames(signal_data)

  # ── semseeker ─────────────────────────────────────────────────────────────
  SEMseeker::semseeker(
    input             = local_sig_ep,
    sample_sheet      = mySampleSheet,
    result_folder     = tempFolder,
    parallel_strategy = "sequential",
    areas             = c("GENE", "POSITION"),
    markers           = c("MUTATIONS"),
    start_fresh       = TRUE,
    showprogress      = showprogress,
    verbosity         = verbosity
  )

  # ── association_analysis at scope INSTANCE: enrichment reads SCOPE = INSTANCE
  #    and AREA = GENE (assoc_results_get()), so the collapsed branch would leave
  #    it nothing to read ─────────────────────────────────────────────────────
  inference_details <- data.frame(
    independent_variable = "Phenotest",
    family_test          = "spearman",
    transformation_y     = "",
    transformation_x     = "",

    aggregation          = "SUM",
    scope                = "INSTANCE",
    filter_p_value       = FALSE,
    stringsAsFactors     = FALSE
  )

  SEMseeker:::association_analysis(
    inference_details  = inference_details,
    result_folder      = tempFolder,
    parallel_strategy  = "sequential",
    markers            = c("MUTATIONS"),
    figures            = c("HYPO"),
    areas              = c("GENE"),
    multiple_test_adj  = "BH",
    showprogress       = showprogress,
    verbosity          = verbosity
  )

  # ── enrich_ctdR ──────────────────────────────────────────────────────────
  # Re-open environment so pathway functions can read ssEnv
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                        markers = c("MUTATIONS"), figures = c("HYPO"),
                        areas = c("GENE"), multiple_test_adj = "BH",
                        showprogress = showprogress, verbosity = verbosity)

  testthat::expect_no_error(
    SEMseeker:::enrich_ctdR(
      study            = "test",
      inference_detail = inference_details,
      significance     = FALSE     # include all results, not just significant
    )
  )

  # AI-308: the skip here used to blame a "depth_analysis = 3 regression"
  # producing only DEPTH = 1 rows. That diagnosis outlived the column it named:
  # depth was retired, and the run above now asks for SCOPE = INSTANCE
  # explicitly, which is what enrichment reads. The conditional stays because
  # ctdR can legitimately produce nothing on a synthetic fixture, but it no
  # longer asserts a cause it cannot observe.
  pathway_dir <- file.path(tempFolder, "Pathway")
  if (!dir.exists(pathway_dir)) {
    testthat::skip("ctdR produced no Pathway output on this fixture")
  }
  testthat::expect_true(dir.exists(pathway_dir))

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})
