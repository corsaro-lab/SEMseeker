# Resuming an enrichment already on disk.
#
# The resume path is the one road that reaches both the result writer,
# enrich_result_save(), and the old "exists('existing_db')" branch of the
# enrichers without calling the remote service. Two defects lived there:
#   - enrich_result_save() read `study` as a free name and failed on every
#     non-empty result (fixed in #49, never tested until now);
#   - the enrichers carried a branch keyed on exists("existing_db"): exists()
#     searches up to the global environment, so a user object of that name
#     switched it on and the next line read `pathway_report`, which nothing
#     defines.
# Expected values come from the fixture (the phenotype flags, written by hand)
# and from stats::p.adjust() called here, never from SEMseeker.

test_that("resuming a ctdR result rewrites it with phenotype flags and BH p-values", {
  skip_if_not_installed("ctdR")

  result_folder <- sem_test_folder()
  ssEnv <- SEMseeker:::core_init_env(result_folder, parallel_strategy = "sequential", markers = "MUTATIONS",
                                     figures = "HYPO", areas = "GENE", multiple_test_adj = "BH")
  on.exit(SEMseeker:::core_close_env(), add = TRUE)

  inference_detail <- data.frame(independent_variable = "Phenotest", family_test = "spearman",
                                 transformation_y = "", transformation_x = "", aggregation = "SUM",
                                 scope = "INSTANCE", filter_p_value = FALSE, stringsAsFactors = FALSE)

  # One result already on disk for every key the enrichment will visit, in
  # the ctdR format: three chemicals, two of them named after the study.
  raw_p_values <- c(0.04, 0.001, 0.2)
  stored_result <- data.frame(ChemicalID = c("C1", "C2", "C3"),
                              ChemicalName = c("Breast cancer agent A", "Tamoxifen for BREAST tumours", "Caffeine"),
                              pvalue = raw_p_values, padj = raw_p_values, foldEnrichment = c(2, 5, 1.1),
                              source = "ctd", stringsAsFactors = FALSE)
  keys <- unique(ssEnv$keys_for_pathway)
  result_paths <- vapply(seq_len(nrow(keys)), function(i) {
    result_name <- SEMseeker:::enrich_phenotype_analysis_name(inference_detail, keys[i, ], prefix = "", suffix = "without_signal_",
                                                               "PVALUE_ADJ_ALL_BH", ssEnv$alpha, TRUE)
    SEMseeker:::io_file_path_build(SEMseeker:::io_enrichment_folder(inference_detail, "ctdR"), result_name, "csv")
  }, character(1))
  for (result_path in result_paths) {
    dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv2(stored_result, result_path, row.names = FALSE)
  }

  # A user object with the name the old branch looked for.
  assign("existing_db", "anything", envir = globalenv())
  on.exit(rm("existing_db", envir = globalenv()), add = TRUE)

  expect_no_error(SEMseeker:::enrich_ctdR(study = "breast", inference_detail = inference_detail))

  expected_order <- order(stats::p.adjust(raw_p_values, method = "BH"))
  for (result_path in result_paths) {
    rewritten <- utils::read.csv2(result_path, stringsAsFactors = FALSE)
    expect_identical(rewritten$ChemicalID, stored_result$ChemicalID[expected_order])
    expect_identical(rewritten$PHENOTYPE, c(TRUE, TRUE, FALSE)[expected_order])
    expect_equal(rewritten$PVALUE_ADJ_BH, stats::p.adjust(raw_p_values, method = "BH")[expected_order], tolerance = 1e-12)
  }
})
