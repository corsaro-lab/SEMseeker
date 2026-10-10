## Concordance benchmark: Illumina manifest vs WGBS pipeline (C-04 validation).
##
## Strategy: take a subset of Illumina K850 probes (ground truth annotation
## from the manifest), run the WGBS pipeline (anno_area_granges_build() +
## findOverlaps) on the same coordinates, and check that the two agree per
## area with sensible thresholds.
##
## Assertion strategy:
##   - bundled areas (CHR_CYTOBAND, DMR_*): strict identity (rate = 1.0).
##     Both pipelines use the same bundled RDA, so concordance must be perfect.
##   - annotationhub (ISLAND_*) and txdb (GENE_*): no pass/fail thresholds.
##     Label encodings differ between manifest (Islands_Name, UCSC RefGene)
##     and the WGBS path (AnnotationHub coords, TxDb knownGene symbols), so a
##     rigid string equality is not meaningful here. Instead the exploratory
##     test below invokes anno_concordance_report() over all areas
##     with assertions only on the structure of the report — the CSV it
##     produces is intended for paper figure `fig:illumina_wgbs_concordance`,
##     where human judgement can compare strict vs intersection rates per
##     area category.
##
## The test is heavy (TxDb load, optional AnnotationHub download). Guarded
## with skip_on_cran() + skip_if_not_installed().

.req_illumina <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"

# -------------------------------------------------------------------------
# Per-test setup/teardown: anno_probe_annotation_build() requires an initialised
# session (ssEnv). Wrap each call in a temp core_init_env()/core_close_env() scope.
# -------------------------------------------------------------------------
.with_session <- function(expr, envir = parent.frame()) {
  tempFolder <- tempfile(pattern = "semseeker_concordance_")
  SEMseeker:::core_init_env(
    result_folder     = tempFolder,
    parallel_strategy = "sequential",
    showprogress      = FALSE,
    verbosity         = 0
  )
  withr::defer({
    try(SEMseeker:::core_close_env(), silent = TRUE)
    unlink(tempFolder, recursive = TRUE)
  }, envir = envir)
  force(expr)
}

# =========================================================================
# Bundled areas — must be identical by construction (rate = 1.0)
# =========================================================================

test_that("annotation_concordance: CHR_CYTOBAND is 100% identical (bundled)", {
  skip_on_cran()
  skip_if_not_installed(.req_illumina)
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")

  report <- .with_session(SEMseeker:::anno_concordance_report(
    tech         = "K850",
    n_probes     = 500L,
    genome_build = "hg19",
    areas        = "CHR_CYTOBAND"
  ))

  expect_equal(report$category, "bundled")
  expect_equal(report$concordance_rate_strict, 1.0)
  expect_equal(report$concordance_rate_intersection, 1.0)
})

test_that("annotation_concordance: DMR_WHOLE is 100% identical (bundled)", {
  skip_on_cran()
  skip_if_not_installed(.req_illumina)
  skip_if_not_installed("GenomicRanges")

  report <- .with_session(SEMseeker:::anno_concordance_report(
    tech         = "K850",
    n_probes     = 500L,
    genome_build = "hg19",
    areas        = "DMR_WHOLE"
  ))

  expect_equal(report$category, "bundled")
  expect_equal(report$concordance_rate_strict, 1.0)
})

test_that("annotation_concordance: DMR_DMR is 100% identical (bundled)", {
  skip_on_cran()
  skip_if_not_installed(.req_illumina)
  skip_if_not_installed("GenomicRanges")

  report <- .with_session(SEMseeker:::anno_concordance_report(
    tech         = "K850",
    n_probes     = 500L,
    genome_build = "hg19",
    areas        = "DMR_DMR"
  ))

  expect_equal(report$category, "bundled")
  expect_equal(report$concordance_rate_strict, 1.0)
})

# =========================================================================
# Report structure + informative exploration over all areas.
# No assertion on rate values for txdb/annotationhub categories.
# =========================================================================

test_that("anno_concordance_report returns correct structure", {
  skip_on_cran()
  skip_if_not_installed(.req_illumina)
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")

  report <- .with_session(SEMseeker:::anno_concordance_report(
    tech         = "K850",
    n_probes     = 100L,
    genome_build = "hg19",
    areas        = c("DMR_WHOLE", "CHR_CYTOBAND")
  ))

  expected_cols <- c(
    "area", "category", "n_probes", "n_both_na", "n_both_labeled",
    "n_only_illumina", "n_only_wgbs",
    "n_label_match_strict", "n_label_match_intersection",
    "concordance_rate_strict", "concordance_rate_intersection"
  )
  expect_setequal(colnames(report), expected_cols)
  expect_equal(nrow(report), 2L)
  expect_true(all(report$n_probes == 100L))

  # Rows counted consistently:
  # n_both_na + n_both_labeled + n_only_illumina + n_only_wgbs == n_probes
  expect_true(all(
    report$n_both_na + report$n_both_labeled +
      report$n_only_illumina + report$n_only_wgbs == report$n_probes
  ))

  # Bundled areas: rate must be 1.0 (NA only if no probe has a label at all,
  # which is impossible for a 100-probe K850 subset on CHR_CYTOBAND/DMR_WHOLE).
  expect_true(all(report$concordance_rate_strict == 1.0, na.rm = TRUE))
})
