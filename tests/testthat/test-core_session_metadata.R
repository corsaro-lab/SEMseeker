## Coverage for session-metadata provenance helpers and quantreg metrics.
##
## Covered:
##   core_session_metadata_write()        session_metadata.json provenance
##   core_pivot_sidecar_write()           per-pivot _meta.json sidecar
##   core_check_session_compatibility()   cross-session build/tech enforcement
##   assoc_quantreg_metrics()             pinball loss / QQ metrics, and its charts
##

# ---------------------------------------------------------------------------
# core_session_metadata_write
# ---------------------------------------------------------------------------

test_that("core_session_metadata_write serialises provenance to JSON", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  meta <- SEMseeker:::core_session_metadata_write(tf, sample_n = 7L)

  jf <- file.path(tf, "session_metadata.json")
  expect_true(file.exists(jf))
  parsed <- jsonlite::fromJSON(jf)
  expect_equal(parsed$sample_n, 7L)
  expect_true(nzchar(parsed$genome_build))
  expect_true(nzchar(parsed$semseeker_version))
  expect_equal(meta$sample_n, 7L)   # returned invisibly
})

# ---------------------------------------------------------------------------
# core_pivot_sidecar_write
# ---------------------------------------------------------------------------

test_that("core_pivot_sidecar_write writes a _meta.json next to the pivot", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  pivot   <- file.path(tf, "MUTATIONS_HYPER_GENE_TSS1500_hg19.parquet")
  SEMseeker:::core_pivot_sidecar_write(pivot)

  sidecar <- file.path(tf, "MUTATIONS_HYPER_GENE_TSS1500_hg19_meta.json")
  expect_true(file.exists(sidecar))
  parsed <- jsonlite::fromJSON(sidecar)
  expect_equal(parsed$pivot_file, "MUTATIONS_HYPER_GENE_TSS1500_hg19.parquet")
  expect_true(nzchar(parsed$genome_build))
})

# ---------------------------------------------------------------------------
# core_check_session_compatibility
# ---------------------------------------------------------------------------

.write_session_meta <- function(dir, genome_build, tech) {
  writeLines(
    jsonlite::toJSON(list(genome_build = genome_build, tech = tech),
                     auto_unbox = TRUE),
    file.path(dir, "session_metadata.json"))
}

test_that("core_check_session_compatibility accepts matching provenance", {
  d1 <- tempfile("studyA"); dir.create(d1)
  d2 <- tempfile("studyB"); dir.create(d2)
  on.exit(unlink(c(d1, d2), recursive = TRUE), add = TRUE)
  .write_session_meta(d1, "hg19", "K450")
  .write_session_meta(d2, "hg19", "K450")

  res <- SEMseeker:::core_check_session_compatibility(c(d1, d2))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2L)
  expect_true(all(res$genome_build == "hg19"))
})

test_that("core_check_session_compatibility stops on differing genome build", {
  d1 <- tempfile("studyA"); dir.create(d1)
  d2 <- tempfile("studyB"); dir.create(d2)
  on.exit(unlink(c(d1, d2), recursive = TRUE), add = TRUE)
  .write_session_meta(d1, "hg19", "K450")
  .write_session_meta(d2, "hg38", "K450")

  expect_error(SEMseeker:::core_check_session_compatibility(c(d1, d2)),
               "genome_build differs")
})

test_that("core_check_session_compatibility warns on differing tech", {
  d1 <- tempfile("studyA"); dir.create(d1)
  d2 <- tempfile("studyB"); dir.create(d2)
  on.exit(unlink(c(d1, d2), recursive = TRUE), add = TRUE)
  .write_session_meta(d1, "hg19", "K450")
  .write_session_meta(d2, "hg19", "K850")

  expect_warning(SEMseeker:::core_check_session_compatibility(c(d1, d2)),
                 "Technologies differ")
})

test_that("core_check_session_compatibility is a no-op for a single session", {
  expect_null(SEMseeker:::core_check_session_compatibility("only_one_folder"))
})

# ---------------------------------------------------------------------------
# assoc_quantreg_metrics
# ---------------------------------------------------------------------------

test_that("assoc_quantreg_metrics computes pinball loss and QQ metrics", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  out <- SEMseeker:::assoc_quantreg_metrics(
    predicted_values     = c(1, 2, 3),
    expected_values      = c(1, 2, 4),
    tau                  = 0.5,
    res                  = data.frame(seed = 1),
    family_test          = "quantreg_0.5",
    independent_variable = "x",
    transformation_y     = "",
    dependent_variable   = "y",
    plot                 = FALSE)

  expect_true(all(c("pinball_loss", "below_quantile", "qq_distance",
                    "qq_correlation") %in% colnames(out)))
  # residuals = c(0,0,1); tau=0.5 -> only the +1 residual contributes 0.5
  expect_equal(out$pinball_loss, 0.5)
  expect_equal(out$below_quantile, 2 / 3)
})

test_that("assoc_quantreg_metrics with plot = TRUE writes both fitted-model charts", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  out <- SEMseeker:::assoc_quantreg_metrics(
    predicted_values     = c(1.1, 2.3, 3.0, 4.5),
    expected_values      = c(1.0, 2.0, 3.0, 4.0),
    tau                  = 0.5,
    res                  = data.frame(seed = 1),
    family_test          = "quantreg_0.5",
    independent_variable = "x",
    transformation_y     = "",
    dependent_variable   = "y",
    permutation_vector   = c(0.1, 0.2, 0.3, 0.4),
    plot                 = TRUE)

  # residuals = c(-0.1,-0.3,0,-0.5), all <= 0: (1 - 0.5) * 0.9 = 0.45;
  # every observed value is at or below its prediction, so 1, and |0.5 - 1|
  expect_equal(out$pinball_loss, 0.45)
  expect_equal(out$below_quantile, 1)
  expect_equal(out$qq_distance, 0.5)

  charts <- list.files(file.path(tf, "Chart", "FITTED_MODEL"))
  expect_setequal(charts, c("QUANTREG_0_5_X_VS_Y.png",
                            "QUANTREG_0_5_X_VS_Y_HISTOGRAM.png"))
})

test_that("assoc_quantreg_metrics with plot = TRUE and no data warns and draws nothing", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  warnings <- character(0)
  withCallingHandlers(
    SEMseeker:::assoc_quantreg_metrics(
      predicted_values     = numeric(0),
      expected_values      = numeric(0),
      tau                  = 0.5,
      res                  = data.frame(seed = 1),
      family_test          = "quantreg_0.5",
      independent_variable = "x",
      transformation_y     = "",
      dependent_variable   = "y",
      permutation_vector   = numeric(0),
      plot                 = TRUE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  expect_true("No data available for plotting observed vs. predicted." %in% warnings)
  expect_true("No data available for plotting histogram of betas." %in% warnings)
  expect_length(list.files(file.path(tf, "Chart", "FITTED_MODEL")), 0L)
})
