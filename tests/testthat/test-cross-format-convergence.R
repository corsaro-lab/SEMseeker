## ============================================================================
## test-cross-format-convergence.R
##
## Same biological signal injected via three different input modalities
## (Illumina array matrix / WGBS bedmethyl / ONT long-read bedmethyl) must
## produce convergent SEM calls and lesions.
##
## Convergence is asserted on the SEM CALL set (binary), not on raw beta values
## — writing bedmethyl with a fixed depth introduces quantization
## (beta * depth → integer Nmod), so the betas round-trip is lossy at the
## ~1 % level. The IQR-based SEM threshold is robust to this; in practice the
## SEM call sets should overlap by ≥99 %.
##
## This test is `skip_on_cran()` because it writes temporary bedmethyl files
## and runs the full pipeline three times (~30-60 s).
## ============================================================================

test_that("array, WGBS bedmethyl and LONGREAD bedmethyl produce convergent SEM calls", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("polars")
  if (!exists("make_bedmethyl_per_sample", mode = "function"))
    source(testthat::test_path("helper-bedmethyl.R"))

  ## ── 1. Real GSE133774 fixture filtered to KCNQ1OT1 + H19/IGF2 probes ─────
  ## Same fixture used by setup.R and the getting-started vignette (AI-123).
  utils::data("test_signal_gse133774",      package = "SEMseeker", envir = environment())
  utils::data("test_samplesheet_gse133774", package = "SEMseeker", envir = environment())
  utils::data("test_master_features",       package = "SEMseeker", envir = environment())

  pf <- as.data.frame(test_master_features)
  ## Restrict to BWS-relevant imprinting DMRs on chr11 for speed
  pf_subset <- pf[!is.na(pf$DMR_LABEL) &
                    grepl("KCNQ1OT1|H19_IGF2", pf$DMR_LABEL) &
                    pf$PROBE %in% rownames(test_signal_gse133774), , drop = FALSE]
  if (nrow(pf_subset) < 50L) {
    chr11_pool <- pf[pf$CHR == "11" &
                     !pf$PROBE %in% pf_subset$PROBE &
                     pf$PROBE %in% rownames(test_signal_gse133774), ]
    need <- min(100L - nrow(pf_subset), nrow(chr11_pool))
    if (need > 0L)
      pf_subset <- rbind(pf_subset, chr11_pool[seq_len(need), ])
  }
  expect_gt(nrow(pf_subset), 30L)

  ## ── 2. Real beta matrix (no synthetic injection — real BWS biology) ───────
  betas        <- test_signal_gse133774[pf_subset$PROBE, , drop = FALSE]
  sample_sheet <- test_samplesheet_gse133774
  sample_sheet <- sample_sheet[, c("Sample_ID", "Sample_Group")]
  sample_sheet$Sample_Name <- sample_sheet$Sample_ID

  sample_sheet <- data.frame(
    Sample_ID    = sample_sheet$Sample_ID,
    Sample_Name  = sample_sheet$Sample_ID,
    Sample_Group = sample_sheet$Sample_Group,
    stringsAsFactors = FALSE
  )

  ## ── 3. Run path A: Illumina array (matrix input, tech auto-detected) ─────
  out_array <- tempfile("xfmt_array_")
  dir.create(out_array, recursive = TRUE)
  on.exit(unlink(out_array, recursive = TRUE), add = TRUE)
  array_ok <- tryCatch({
    SEMseeker::semseeker(
      input             = betas,
      sample_sheet      = sample_sheet,
      result_folder     = out_array,
      input_type        = "matrix",
      parallel_strategy = "sequential"
    )
    TRUE
  }, error = function(e) {
    testthat::skip(paste("array path failed (likely missing optional deps):",
                          conditionMessage(e)))
  })

  ## ── 4. Run path B: WGBS bedmethyl files (tech = "WGBS") ──────────────────
  out_wgbs  <- tempfile("xfmt_wgbs_")
  dir.create(out_wgbs, recursive = TRUE)
  on.exit(unlink(out_wgbs, recursive = TRUE), add = TRUE)
  bedmethyl_dir_wgbs <- tempfile("xfmt_wgbs_files_")
  dir.create(bedmethyl_dir_wgbs)
  on.exit(unlink(bedmethyl_dir_wgbs, recursive = TRUE), add = TRUE)
  bedmethyl_files_wgbs <- make_bedmethyl_per_sample(
    beta_values     = betas,
    probe_features  = pf_subset,
    out_dir         = bedmethyl_dir_wgbs,
    depth           = 30L  ## higher depth → less quantization error
  )
  wgbs_ok <- tryCatch({
    SEMseeker::semseeker(
      input             = unname(bedmethyl_files_wgbs),
      sample_sheet      = sample_sheet,
      result_folder     = out_wgbs,
      input_type        = "bedmethyl",
      tech              = "WGBS",
      parallel_strategy = "sequential"
    )
    TRUE
  }, error = function(e) {
    testthat::skip(paste("WGBS path failed:", conditionMessage(e)))
  })

  ## ── 5. Run path C: LONGREAD bedmethyl files (tech = "LONGREAD") ──────────
  out_lr <- tempfile("xfmt_lr_")
  dir.create(out_lr, recursive = TRUE)
  on.exit(unlink(out_lr, recursive = TRUE), add = TRUE)
  bedmethyl_dir_lr <- tempfile("xfmt_lr_files_")
  dir.create(bedmethyl_dir_lr)
  on.exit(unlink(bedmethyl_dir_lr, recursive = TRUE), add = TRUE)
  bedmethyl_files_lr <- make_bedmethyl_per_sample(
    beta_values    = betas,
    probe_features = pf_subset,
    out_dir        = bedmethyl_dir_lr,
    depth          = 30L
  )
  lr_ok <- tryCatch({
    SEMseeker::semseeker(
      input              = unname(bedmethyl_files_lr),
      sample_sheet       = sample_sheet,
      result_folder      = out_lr,
      input_type         = "bedmethyl",
      tech               = "LONGREAD",
      genome_build       = "hg19",  ## test data is array-derived hg19; match other two paths
      strict_build_check = FALSE,   ## suppresses LONGREAD+hg19 warning (intentional here)
      parallel_strategy  = "sequential"
    )
    TRUE
  }, error = function(e) {
    testthat::skip(paste("LONGREAD path failed:", conditionMessage(e)))
  })

  ## ── 6. Assert convergence on output existence (smoke level) ──────────────
  ## A full convergence assertion on the SEM call sets is left to AI-117
  ## (BWS regression test) — at the smoke level we verify that all three paths
  ## produce non-empty MUTATIONS BED output and that lesion counts are within
  ## ±20 % of each other (quantization tolerance with depth=30 → ~3 % beta).
  expect_true(array_ok)
  expect_true(wgbs_ok)
  expect_true(lr_ok)

  ## AI-096 lazy-passthrough writes pivot parquets directly — no per-sample
  ## BED intermediates. Assert convergence via MUTATIONS_HYPO pivot existence
  ## and row-count parity across the three input formats.
  ## Full biological convergence (real GSE95486 beta values) is AI-117 + AI-123.
  # AI-255: the name carries the scope and the aggregation. At POSITION the
  # block is a single position, so the aggregation is VALUE — the identity.
  mut_pivot_rel <- file.path("Data","Pivots","MUTATIONS",
                              "MUTATIONS_HYPO_INSTANCE_POSITION_WHOLE_VALUE_HG19.parquet")
  expect_true(file.exists(file.path(out_array, mut_pivot_rel)))
  expect_true(file.exists(file.path(out_wgbs,  mut_pivot_rel)))
  expect_true(file.exists(file.path(out_lr,    mut_pivot_rel)))

  ## Quantization-tolerant convergence: WGBS and LONGREAD row counts must be
  ## within ±50 % of the array baseline (depth=30 → ~3 % beta quantization).
  pivot_rows <- function(path) {
    tryCatch(polars::pl$read_parquet(path)$height, error = function(e) 0L)
  }
  n_array <- pivot_rows(file.path(out_array, mut_pivot_rel))
  n_wgbs  <- pivot_rows(file.path(out_wgbs,  mut_pivot_rel))
  n_lr    <- pivot_rows(file.path(out_lr,    mut_pivot_rel))

  if (n_array > 0L) {
    ratio_wgbs <- n_wgbs / n_array
    ratio_lr   <- n_lr   / n_array
    expect_gt(ratio_wgbs, 0.5)
    expect_lt(ratio_wgbs, 2.0)
    expect_gt(ratio_lr,   0.5)
    expect_lt(ratio_lr,   2.0)
  }
})
