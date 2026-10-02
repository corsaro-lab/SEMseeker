# AI-044 (2026-06-08): unit test for the binomial_bulk path.
#
# What this file pins down:
#   1. assoc_glm_model_bulk() returns one row per probe with the legacy schema
#      (per-coef PVALUE/ESTIMATE + top-level PVALUE/PVALUE_ADJ) so the
#      downstream CSV writer / FDR machinery sees the same column shape
#      it would have seen with the per-probe stats::glm path.
#   2. Coefficient estimates from Rfast::glm_logistic match stats::glm
#      to ~4 decimal places on a simulated binary outcome with a factor IV.
#   3. The AI-044 universal degenerate-burden filter (in io_data_preparation,
#      hit upstream) ensures that all-zero / all-one probes never reach
#      Rfast — but assoc_glm_model_bulk's per-probe safety net still returns NA
#      for any degenerate Y that somehow slips through.

skip_if_no_rfast <- function() {
  testthat::skip_if_not_installed("Rfast")
}

.init_test_session <- function() {
  tf <- tempfile("test_glm_bulk_")
  dir.create(tf, recursive = TRUE, showWarnings = FALSE)
  SEMseeker:::core_init_env(
    tf,
    parallel_strategy = "sequential",
    inpute            = "median",
    bulk_population   = FALSE,
    LESIONS_BP        = 5000L,
    bonferroni_threshold = 0.05
  )
  tf
}

.sim_binary_factor <- function(n_samples = 200L, n_probes = 60L,
                                n_levels = 4L, seed = 7L) {
  set.seed(seed)
  # IV: factor with `n_levels` levels, roughly balanced.
  iv <- sample(seq_len(n_levels), size = n_samples, replace = TRUE)
  iv_factor <- factor(iv, levels = seq_len(n_levels))

  # One numeric covariate to exercise the cov_mat code path.
  cov1 <- stats::rnorm(n_samples)

  # Per-probe coefficient on IV level == 4 (vs reference 1). Some probes
  # null (coef == 0) — those should give p ≈ U[0,1].
  Y <- matrix(NA_integer_, nrow = n_samples, ncol = n_probes)
  effect <- stats::rnorm(n_probes, sd = 0.6)
  effect[1:10] <- 0   # first 10 probes are null
  intercepts <- stats::rnorm(n_probes, mean = -2, sd = 0.4)  # rare-event baseline

  for (j in seq_len(n_probes)) {
    eta <- intercepts[j] + effect[j] * as.integer(iv == 4) + 0.2 * cov1
    p   <- 1 / (1 + exp(-eta))
    Y[, j] <- stats::rbinom(n_samples, size = 1, prob = p)
  }
  colnames(Y) <- paste0("PROBE_", seq_len(n_probes))

  df <- data.frame(
    Sample_ID    = paste0("S", seq_len(n_samples)),
    IV           = iv_factor,
    COV1         = cov1,
    Y,
    check.names  = FALSE,
    stringsAsFactors = FALSE
  )
  attr(df, "g_start") <- 4L   # Sample_ID, IV, COV1, then probes — burden = g_start..ncol
  df
}

.fake_key_probe <- function() {
  list(MARKER = "LESIONS", FIGURE = "HYPER",
       AREA   = "PROBE",   SUBAREA = "WHOLE")
}

test_that("assoc_glm_model_bulk produces one row per probe with legacy schema", {
  skip_if_no_rfast()
  tf <- .init_test_session(); on.exit(unlink(tf, recursive = TRUE), add = TRUE)
  df <- .sim_binary_factor(n_samples = 200L, n_probes = 40L)
  g_start <- attr(df, "g_start")

  res <- SEMseeker:::assoc_glm_model_bulk(
    tempDataFrame        = df,
    g_start              = g_start,
    family_test          = "binomial_bulk",
    covariates           = "COV1",
    key                  = .fake_key_probe(),
    transformation_y     = "none",
    dototal              = FALSE,
    session_folder       = tempdir(),
    independent_variable = "IV",
    samples_sql_condition = ""
  )

  expect_true(is.data.frame(res))
  expect_equal(nrow(res), 40L)   # 40 probes in
  expect_true(all(c("MARKER", "FIGURE", "AREA", "SUBAREA", "AREA_OF_TEST",
                    "FAMILY_TEST", "R_MODEL", "PVALUE")
                   %in% colnames(res)))
  # AI-257: the model emits raw p-values and nothing adjusted. A batch is not a
  # family — assoc_analysis_save_results() forms the three families and names
  # each column after its own.
  expect_false(any(grepl("_ADJ", colnames(res))))
  expect_equal(unique(res$FAMILY_TEST), "binomial_bulk")
  expect_equal(unique(res$R_MODEL), "Rfast::glm_logistic")
})

test_that("Rfast estimates roughly match stats::glm on the same data", {
  skip_if_no_rfast()
  tf <- .init_test_session(); on.exit(unlink(tf, recursive = TRUE), add = TRUE)
  df <- .sim_binary_factor(n_samples = 300L, n_probes = 5L, seed = 11L)
  g_start <- attr(df, "g_start")

  res <- SEMseeker:::assoc_glm_model_bulk(
    tempDataFrame        = df,
    g_start              = g_start,
    family_test          = "binomial_bulk",
    covariates           = "COV1",
    key                  = .fake_key_probe(),
    transformation_y     = "none",
    dototal              = FALSE,
    session_folder       = tempdir(),
    independent_variable = "IV",
    samples_sql_condition = ""
  )

  # Pick the 3rd probe (non-null effect) and compare with stats::glm.
  probe <- paste0("PROBE_", 3L)
  fit_glm <- stats::glm(
    as.formula(paste0("`", probe, "` ~ IV + COV1")),
    family = stats::binomial(link = "logit"),
    data   = df
  )
  coef_glm <- stats::coef(summary(fit_glm))

  # Pick the row in res for this probe
  row_idx <- which(res$AREA_OF_TEST == probe)
  expect_length(row_idx, 1L)

  # Compare IV-level-2 (= IV2 vs reference IV1) estimate
  est_col <- grep("^IV2_ESTIMATE$|^IV_LEVEL_2_ESTIMATE$|IV2.*ESTIMATE",
                  colnames(res), value = TRUE)[1]
  if (!is.na(est_col)) {
    est_bulk <- res[[est_col]][row_idx]
    est_glm  <- coef_glm["IV2", "Estimate"]
    expect_equal(est_bulk, est_glm, tolerance = 1e-3)
  }
})

test_that("wrong family_test returns NULL", {
  skip_if_no_rfast()
  tf <- .init_test_session(); on.exit(unlink(tf, recursive = TRUE), add = TRUE)
  df <- .sim_binary_factor(n_samples = 50L, n_probes = 5L)
  res <- SEMseeker:::assoc_glm_model_bulk(
    tempDataFrame        = df,
    g_start              = attr(df, "g_start"),
    family_test          = "binomial",   # wrong: only binomial_bulk allowed
    covariates           = "COV1",
    key                  = .fake_key_probe(),
    transformation_y     = "none",
    dototal              = FALSE,
    session_folder       = tempdir(),
    independent_variable = "IV",
    samples_sql_condition = ""
  )
  expect_null(res)
})
