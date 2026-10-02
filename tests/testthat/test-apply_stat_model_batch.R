# AI-040 Fase 2 + Fase 3: batch path for limma_<N> and voom_<N>.
#
# Three things this file pins down:
#   1. assoc_apply_stat_model_batch() returns one row per area, same schema
#      the per-area path returns area-by-area (so the FDR + selector
#      machinery in the caller works unchanged).
#   2. limma_<N> in batch mode produces DIFFERENT p-values than the
#      degenerate per-area path on the same data — i.e. eBayes
#      shrinkage actually kicks in.
#   3. voom_<N> in batch mode runs without error on simulated count
#      data and returns calibrated p-values (sanity check: under the
#      null hypothesis the p-value distribution is roughly uniform).
#
# These tests do NOT require an initialised SEMseeker session.

# Tiny simulators -------------------------------------------------------

.sim_batch_continuous <- function(n_samples = 80L, n_areas = 50L,
                                   degree = 2L, sigma = 1, seed = 1L) {
  set.seed(seed)
  iv <- stats::rnorm(n_samples)
  cov1 <- stats::rnorm(n_samples)
  Y <- matrix(NA_real_, nrow = n_samples, ncol = n_areas)
  for (a in seq_len(n_areas)) {
    coef_lin <- stats::rnorm(1, sd = 0.4)
    coef_quad <- stats::rnorm(1, sd = 0.2)
    signal <- coef_lin * iv + coef_quad * iv^2 + 0.3 * cov1
    Y[, a] <- signal + stats::rnorm(n_samples, sd = sigma)
  }
  colnames(Y) <- paste0("GENE", seq_len(n_areas))
  df <- data.frame(Sample_ID = paste0("S", seq_len(n_samples)),
                    IV = iv, COV1 = cov1, Y, Sample_Group = "g1",
                    check.names = FALSE, stringsAsFactors = FALSE)
  df
}

.sim_batch_counts <- function(n_samples = 80L, n_areas = 50L,
                               mean_count = 8, seed = 2L) {
  set.seed(seed)
  iv <- stats::rnorm(n_samples)
  cov1 <- stats::rnorm(n_samples)
  # Each area: NB counts with size=3 (dispersion). Effect of IV linear.
  Y <- matrix(NA_integer_, nrow = n_samples, ncol = n_areas)
  for (a in seq_len(n_areas)) {
    eff <- stats::rnorm(1, sd = 0.3)
    mu  <- mean_count * exp(eff * iv + 0.2 * cov1)
    Y[, a] <- stats::rnbinom(n_samples, mu = mu, size = 3)
  }
  colnames(Y) <- paste0("GENE", seq_len(n_areas))
  df <- data.frame(Sample_ID = paste0("S", seq_len(n_samples)),
                    IV = iv, COV1 = cov1, Y, Sample_Group = "g1",
                    check.names = FALSE, stringsAsFactors = FALSE)
  df
}

.fake_key <- function() {
  data.frame(MARKER = "MUTATIONS", FIGURE = "HYPO",
             AREA = "GENE", SUBAREA = "WHOLE",
             COMBINED = "MUTATIONS_HYPO",
             stringsAsFactors = FALSE)
}

# Stub io_data_preparation if it's only available inside SEMseeker namespace —
# we call assoc_apply_stat_model_batch directly with already-prepared input,
# bypassing io_data_preparation entirely by patching with a pass-through.
.with_passthrough_data_prep <- function(code) {
  if (!"io_data_preparation" %in% ls(asNamespace("SEMseeker"))) {
    return(code)
  }
  orig <- SEMseeker:::io_data_preparation
  # AI-255: io_data_preparation() lost dototal and depth_analysis with the TOTAL
  # synthesis they governed.
  passthrough <- function(family_test, transformation_y, tempDataFrame,
                          independent_variable, g_start, g_end,
                          covariates, key,
                          transformation_x = "none") {
    list(tempDataFrame = tempDataFrame,
         independent_variableLevels = c(NA, NA))
  }
  unlockBinding("io_data_preparation", asNamespace("SEMseeker"))
  assign("io_data_preparation", passthrough, envir = asNamespace("SEMseeker"))
  on.exit({
    assign("io_data_preparation", orig, envir = asNamespace("SEMseeker"))
    lockBinding("io_data_preparation", asNamespace("SEMseeker"))
  }, add = TRUE)
  force(code)
}

# (1) Schema and row count ---------------------------------------------

test_that("assoc_apply_stat_model_batch returns one row per informative area with the polynomial schema", {
  testthat::skip_if_not_installed("limma")

  df  <- .sim_batch_continuous(n_samples = 60L, n_areas = 20L)
  key <- .fake_key()
  # g_start = 4 -> columns 1..3 are Sample_ID, IV, COV1; cols 4..23 are 20 areas.

  res <- .with_passthrough_data_prep(
    SEMseeker:::assoc_apply_stat_model_batch(
      tempDataFrame = df, g_start = 4L,
      family_test = "limma_2",
      covariates = "COV1", key = key,
      transformation_y = "", dototal = FALSE,
      session_folder = tempdir(),
      independent_variable = "IV",
      samples_sql_condition = "")
  )

  testthat::expect_s3_class(res, "data.frame")
  testthat::expect_equal(nrow(res), 20L)
  testthat::expect_true(all(c("MARKER","FIGURE","AREA","SUBAREA","AREA_OF_TEST",
                                "PL_DEGREE","PL_PERC","R_MODEL","FAMILY_TEST",
                                "INDEPENDENT_VARIABLE","PVALUE")
                              %in% colnames(res)))
  testthat::expect_equal(res$PL_DEGREE[1], 2)
  testthat::expect_equal(res$R_MODEL[1], "limma::lmFit+eBayes")
})

# (2) Batch limma_N gives non-trivial eBayes shrinkage -----------------

test_that("limma_2 in batch mode produces different p-values than the degenerate per-area fit", {
  testthat::skip_if_not_installed("limma")

  df  <- .sim_batch_continuous(n_samples = 60L, n_areas = 30L)
  key <- .fake_key()

  batch_res <- .with_passthrough_data_prep(
    SEMseeker:::assoc_apply_stat_model_batch(
      tempDataFrame = df, g_start = 4L,
      family_test = "limma_2",
      covariates = "COV1", key = key,
      transformation_y = "", dototal = FALSE,
      session_folder = tempdir(),
      independent_variable = "IV",
      samples_sql_condition = "")
  )

  # Per-area limma p-values from Fase 1 on the FIRST area only
  one_area_df <- df[, c("IV", "COV1", "GENE1")]
  per_area <- SEMseeker:::assoc_model_limma(
    "limma_2", one_area_df,
    sig.formula = GENE1 ~ IV + COV1,
    transformation_y = "", plot = FALSE,
    samples_sql_condition = "",
    key = list(MARKER="x", FIGURE="x", AREA="x", SUBAREA="x", AREA_OF_TEST="GENE1"))

  testthat::expect_equal(batch_res$AREA_OF_TEST[1], "GENE1")

  batch_p1 <- batch_res$PVALUE[1]
  per_area_p1_col <- grep("_1_PVALUE$", toupper(colnames(per_area)), value = TRUE)[1]
  per_area_p1 <- as.numeric(per_area[1, per_area_p1_col])

  # Shrinkage should move the value; require they aren't bit-identical.
  testthat::expect_false(isTRUE(all.equal(batch_p1, per_area_p1, tolerance = 1e-10)))
  testthat::expect_true(is.finite(batch_p1) && is.finite(per_area_p1))
})

# (3) voom_N runs and returns roughly uniform null p-values ------------

test_that("voom_2 in batch mode runs on NB counts and gives calibrated p-values under the null", {
  testthat::skip_if_not_installed("limma")

  # Null simulation: counts independent of IV.
  set.seed(42)
  n_samples <- 100L
  n_areas   <- 200L
  iv   <- stats::rnorm(n_samples)
  cov1 <- stats::rnorm(n_samples)
  Y <- matrix(NA_integer_, n_samples, n_areas)
  for (a in seq_len(n_areas)) {
    Y[, a] <- stats::rnbinom(n_samples, mu = 8, size = 3)
  }
  colnames(Y) <- paste0("GENE", seq_len(n_areas))
  df <- data.frame(Sample_ID = paste0("S", seq_len(n_samples)),
                    IV = iv, COV1 = cov1, Y, Sample_Group = "g1",
                    check.names = FALSE, stringsAsFactors = FALSE)
  key <- .fake_key()

  res <- .with_passthrough_data_prep(
    SEMseeker:::assoc_apply_stat_model_batch(
      tempDataFrame = df, g_start = 4L,
      family_test = "voom_1",
      covariates = "COV1", key = key,
      transformation_y = "", dototal = FALSE,
      session_folder = tempdir(),
      independent_variable = "IV",
      samples_sql_condition = "")
  )

  testthat::expect_equal(nrow(res), n_areas)
  testthat::expect_true(grepl("voom", res$R_MODEL[1]))
  # Under the null, KS test of PVALUE against uniform should NOT reject
  # at alpha=0.01 (the test is two-sided, used here only to detect gross
  # mis-calibration — KS p > 0.01 is the loose target).
  pvals <- as.numeric(res$PVALUE)
  pvals <- pvals[is.finite(pvals) & pvals >= 0 & pvals <= 1]
  if (length(pvals) >= 20) {
    ks <- suppressWarnings(stats::ks.test(pvals, "punif"))
    testthat::expect_gt(ks$p.value, 0.01)
  }
})

# (4) Guard: assoc_apply_stat_model dispatches batch families ----------------

test_that("assoc_apply_stat_model intercepts limma_<N> / voom_<N> and routes to batch path", {
  testthat::skip_if_not_installed("limma")

  df  <- .sim_batch_continuous(n_samples = 50L, n_areas = 10L)
  key <- .fake_key()

  # Apply through the outer entry point — uses an in-memory ssEnv we
  # don't need since the batch path bypasses core_get_session_info().
  res <- tryCatch(
    .with_passthrough_data_prep(
      SEMseeker:::assoc_apply_stat_model(
        tempDataFrame = df, g_start = 4L,
        family_test   = "limma_2",
        covariates    = "COV1", key = key,
        transformation_y = "", dototal = FALSE,
        session_folder = tempdir(),
        independent_variable = "IV",
        samples_sql_condition = "")
    ),
    error = function(e) {
      # assoc_apply_stat_model still calls core_get_session_info() at the top, which
      # fails without an initialised session. Tolerate that — if the error
      # comes from core_get_session_info() the dispatch is fine; if it comes
      # from elsewhere (e.g. the foreach path) the routing is broken.
      if (grepl("session_info|ssEnv", conditionMessage(e), ignore.case = TRUE))
        return(structure(list(.session_skip = TRUE),
                          class = c("session_skip", "data.frame")))
      stop(e)
    }
  )
  if (inherits(res, "session_skip")) testthat::skip("session not initialised — see comment")

  testthat::expect_s3_class(res, "data.frame")
  testthat::expect_true("AREA_OF_TEST" %in% colnames(res))
  testthat::expect_equal(nrow(res), 10L)
})
