# 2026-06-09: parity test: limma chunked-per-chr produces
# numerically identical results to a monolithic limma fit + eBayes
# on the same input.
#
# This is the load-bearing claim of the chunking design: cross-gene
# eBayes shrinkage is preserved when we lmFit per chr, concat the fit
# objects, and run eBayes once globally. If this test ever drifts past
# the 1e-9 tolerance, the chunking path has diverged from the
# monolithic and the dispatcher's automatic choice is no longer safe.

test_that("limma chunked per chr is bit-equal to monolithic (+ global eBayes)", {
  skip_if_not_installed("limma")
  skip_if_not_installed("polars")

  set.seed(20260609L)
  n_probes  <- 100L
  n_samples <- 30L
  n_chr     <- 5L   # 20 probes per chr

  # Simulate a pivot: 100 probes × 30 samples with mild differential
  # signal driven by a continuous IV. Probes 1..20 chr1, 21..40 chr2, …
  iv <- stats::rnorm(n_samples)
  signal_per_probe <- stats::rnorm(n_probes, sd = 0.4)
  Y <- matrix(NA_real_, nrow = n_probes, ncol = n_samples)
  for (g in seq_len(n_probes)) {
    Y[g, ] <- signal_per_probe[g] * iv + stats::rnorm(n_samples, sd = 0.6)
  }
  probe_ids <- paste0("cg", sprintf("%07d", seq_len(n_probes)))
  chr_per_probe <- rep(paste0("chr", seq_len(n_chr)),
                       each = n_probes %/% n_chr)
  rownames(Y) <- probe_ids
  colnames(Y) <- paste0("S", sprintf("%03d", seq_len(n_samples)))

  probe_features <- data.frame(
    PROBE = probe_ids, CHR = chr_per_probe,
    stringsAsFactors = FALSE
  )

  # Lazy pivot in the schema: AREA + sample cols.
  df <- data.frame(AREA = probe_ids, Y, check.names = FALSE,
                    stringsAsFactors = FALSE)
  pivot_lazy <- polars::as_polars_df(df)$lazy()

  # Same design for both paths.
  design <- cbind(`(Intercept)` = 1, IV_1 = iv)

  fit_mono <- SEMseeker:::assoc_lmfit_monolithic_lazy(
    pivot_lazy        = pivot_lazy,
    sample_cols_kept  = colnames(Y),
    design            = design,
    engine            = "limma",
    key               = list(MARKER = "TEST", FIGURE = "MEAN",
                             AREA = "PROBE", SUBAREA = "WHOLE"),
    family_test       = "limma_1"
  )
  expect_true(inherits(fit_mono, "MArrayLM"))

  fit_chunk <- SEMseeker:::assoc_lmfit_chunked_by_chr(
    pivot_lazy        = pivot_lazy,
    sample_cols_kept  = colnames(Y),
    design            = design,
    engine            = "limma",
    key               = list(MARKER = "TEST", FIGURE = "MEAN",
                             AREA = "PROBE", SUBAREA = "WHOLE"),
    family_test       = "limma_1",
    probe_features    = probe_features,
    tech_is_longread  = FALSE
  )
  expect_true(inherits(fit_chunk, "MArrayLM"))

  # The chunked fit returns probes in CHR-block order — same set,
  # different order than the monolithic one. Reorder by row name
  # before comparing.
  ord <- match(rownames(fit_mono$coefficients),
                rownames(fit_chunk$coefficients))
  fit_chunk$coefficients   <- fit_chunk$coefficients[ord, , drop = FALSE]
  fit_chunk$stdev.unscaled <- fit_chunk$stdev.unscaled[ord, , drop = FALSE]
  fit_chunk$sigma          <- fit_chunk$sigma[ord]
  fit_chunk$Amean          <- fit_chunk$Amean[ord]
  fit_chunk$df.residual    <- fit_chunk$df.residual[ord]

  # ---- the load-bearing assertions ----
  expect_equal(fit_chunk$coefficients,   fit_mono$coefficients,
                tolerance = 1e-9)
  expect_equal(fit_chunk$sigma,          fit_mono$sigma,
                tolerance = 1e-9)
  expect_equal(fit_chunk$Amean,          fit_mono$Amean,
                tolerance = 1e-9)
  expect_equal(fit_chunk$stdev.unscaled, fit_mono$stdev.unscaled,
                tolerance = 1e-9)
  expect_equal(as.numeric(fit_chunk$df.residual),
                as.numeric(fit_mono$df.residual), tolerance = 1e-9)

  # Now run global eBayes on both — the moderated t-stats / p-values
  # must match too, which is the whole reason we concat instead of
  # running per-chunk eBayes.
  e_mono  <- limma::eBayes(fit_mono)
  e_chunk <- limma::eBayes(fit_chunk)

  expect_equal(e_chunk$t,        e_mono$t,        tolerance = 1e-9)
  expect_equal(e_chunk$p.value,  e_mono$p.value,  tolerance = 1e-9)
  expect_equal(e_chunk$s2.post,  e_mono$s2.post,  tolerance = 1e-9)
})
