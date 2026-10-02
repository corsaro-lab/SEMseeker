# AI-061+ (2026-06-09): memory gate for the batch-lazy lmFit dispatcher.
#
# Decides whether the limma::lmFit call should run monolithic (single
# y_mat materialisation) or chunked-per-chromosome (one chr at a time,
# concat fit objects, then a single global eBayes call). Returns also
# the breakdown so the dispatcher can log it and the caller can fail
# fast with an actionable error when even per-chr chunking would
# exceed the budget.
#
# Pattern mirrors `.sem_knn_memory_gate()` (AI-096 Phase 2): explicit
# estimates of every contributor, fraction-of-RAM budget tunable via
# env var, "fail-fast with all four numbers" diagnostic on overrun.
#
# Decision logic:
#   - monolithic  : monolithic y_mat + lmFit fits in budget (cheapest path)
#   - chunked     : monolithic would OOM, but the biggest chunk's y_mat +
#                   lmFit fits in budget - proceed with chunked dispatch
#   - abort       : even the biggest chunk's y_mat would OOM - there is no
#                   safe path on this machine; raise the env var or move
#                   to a bigger node.
#
# Env var: SEMSEEKER_BULK_MODEL_MEM_FRACTION (default 0.6). Sets the fraction
# of total system RAM that limma fitting is allowed to use as its peak.
# 0.6 leaves room for the OS, Polars cache, R working set and other
# processes - matches the SEMSEEKER_KNN_MEM_FRACTION convention.

#' Memory gate for the AI-061 batch-lazy lmFit path
#'
#' @param n_probes  integer. Total probe/area count of the input pivot
#'   (after the AI-044 degenerate-row filter has been applied lazily).
#' @param n_samples integer. Sample count after dropping rows with NA in
#'   the IV / covariates.
#' @param n_coef    integer. Number of columns in the lmFit design matrix
#'   - typically `1 + degree + length(covariates)`.
#' @param max_chunk_probes integer or NULL. Size of the largest chunk
#'   under the chunking plan. For Illumina + per-chr chunking this is
#'   `max(probes_per_chr)` (~16k-20k). For long-reads it is the maximum
#'   positions found on a single chromosome (variable). When NULL the
#'   gate assumes a default of `max(50000L, n_probes / 22L)` so the
#'   chunked-path estimate is conservative.
#'
#' @return Named list with the decision and the numeric breakdown:
#'   - `decision`     : `"monolithic"` / `"chunked"` / `"abort"`
#'   - `y_mat_GB`     : peak GB of the monolithic Y matrix
#'   - `lmfit_GB`     : peak GB of the monolithic lmFit working set
#'   - `mono_peak_GB` : sum of the two above
#'   - `chunk_y_GB`   : peak GB of the biggest chunk's Y matrix
#'   - `chunk_peak_GB`: sum of chunk_y_GB + lmFit working set
#'   - `fit_object_GB`: GB held by the concatenated MArrayLM after fitting
#'                     (`n_probes × n_coef × 8 × 6`, includes eBayes outputs)
#'   - `total_RAM_GB` : detected system RAM (`.sem_total_ram_GB()`)
#'   - `mem_frac`     : the env-var fraction in effect (default 0.6)
#'   - `available_GB` : `total_RAM_GB × mem_frac`
#'
#' @keywords internal
#' @noRd
.core_bulk_model_memory_gate <- function(n_probes,
                                     n_samples,
                                     n_coef,
                                     max_chunk_probes = NULL) {

  n_probes  <- as.numeric(n_probes)
  n_samples <- as.numeric(n_samples)
  n_coef    <- as.numeric(n_coef)
  byte_per_double <- 8

  # Monolithic estimates.
  y_mat_GB    <- (n_probes * n_samples * byte_per_double) / (1024^3)
  # Empirical: limma::lmFit's working set on a (probe × sample) matrix
  # peaks at ~1.5× the response matrix (residuals + per-gene sigma).
  lmfit_GB    <- 1.5 * y_mat_GB
  mono_peak_GB <- y_mat_GB + lmfit_GB

  # Chunked estimates. Default: max(50k probes per chunk, n_probes/22)
  # - a conservative ceiling that absorbs the biggest chromosome on
  # both Illumina (chr1 ≈ 22k on K850) and long-reads (chr1 can carry
  # 10⁵ - 10⁶ positions).
  if (is.null(max_chunk_probes) || !is.finite(max_chunk_probes) ||
      max_chunk_probes <= 0) {
    max_chunk_probes <- max(50000L, ceiling(n_probes / 22L))
  }
  max_chunk_probes <- as.numeric(max_chunk_probes)
  chunk_y_GB    <- (max_chunk_probes * n_samples * byte_per_double) / (1024^3)
  chunk_peak_GB <- chunk_y_GB + 1.5 * chunk_y_GB

  # Fit object after concatenation + eBayes:
  # 6 dense (n_probes × n_coef) matrices - coefficients, stdev.unscaled,
  # t-stat, p.value, lods, s2.post stretched on coef axis - plus the
  # sigma/Amean/df.residual vectors (small, rolled in). 6 is a tight
  # over-estimate.
  fit_object_GB <- (n_probes * n_coef * byte_per_double * 6) / (1024^3)

  # Memory budget.
  mem_frac <- suppressWarnings(as.numeric(
    Sys.getenv("SEMSEEKER_BULK_MODEL_MEM_FRACTION", "0.6")))
  if (is.na(mem_frac) || mem_frac <= 0 || mem_frac > 1) mem_frac <- 0.6
  total_GB <- .sem_total_ram_GB()
  if (is.na(total_GB) || total_GB <= 0) {
    # Without RAM detection we can't decide. Conservative default:
    # assume monolithic is fine (legacy behaviour) so the gate doesn't
    # block users on a portability issue. Operator can detect a true
    # OOM kill via dmesg / jetsam and override the env var.
    return(list(
      decision     = "monolithic",
      y_mat_GB     = y_mat_GB,
      lmfit_GB     = lmfit_GB,
      mono_peak_GB = mono_peak_GB,
      chunk_y_GB   = chunk_y_GB,
      chunk_peak_GB= chunk_peak_GB,
      fit_object_GB= fit_object_GB,
      total_RAM_GB = NA_real_,
      mem_frac     = mem_frac,
      available_GB = NA_real_
    ))
  }
  available_GB <- total_GB * mem_frac

  # AI-061+ empirical safety factor on the gate threshold.
  # The 1.5× lmfit_GB factor (line 71) under-estimates real peak because:
  #   - polars→R hand-off materialises 1-2 extra wide-frame copies during
  #     as.data.frame(collect()) before lmFit even starts
  #   - limma internals (residuals + sigma2 + Amean + qr.Q ...) routinely
  #     run another ~0.5× on top of the response matrix
  #   - eBayes post-fit objects (lods, t, p, B, var.post) add up to fit_object_GB
  #     ALREADY accounted, but are produced WHILE residuals are still alive
  #   - the R process carries a non-trivial baseline (polars cache + R
  #     working set accumulated by previous dispatches): mem available for
  #     the *next* lmFit peak is < available_GB by ~current RSS
  #
  # Observed run/factor pairs (refine on each over/under-shoot):
  #   - SAFETY_FACTOR=1.0 (no factor, original commit):
  #       v49 ewas, limma_2 DELTARQ HYPER PROBE WHOLE, mono=25.1 / avail=38.4
  #       → gate=monolithic → SIGKILL by macOS Jetsam.
  #   - SAFETY_FACTOR=1.5 (commit 98fc2a0):
  #       v50 ewas, limma_2 DELTARP HYPER PROBE WHOLE, mono=25.1 / avail=38.4
  #       → gate=monolithic (1.5 × 25.23 = 37.85 ≤ 38.4 - passes by 0.5 GB)
  #       → SIGKILL. The baseline RSS (~28 GB before fit started) ate the
  #       margin: real peak ≈ 53 GB > 64 GB - swap pressure.
  #   - SAFETY_FACTOR=2.0 (this revision):
  #       2.0 × 25.23 = 50.5 GB > 38.4 → chunked (was monolithic). Should
  #       hold until the next miss; document the next OOM here.
  SAFETY_FACTOR <- 2.0

  decision <- if (SAFETY_FACTOR * (mono_peak_GB  + fit_object_GB) <= available_GB) {
    "monolithic"
  } else if (SAFETY_FACTOR * (chunk_peak_GB + fit_object_GB) <= available_GB) {
    "chunked"
  } else {
    "abort"
  }

  list(
    decision     = decision,
    y_mat_GB     = y_mat_GB,
    lmfit_GB     = lmfit_GB,
    mono_peak_GB = mono_peak_GB,
    chunk_y_GB   = chunk_y_GB,
    chunk_peak_GB= chunk_peak_GB,
    fit_object_GB= fit_object_GB,
    total_RAM_GB = total_GB,
    mem_frac     = mem_frac,
    available_GB = available_GB
  )
}
