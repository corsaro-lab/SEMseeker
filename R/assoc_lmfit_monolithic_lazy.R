# 2026-06-09: monolithic lmFit path for the batch-lazy
# dispatcher. Extracted from `apply_stat_model_batch_lazy.R` so the
# dispatcher can pick this path or the chunked-per-chr variant based on
# the memory gate decision.
#
# Memory profile (the reason the gate gates this function):
#   y_mat        ≈ n_probes × n_samples × 8 B   (single biggest alloc)
#   lmFit work   ≈ 1.5 × y_mat                  (per-gene residuals)
#   fit object   ≈ n_probes × n_coef × 8 × 6 B  (after eBayes)
#
# Returns the pre-eBayes MArrayLM. The caller (dispatcher) is
# responsible for running `limma::eBayes()` globally so this function
# composes with the chunked path without behavioural divergence.

#' Run a single monolithic lmFit (optionally with voom) on a lazy pivot
#'
#' @param pivot_lazy `polars_lazy_frame` already filtered for
#'   `area_to_remove` and transformed by `io_data_preparation_lazy()`.
#'   First column MUST be `AREA` (probe / position identifier); the
#'   remaining columns are samples.
#' @param sample_cols_kept Character vector - sample column names that
#'   survived the IV / covariates complete-cases filter. Used to subset
#'   the pivot lazily before materialisation.
#' @param design Numeric matrix - the design passed to `limma::lmFit()`.
#' @param engine Character - `"limma"` or `"voom"`. Selects whether to
#'   run `limma::voom()` first.
#' @param key List with `MARKER`/`FIGURE`/`AREA`/`SUBAREA` - for log lines.
#' @param family_test Character - the inference_detail family identifier
#'   (`"limma_2"`, `"voom_2"`, ...). Passed through to log lines.
#'
#' @return The pre-eBayes `MArrayLM` fit object, or NULL on failure
#'   (NULL is silently propagated by the dispatcher).
#'
#' @keywords internal
#' @noRd
assoc_lmfit_monolithic_lazy <- function(pivot_lazy,
                                   sample_cols_kept,
                                   design,
                                   engine,
                                   key,
                                   family_test) {

  # Subset the lazy pivot to the kept sample columns (+ AREA) before
  # materialising. With monolithic + lazy + already-filtered rows, this
  # is the canonical place to allocate y_mat.
  pivot_lazy <- pivot_lazy$select(c("AREA", sample_cols_kept))
  pivot_df   <- pivot_lazy$collect()
  n_genes    <- pivot_df$height
  if (n_genes == 0L) {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_lmfit_monolithic_lazy: empty pivot after filter [",
              key$MARKER, "/", key$FIGURE, "/", key$AREA, "/",
              key$SUBAREA, "].")
    return(NULL)
  }

  gene_names <- as.character(pivot_df$select("AREA")$to_series())
  pivot_vals <- pivot_df$drop("AREA")
  sample_cols <- names(pivot_vals)
  rm(pivot_df)

  y_mat <- matrix(NA_real_, nrow = n_genes, ncol = length(sample_cols),
                   dimnames = list(gene_names, sample_cols))
  for (i in seq_along(sample_cols)) {
    y_mat[, i] <- as.vector(pivot_vals$select(sample_cols[i])$to_series())
  }
  rm(pivot_vals)
  gc(verbose = FALSE)
  y_mat[is.na(y_mat)] <- 0

  core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
            " assoc_lmfit_monolithic_lazy [", key$MARKER, "/", key$FIGURE, "/",
            key$AREA, "/", key$SUBAREA, "] y_mat=", n_genes, "x",
            ncol(y_mat), " engine=", engine)

  fit <- if (engine == "voom") {
    voom_obj <- tryCatch(limma::voom(y_mat, design),
                          error = function(e) {
                            core_log_event("ERROR: ",
                                      format(Sys.time(), "%a %b %d %X %Y"),
                                      " assoc_lmfit_monolithic_lazy voom failed: ",
                                      conditionMessage(e))
                            NULL
                          })
    if (is.null(voom_obj)) return(NULL)
    f <- limma::lmFit(voom_obj, design)
    rm(voom_obj); gc(verbose = FALSE)
    f
  } else {
    tryCatch(limma::lmFit(y_mat, design),
              error = function(e) {
                core_log_event("ERROR: ",
                          format(Sys.time(), "%a %b %d %X %Y"),
                          " assoc_lmfit_monolithic_lazy lmFit failed: ",
                          conditionMessage(e))
                NULL
              })
  }
  rm(y_mat); gc(verbose = FALSE)
  fit
}
