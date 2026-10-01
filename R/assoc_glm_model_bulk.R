# AI-044 (2026-06-08): bulk path for logistic regression - family_test
# "binomial_bulk". Mirror of apply_stat_model_batch.R (limma path) but
# for binomial GLM: per-probe Rfast::glm_logistic with shared design
# matrix, parallelised via foreach %dorng%.
#
# Why a bulk path:
#   - Per-probe stats::glm via foreach has heavy R-level overhead
#     (~5 ms per fit × 600k  probes x 4 inference cycles ≈ hours).
#   - Rfast::glm_logistic is a C++ Newton-Raphson implementation,
#     ~10-20× faster than stats::glm; combined with the AI-044
#     degenerate-burden filter in io_data_preparation (~92% of LESIONS
#     probes removed before reaching here) we get ~2-3 min total
#     instead of ~60 min for the scheda 3a-NEW binomial dispatch.
#
# Architectural placement: apply_stat_model.R early-returns into this
# function when family_test == "binomial_bulk" (analog of the limma/voom
# early-return). Per-probe binomial via family_test == "binomial" stays
# on the legacy foreach + stats::glm path for backward compatibility.
#
# Output schema matches the legacy binomial path exactly so the
# downstream CSV writer / FDR / selector machinery in assoc_apply_stat_model's
# caller works without changes:
#   INDIPENDENT_VARIABLE, MARKER, FIGURE, AREA, SUBAREA, AREA_OF_TEST,
#   FAMILY_TEST, transformation_y, COVARIATES, R_MODEL,
#   <coef>_PVALUE / <coef>_ESTIMATE per design column,
#   PVALUE (= first non-intercept p-value), PVALUE_ADJ (BH).

#' Bulk logistic regression across many probes with a shared design.
#'
#' @param tempDataFrame data.frame; rows = samples, columns = IV +
#'   covariates + one column per probe (the 0/1 burden).
#' @param g_start integer; index of the first probe column (everything
#'   before is sample-level metadata).
#' @param family_test character; must equal "binomial_bulk".
#' @param covariates character vector of covariate column names.
#' @param key list / row carrying MARKER/FIGURE/AREA/SUBAREA.
#' @param transformation_y character; passed through to the output.
#' @param session_folder character; ignored (kept for caller-symmetry).
#' @param independent_variable character; single column name (factor IV).
#' @param samples_sql_condition character; passed to io_data_preparation().
#' @param ... ignored.
#'
#' @return data.frame with one row per kept probe, schema as above.
#'   Returns NULL when Rfast is missing, the design can't be built, or
#'   no probes survive the degenerate-burden filter in data_preparation.
#'
#' @keywords internal
#' @noRd
assoc_glm_model_bulk <- function(tempDataFrame, g_start, family_test,
                            covariates = NULL, key,
                            transformation_y,
                            session_folder,
                            independent_variable,
                            samples_sql_condition,
                            inference_detail = NULL, ...) {

  if (!requireNamespace("Rfast", quietly = TRUE)) {
    stop("family_test='", family_test,
         "' requires the 'Rfast' package. Install it with:\n",
         "  install.packages('Rfast')")
  }

  if (family_test != "binomial_bulk") {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_glm_model_bulk: unexpected family_test='", family_test,
              "' (only 'binomial_bulk' supported).")
    return(NULL)
  }

  # 1. io_data_preparation: factors the IV (assoc_is_family_dicotomic branch),
  # then runs the AI-044 universal degenerate-burden filter - so by the
  # time we get back, tempDataFrame only contains informative probes.
  # inference_detail now arrives as an argument. It used to be looked up with
  # exists("inference_detail", inherits = TRUE), which from inside a function
  # searches the LEXICAL chain (namespace, imports, base, then the user's
  # global environment) and never the caller's frame. So the transformation_x
  # asked for in the request was silently dropped and the bulk fit ran on
  # untransformed data, while a stray object of that name in the workspace
  # would have been picked up instead.
  transformation_x_local <- if (!is.null(inference_detail$transformation_x))
    as.character(inference_detail$transformation_x) else "none"
  prepared <- io_data_preparation(family_test, transformation_y, tempDataFrame,
                                independent_variable, g_start, ncol(tempDataFrame),
                                covariates, key,
                                transformation_x = transformation_x_local)
  tempDataFrame <- prepared$tempDataFrame
  iv_levels <- prepared$independent_variableLevels

  cols <- colnames(tempDataFrame)
  g_end <- length(cols)
  if (g_start > g_end) {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_glm_model_bulk: no probes survived io_data_preparation - returning NULL.")
    return(NULL)
  }
  probe_cols <- cols[g_start:g_end]

  # 2. Drop samples with NA in IV or covariates (Rfast::glm_logistic
  # doesn't tolerate NAs in the design / response).
  use_iv <- c(independent_variable, covariates)
  use_iv <- use_iv[nzchar(use_iv) & use_iv %in% colnames(tempDataFrame)]
  keep_rows <- stats::complete.cases(tempDataFrame[, use_iv, drop = FALSE])
  td <- tempDataFrame[keep_rows, , drop = FALSE]
  if (nrow(td) < 5L) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_glm_model_bulk: too few complete samples (", nrow(td), " < 5) - skip.")
    return(NULL)
  }

  # 3. Design matrix: factor-expanded IV + covariates (no intercept here -
  # Rfast::glm_logistic adds it). One reference level dropped so each
  # remaining IV coefficient is a log-OR vs reference.
  iv_factor <- as.factor(td[, independent_variable])
  iv_factor <- droplevels(iv_factor)
  if (nlevels(iv_factor) < 2L) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_glm_model_bulk: IV has < 2 levels after droplevels - skip.")
    return(NULL)
  }
  iv_dummies <- stats::model.matrix(~iv_factor)[, -1L, drop = FALSE]
  iv_coef_names <- paste0(independent_variable, levels(iv_factor)[-1L])
  colnames(iv_dummies) <- iv_coef_names

  cov_used <- character(0)
  if (length(covariates) > 0L && any(nzchar(covariates))) {
    cov_used <- intersect(covariates[nzchar(covariates)], colnames(td))
    if (length(cov_used) > 0L) {
      cov_mat <- as.matrix(vapply(td[, cov_used, drop = FALSE], as.numeric, numeric(nrow(td))))
      if (!is.matrix(cov_mat)) cov_mat <- matrix(cov_mat, ncol = length(cov_used))
      colnames(cov_mat) <- cov_used
      design_no_int <- cbind(iv_dummies, cov_mat)
    } else {
      design_no_int <- iv_dummies
    }
  } else {
    design_no_int <- iv_dummies
  }
  storage.mode(design_no_int) <- "numeric"

  # 4. Response matrix:  probes x samples → as integer 0/1
  y_mat <- as.matrix(vapply(td[, probe_cols, drop = FALSE], as.integer, integer(nrow(td))))
  if (!is.matrix(y_mat)) y_mat <- matrix(y_mat, ncol = length(probe_cols))
  colnames(y_mat) <- probe_cols
  y_mat[is.na(y_mat)] <- 0L

  n_probes <- length(probe_cols)
  ncoef <- ncol(design_no_int) + 1L  # +1 = (Intercept)
  coef_names_full <- c("(Intercept)", colnames(design_no_int))

  core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
            " assoc_glm_model_bulk: fitting ", n_probes, "  probes x ",
            ncoef, " coefs via Rfast::glm_logistic (parallel foreach).")

  # 5. Per-probe fit in parallel. Each worker only needs Rfast +
  # the design matrix + one Y column → small payload.
  # Rfast::glm_logistic returns only estimates ($be) - NO standard errors.
  # We compute SEs from the Fisher information at the MLE:
  #   I(β) = X' diag(p(1-p)) X  →  Var(β̂) = I(β̂)^{-1}
  # This is the standard logistic SE (same formula stats::glm uses).
  X_full <- cbind(`(Intercept)` = 1, design_no_int)
  X_design <- design_no_int  # captured by foreach closure
  ncoef_local <- ncoef        # captured by foreach closure
  # Per-probe output width: ncoef estimates + ncoef p-values + 4 metrics
  # (MCFADDEN_R2, NAGELKERKE_R2, C_STATISTIC_AUC, DEVIANCE_RATIO).
  n_metrics <- 4L
  j <- NULL                   # quiet R CMD check note
  fits <- foreach::foreach(
    j = seq_len(n_probes),
    .combine  = rbind,
    .packages = c("Rfast", "stats"),
    .export   = c("X_design", "X_full", "y_mat", "ncoef_local", "n_metrics")
  ) %dorng% {
    na_vec <- c(rep(NA_real_, ncoef_local), rep(NA_real_, ncoef_local),
                rep(NA_real_, n_metrics))
    y <- y_mat[, j]
    # Skip degenerate Y (safety net - io_data_preparation should have caught it).
    if (length(unique(y)) < 2L) return(na_vec)
    f <- tryCatch(
      Rfast::glm_logistic(x = X_design, y = y),
      error = function(e) NULL
    )
    if (is.null(f) || is.null(f$be)) return(na_vec)
    est <- as.numeric(f$be)
    # Fisher info → SE
    eta <- as.numeric(X_full %*% est)
    p   <- 1 / (1 + exp(-eta))
    w   <- p * (1 - p)
    # Guard against numeric blow-up (separation): w == 0 for any row breaks I
    if (any(!is.finite(w)) || any(w <= 0)) {
      return(c(est, rep(NA_real_, ncoef_local), rep(NA_real_, n_metrics)))
    }
    XtWX <- crossprod(X_full * sqrt(w))
    vcov_mat <- tryCatch(solve(XtWX), error = function(e) NULL)
    if (is.null(vcov_mat)) {
      return(c(est, rep(NA_real_, ncoef_local), rep(NA_real_, n_metrics)))
    }
    se   <- sqrt(diag(vcov_mat))
    z    <- est / se
    pval <- 2 * stats::pnorm(-abs(z))

    # AI-044 (2026-06-09): goodness-of-fit metrics per probe. Registered
    # in metrics_properties.rda. Rationale: R²/R²_adj don't apply to
    # logistic - we report McFadden + Nagelkerke pseudo-R² (variance
    # explained analogs), C-statistic (= AUC, discrimination), and the
    # deviance ratio (devi/null_devi, lower = better fit).
    n <- length(y)
    p1 <- mean(y)
    null_devi <- if (p1 > 0 && p1 < 1) {
      -2 * (sum(y) * log(p1) + (n - sum(y)) * log(1 - p1))
    } else NA_real_
    metrics_vec <- c(NA_real_, NA_real_, NA_real_, NA_real_)
    if (is.finite(null_devi) && null_devi > 0 && !is.null(f$devi)) {
      devi <- f$devi
      mcfadden <- 1 - (devi / null_devi)
      cox_snell <- 1 - exp((devi - null_devi) / n)
      max_cs   <- 1 - exp(-null_devi / n)
      nagelkerke <- if (max_cs > 0) cox_snell / max_cs else NA_real_
      # C-stat = AUC via Mann-Whitney U
      n1 <- sum(y == 1L); n0 <- n - n1
      auc <- if (n1 > 0L && n0 > 0L) {
        rk <- rank(p)
        (sum(rk[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
      } else NA_real_
      metrics_vec <- c(mcfadden, nagelkerke, auc, devi / null_devi)
    }
    c(est, pval, metrics_vec)
  }

  if (is.null(fits) || nrow(fits) == 0L) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_glm_model_bulk: all fits failed - returning NULL.")
    return(NULL)
  }

  est_mat   <- fits[, seq_len(ncoef),                drop = FALSE]
  pval_mat  <- fits[, (ncoef + 1):(2 * ncoef),       drop = FALSE]
  # AI-044 (2026-06-09): goodness-of-fit metrics block (4 columns) sits
  # after the est/pval blocks. Order MUST match the c(est, pval, metrics_vec)
  # return in fit_one above - MCFADDEN_R2, NAGELKERKE_R2, C_STATISTIC_AUC,
  # DEVIANCE_RATIO. See metrics_properties.rda for direction.
  metrics_mat <- fits[, (2 * ncoef + 1):(2 * ncoef + 4), drop = FALSE]
  colnames(est_mat)     <- coef_names_full
  colnames(pval_mat)    <- coef_names_full
  colnames(metrics_mat) <- c("MCFADDEN_R2", "NAGELKERKE_R2",
                              "C_STATISTIC_AUC", "DEVIANCE_RATIO")

  # 6. Build result data.frame in the legacy schema. Coefficient column
  # names are sanitised the same way the per-probe path does it (toupper
  # + core_name_cleaning), so any downstream selector / FDR machinery sees
  # the same column shape it would have seen with stats::glm.
  cov_label <- if (length(cov_used) > 0L) paste(cov_used, collapse = " ") else NA

  result <- data.frame(
    INDIPENDENT_VARIABLE = rep(independent_variable, n_probes),
    MARKER               = rep(as.character(key$MARKER),  n_probes),
    FIGURE               = rep(as.character(key$FIGURE),  n_probes),
    AREA                 = rep(as.character(key$AREA),    n_probes),
    SUBAREA              = rep(as.character(key$SUBAREA), n_probes),
    AREA_OF_TEST         = probe_cols,
    FAMILY_TEST          = family_test,
    transformation_y     = transformation_y,
    COVARIATES           = cov_label,
    R_MODEL              = "Rfast::glm_logistic",
    stringsAsFactors     = FALSE
  )

  for (i in seq_along(coef_names_full)) {
    cn <- coef_names_full[i]
    pname <- core_name_cleaning(paste0(cn, "_PVALUE"))
    ename <- core_name_cleaning(paste0(cn, "_ESTIMATE"))
    result[[pname]] <- pval_mat[, i]
    result[[ename]] <- est_mat[, i]
  }

  # AI-044 (2026-06-09): goodness-of-fit metrics block - names canonical
  # (uppercase, registered in metrics_properties.rda).
  result$MCFADDEN_R2     <- metrics_mat[, "MCFADDEN_R2"]
  result$NAGELKERKE_R2   <- metrics_mat[, "NAGELKERKE_R2"]
  result$C_STATISTIC_AUC <- metrics_mat[, "C_STATISTIC_AUC"]
  result$DEVIANCE_RATIO  <- metrics_mat[, "DEVIANCE_RATIO"]

  # Top-level PVALUE = first non-intercept coefficient p-value, mirroring
  # the existing pattern (apply_stat_model_batch.R uses first poly term).
  if (ncoef >= 2L) {
    result$PVALUE <- pval_mat[, 2L]
  }

  colnames(result) <- toupper(colnames(result))
  colnames(result) <- core_name_cleaning(colnames(result))

  # AI-257: no adjustment here. This function sees the probes of one batch, and
  # a family made of whatever happened to be in a batch is a memory parameter.
  # The three levels - key, scope, file - are computed in
  # assoc_analysis_save_results(), the one place where every row of a family is
  # together, and each is named after the family it controls.

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " assoc_glm_model_bulk: produced ", nrow(result), " probe-level rows.")
  result
}
