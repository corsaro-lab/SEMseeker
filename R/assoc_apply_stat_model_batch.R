# batch path for limma_<N> and voom_<N>.
#
# Why a separate batch path:
#   - Per-area limma (one area at a time) feeds lmFit a 1-row matrix, so
#     eBayes shrinkage collapses to OLS t-stat - same answer as
#     stats::lm(), no statistical gain.
#   - voom literally cannot work per-area: the precision weights come
#     from the empirical mean-variance trend estimated ACROSS areas. A
#     1-row matrix has no trend to fit.
# Batch mode calls lmFit (or voom + lmFit) once on the whole M x N
# response matrix (M genes / areas, N samples), then eBayes, then
# extracts per-area results in the same schema assoc_apply_stat_model would
# have returned area-by-area.
#
# Architectural placement: apply_stat_model.R early-returns into this
# function when family_test matches '^(limma|voom)_'. The per-area
# assoc_model_limma() from Fase 1 stays in place for the direct
# assoc_execute_model() entry point used by unit tests; it is NOT reached
# during a normal pipeline run because the dispatch happens upstream
# in apply_stat_model.R.

#' Fit a batch limma / voom model across all genomic areas at once.
#'
#' @param tempDataFrame data.frame after io_data_preparation(); rows are
#'   samples, columns are: independent variable + covariates + one
#'   column per genomic area (the burden for that area).
#' @param g_start integer index of the first area column (everything
#'   before it is sample-level metadata: IV + covariates).
#' @param family_test character, of the form 'limma_<degree>' or
#'   'voom_<degree>' (optional '_<partition>' suffix is ignored -
#'   eBayes shrinkage replaces train/test holdout).
#' @param covariates character vector of covariate column names.
#' @param key data.frame row carrying MARKER/FIGURE/AREA/SUBAREA.
#' @param transformation_y character label for the y-side transformation
#'   (currently passed through unchanged in the output row).
#' @param independent_variable character, single column name.
#' @param ... ignored (kept for caller-symmetry with assoc_apply_stat_model).
#'
#' @return data.frame with one row per kept area. Returns NULL when
#'   the design matrix can't be built or limma::lmFit fails on every
#'   area.
#'
#' @keywords internal
#' @noRd
assoc_apply_stat_model_batch <- function(tempDataFrame, g_start, family_test,
                                    covariates = NULL, key,
                                    transformation_y,
                                    session_folder,
                                    independent_variable,
                                    samples_sql_condition,
                                    inference_detail = NULL, ...) {

  # Parser: <engine>_<degree>[_<partition>]
  parts  <- unlist(strsplit(as.character(family_test), "_"))
  engine <- parts[1]
  if (!engine %in% c("limma", "voom") || length(parts) < 2L) {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_apply_stat_model_batch: malformed family_test='", family_test, "'")
    return(NULL)
  }
  degree <- suppressWarnings(as.integer(parts[2]))
  if (is.na(degree) || degree < 1L) {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_apply_stat_model_batch: family_test='", family_test,
              "' has invalid degree.")
    return(NULL)
  }
  partition_percentage <- if (length(parts) >= 3L)
    suppressWarnings(as.numeric(parts[3])) else 1
  if (is.na(partition_percentage)) partition_percentage <- 1

  # Same pre-processing the per-area path runs
  transformation_x_local <- if (!is.null(inference_detail$transformation_x)) as.character(inference_detail$transformation_x) else "none"
  level_order_local <- if (!is.null(inference_detail$independent_variable_order))
    as.character(inference_detail$independent_variable_order) else NULL
  prepared <- io_data_preparation(family_test, transformation_y, tempDataFrame,
                                independent_variable, g_start, ncol(tempDataFrame),
                                covariates, key,
                                transformation_x = transformation_x_local,
                                independent_variable_order = level_order_local)
  tempDataFrame <- prepared$tempDataFrame

  # The burden columns by name, as io_data_preparation() found and checked them.
  area_cols <- prepared$burden_columns
  if (length(area_cols) == 0L) return(NULL)

  # Drop areas with no variance - uninformative, would crash lmFit.
  area_keep <- vapply(area_cols, function(a) {
    v <- tempDataFrame[, a]
    length(unique(v[!is.na(v)])) >= 2L
  }, logical(1))
  area_cols <- area_cols[area_keep]
  if (length(area_cols) == 0L) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_apply_stat_model_batch: no informative areas after constant-burden filter.")
    return(NULL)
  }

  # Drop samples with NA in IV or covariates
  use_iv <- c(independent_variable, covariates)
  use_iv <- use_iv[nzchar(use_iv) & use_iv %in% colnames(tempDataFrame)]
  keep_rows <- stats::complete.cases(tempDataFrame[, use_iv, drop = FALSE])
  td <- tempDataFrame[keep_rows, , drop = FALSE]
  min_n <- degree + length(covariates) + 2L
  if (nrow(td) < min_n) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_apply_stat_model_batch: too few complete samples (",
              nrow(td), " < ", min_n, ") - skip.")
    return(NULL)
  }

  # Design matrix: poly(IV, degree) + covariates. Use 'I_<IV>_<n>'
  # naming to match the assoc_model_polynomial CSV schema.
  iv_vec <- as.numeric(td[, independent_variable])
  poly_mat <- stats::poly(iv_vec, degree, raw = TRUE)
  colnames(poly_mat) <- paste0("I_", independent_variable, "_", seq_len(degree))

  cov_used <- character(0)
  if (length(covariates) > 0L && any(nzchar(covariates))) {
    cov_used <- intersect(covariates[nzchar(covariates)], colnames(td))
    if (length(cov_used) > 0L) {
      cov_mat <- as.matrix(td[, cov_used, drop = FALSE])
      design  <- cbind(`(Intercept)` = 1, poly_mat, cov_mat)
    } else {
      design <- cbind(`(Intercept)` = 1, poly_mat)
    }
  } else {
    design <- cbind(`(Intercept)` = 1, poly_mat)
  }

  # Response matrix: M areas x N samples. An empty burden cell follows the
  # one rule for pivot values (util_absent_as_zero).
  y_mat <- t(as.matrix(td[, area_cols, drop = FALSE]))
  y_mat <- util_absent_as_zero(y_mat, key$MARKER)
  rownames(y_mat) <- area_cols
  mode(y_mat) <- "numeric"

  # Fit
  if (engine == "voom") {
    voom_obj <- tryCatch(limma::voom(y_mat, design),
                          error = function(e) {
                            core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
                                      " voom failed: ", conditionMessage(e))
                            NULL
                          })
    if (is.null(voom_obj)) return(NULL)
    fit <- limma::lmFit(voom_obj, design)
    r_model_label <- "limma::voom+lmFit+eBayes"
  } else {
    fit <- tryCatch(limma::lmFit(y_mat, design),
                     error = function(e) {
                       core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
                                 " lmFit failed: ", conditionMessage(e))
                       NULL
                     })
    if (is.null(fit)) return(NULL)
    r_model_label <- "limma::lmFit+eBayes"
  }
  fit <- limma::eBayes(fit)

  # Extract per-area results matching the schema assoc_apply_stat_model
  # would have returned area-by-area through assoc_model_polynomial:
  #   MARKER, FIGURE, AREA, SUBAREA, AREA_OF_TEST,
  #   PL_DEGREE, PL_PERC, R_MODEL,
  #   FAMILY_TEST, TRANSFORMATION_Y, INDEPENDENT_VARIABLE, COVARIATES,
  #   then per-coefficient PVALUE + ESTIMATE columns with the same
  #   name-cleaning the polynomial path uses, finally a top-level
  #   PVALUE column equal to the first poly-term coefficient pvalue
  #   (so the existing FDR + selector machinery in assoc_apply_stat_model's
  #   caller works without changes).
  coef_names <- colnames(fit$coefficients)

  build_pname <- function(cn) {
    pn <- core_name_cleaning(paste0(cn, "_pvalue"))
    pn <- core_name_cleaning(gsub("_STATS_POLY_EVAL_PARSE_TEXT_EQ", "", pn))
    pn <- core_name_cleaning(gsub("_RAW_EQ_TRUE", "", pn))
    pn <- core_name_cleaning(gsub("INDEPENDENT_VARIABLE", independent_variable, pn))
    pn
  }
  build_ename <- function(cn) {
    en <- core_name_cleaning(paste0(cn, "_estimate"))
    en <- core_name_cleaning(gsub("_STATS_POLY_EVAL_PARSE_TEXT_EQ", "", en))
    en <- core_name_cleaning(gsub("_RAW_EQ_TRUE", "", en))
    en <- core_name_cleaning(gsub("INDEPENDENT_VARIABLE", independent_variable, en))
    en
  }
  pnames <- vapply(coef_names, build_pname, character(1))
  enames <- vapply(coef_names, build_ename, character(1))

  cov_label <- if (length(cov_used) > 0L) paste(cov_used, collapse = "+") else ""

  result_temp <- data.frame(
    MARKER               = rep(as.character(key$MARKER),  length(area_cols)),
    FIGURE               = rep(as.character(key$FIGURE),  length(area_cols)),
    AREA                 = rep(as.character(key$AREA),    length(area_cols)),
    SUBAREA              = rep(as.character(key$SUBAREA), length(area_cols)),
    AREA_OF_TEST         = area_cols,
    PL_DEGREE            = degree,
    PL_PERC              = partition_percentage,
    R_MODEL              = r_model_label,
    FAMILY_TEST          = as.character(family_test),
    TRANSFORMATION_Y     = as.character(transformation_y),
    INDEPENDENT_VARIABLE = as.character(independent_variable),
    COVARIATES           = cov_label,
    stringsAsFactors     = FALSE
  )
  for (i in seq_along(coef_names)) {
    result_temp[[pnames[i]]] <- fit$p.value[, i]
    result_temp[[enames[i]]] <- fit$coefficients[, i]
  }

  # 2026-06-09: emit eBayes-derived goodness-of-fit metrics
  # registered in metrics_properties.rda. These replace R²/R²_adj for the
  # limma/voom batch path (R² isn't natively returned by lmFit and is
  # not the canonical limma diagnostic anyway).
  #   T_STAT_MODERATED: moderated t-stat for the first non-intercept coef
  #     (per-area; positive = positive effect direction). Higher |t| = stronger
  #     evidence per row.
  #   B_STATISTIC: log-odds posterior of differential expression (limma's lods),
  #     same coef. Pseudo-R² analog: Higher = stronger evidence.
  #   F_STAT_MODERATED: moderated F for the joint contrast of ALL non-intercept
  #     coefs (e.g. poly_1 + poly_2 → "parabolic effect present?"). Per row.
  #   POSTERIOR_RESIDUAL_VAR: s2.post - posterior residual variance per area.
  #     Lower = better fit; useful as a diagnostic alongside the t-stat.
  if (ncol(fit$coefficients) >= 2L) {
    result_temp$T_STAT_MODERATED <- as.numeric(fit$t[, 2L])
    if (!is.null(fit$lods)) {
      result_temp$B_STATISTIC <- as.numeric(fit$lods[, 2L])
    }
  }
  # Joint F-statistic across ALL non-intercept coefs (i.e. excluding intercept).
  # Computed lazily - only when degree >= 2 since for degree 1 it equals t².
  if (ncol(fit$coefficients) >= 3L) {
    coef_idx_no_intercept <- seq_len(ncol(fit$coefficients))[-1L]
    f_obj <- tryCatch(
      limma::topTable(fit, number = Inf, sort.by = "none"),
      error = function(e) NULL
    )
    if (!is.null(f_obj) && "F" %in% colnames(f_obj)) {
      f_obj <- f_obj[rownames(fit$coefficients), , drop = FALSE]
      result_temp$F_STAT_MODERATED <- as.numeric(f_obj[, "F"])
    }
  }
  if (!is.null(fit$s2.post)) {
    result_temp$POSTERIOR_RESIDUAL_VAR <- as.numeric(fit$s2.post)
  }

  # Top-level PVALUE = first non-intercept (= first poly term) pvalue,
  # so the BH adjustment + significativity selector in the assoc_apply_stat_model
  # caller path can hook into it the same way it does for polynomial.
  first_poly_pcol <- pnames[2L]  # coef 1 methods::is (Intercept), coef 2 is poly_1
  if (!is.null(first_poly_pcol) && first_poly_pcol %in% colnames(result_temp)) {
    result_temp$PVALUE <- result_temp[[first_poly_pcol]]
  }

  colnames(result_temp) <- core_name_cleaning(colnames(result_temp))
  result_temp
}
