#' Limma-based polynomial association on a single genomic area.
#'
#' Mirror of \code{assoc_model_polynomial()} but the backend is
#' \code{limma::lmFit + limma::eBayes} instead of \code{stats::lm}.
#' Parses \code{family_test} as \code{limma_<degree>[_<partition>]} -
#' partition is optional and defaults to 1 because empirical-bayes
#' shrinkage from \code{eBayes()} replaces the train/test split that
#' polynomial uses.
#'
#' Per-area mode caveat: this function is called once per area by
#' \code{assoc_apply_stat_model()}, so lmFit gets a 1-row response matrix and
#' eBayes shrinkage degenerates to ordinary OLS t-stat. The full
#' cross-area shrinkage advantage of limma only materialises in batch
#' mode (planned in a later phase of AI-040).
#'
#' Dispatcher in \code{execute_model.R} checks
#' \code{requireNamespace("limma")} before calling this - see AI-038's
#' dispatch=guard pattern.
#'
#' @keywords internal
#' @noRd
assoc_model_limma <- function(family_test, tempDataFrame, sig.formula,
                                     transformation_y, plot,
                                     samples_sql_condition = samples_sql_condition,
                                     key) {

  # Note: unlike assoc_model_polynomial() this function does not
  # need ssEnv - there's no plotting branch yet and no folder lookup.
  # Adding core_get_session_info() here breaks unit tests that call the
  # model in isolation without a materialised session.

  limma_params <- unlist(strsplit(as.character(family_test), "_"))
  if (length(limma_params) < 2 || length(limma_params) > 3) {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
              " family_test='", family_test,
              "' is malformed. Use limma_<degree> or limma_<degree>_<partition>,",
              " e.g. limma_2 or limma_2_1.")
    return(data.frame())
  }

  degree <- suppressWarnings(as.numeric(limma_params[2]))
  if (is.na(degree) || degree < 1 || degree != as.integer(degree)) {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
              " family_test='", family_test,
              "': degree must be a positive integer.")
    return(data.frame())
  }
  degree <- as.integer(degree)
  partition_percentage <- if (length(limma_params) == 3)
    suppressWarnings(as.numeric(limma_params[3])) else 1
  if (is.na(partition_percentage)) partition_percentage <- 1

  res <- data.frame(PL_DEGREE = degree, PL_PERC = partition_percentage)
  res$r_model <- "limma::lmFit+eBayes"

  tempDataFrame <- as.data.frame(tempDataFrame)
  dep_var <- assoc_sig_formula_vars(sig.formula)
  dependent_variable <- dep_var$dependent_variable
  independent_variable <- dep_var$independent_variable
  covariates <- dep_var$covariates

  if (nrow(tempDataFrame) == 0) return(res)

  # Drop rows with NA in any column used by the model - limma::lmFit
  # accepts NAs but the degenerate 1-row matrix doesn't benefit from
  # observation weights, so we just exclude.
  use_cols <- c(dependent_variable, independent_variable, covariates)
  use_cols <- use_cols[nzchar(use_cols)]
  use_cols <- intersect(use_cols, colnames(tempDataFrame))
  td <- tempDataFrame[stats::complete.cases(tempDataFrame[, use_cols, drop = FALSE]),
                       , drop = FALSE]
  if (nrow(td) < (degree + length(covariates) + 2L)) return(res)

  iv_vec <- as.numeric(td[, independent_variable])
  poly_mat <- stats::poly(iv_vec, degree, raw = TRUE)
  colnames(poly_mat) <- paste0("STATS_POLY_EVAL_PARSE_TEXT_EQ_",
                                independent_variable, "_EQ_RAW_EQ_TRUE_",
                                seq_len(degree))

  if (length(covariates) > 0L && any(nzchar(covariates))) {
    cov_mat <- as.matrix(td[, covariates[nzchar(covariates)], drop = FALSE])
    design <- cbind(`(Intercept)` = 1, poly_mat, cov_mat)
  } else {
    design <- cbind(`(Intercept)` = 1, poly_mat)
  }

  y_vec <- as.numeric(td[, dependent_variable])
  y_mat <- matrix(y_vec, nrow = 1L,
                   dimnames = list(as.character(key$AREA_OF_TEST %||% "area"),
                                    NULL))

  fit <- tryCatch(limma::lmFit(y_mat, design),
                   error = function(e) NULL)
  if (is.null(fit)) return(res)
  fit <- tryCatch(limma::eBayes(fit),
                   error = function(e) NULL)
  if (is.null(fit)) return(res)

  # Coefficient-wise p-values and estimates. This was written to produce "the
  # same column-naming scheme as assoc_model_polynomial", and it never did:
  # measured 2026-10-03, this path gives STATS_POLY_..._<IV>_EQ_1_PVALUE and that
  # one gave STATS_POLY_..._<IV>_DEGREE_1_PVALUE. The two have never matched, so
  # the claim was a comment and not a contract.
  #
  # The polynomial path has since moved to I_<IV>_<k>, which is what the formula
  # it builds cleans to on its own. Aligning this path to it means changing the
  # names limma writes, which is a decision about this family's results and not
  # a side effect of fixing another one, so it is tracked rather than done here.
  # Until then the names of the two families differ, as they always have.
  coef_names <- colnames(fit$coefficients)
  for (i in seq_along(coef_names)) {
    row_name <- coef_names[i]
    pval_name <- core_name_cleaning(paste0(row_name, "_pvalue"))
    pval_name <- core_name_cleaning(gsub("_STATS_POLY_EVAL_PARSE_TEXT_EQ", "", pval_name))
    pval_name <- core_name_cleaning(gsub("_RAW_EQ_TRUE", "", pval_name))
    pval_name <- core_name_cleaning(gsub("INDEPENDENT_VARIABLE", independent_variable, pval_name))
    p_value <- data.frame(p_value = fit$p.value[1L, i])
    colnames(p_value) <- pval_name
    res <- cbind(res, p_value)

    estimate_name <- core_name_cleaning(paste0(row_name, "_estimate"))
    estimate <- data.frame(estimate = fit$coefficients[1L, i])
    colnames(estimate) <- core_name_cleaning(estimate_name)
    res <- cbind(res, estimate)
  }

  rownames(res) <- NULL
  res
}

# Local null-coalescing helper. Avoids importing rlang just for this.
`%||%` <- function(a, b) if (is.null(a)) b else a
