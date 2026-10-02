# AI-061: lazy / low-memory batch path for limma_<N> and voom_<N>.
#
# Why a separate function: apply_stat_model_batch.R receives a fully
# materialised, transposed, sample-merged tempDataFrame that the caller
# (assoc_run_marker) builds with the standard
#   tempDataFrame <- as.data.frame(pivot_lazy$collect())
#   ... t() ... merge(sample_sheet, ...)
# pipeline. Even at GENE/WHOLE scale (~17 k genes × ~500 samples
# ≈ 70 MB raw) that path peaks around 4×-5× the raw size in R memory
# alone, plus the polars Rust DataFrame that stays alive in parallel,
# plus Arrow IPC buffers from the collect step.  Smoke SIGNAL@PROBE
# (366 k × 500 ≈ 1.4 GB raw) ran to 45 GB before lmFit was even
# called and was OOM-killed at ~60 GB.
#
# This function takes the pivot's polars LazyFrame directly, applies
# the AI-043 area_to_remove filter LAZILY, materialises only ONE
# R matrix (genes × samples) without an intervening data.frame, and
# explicitly rm()+gc()'s the polars-side artifacts before lmFit so
# the peak working set drops to roughly:
#   polars collect      ≈ 1×  raw  (in Rust heap, drops after rm+gc)
#   R y_mat             ≈ 1×  raw  (the single matrix lmFit sees)
#   design + lmFit work  small
# i.e. ~3×-4× the raw size at peak vs ~30× on the legacy path.

#' Polars-lazy batch fit for limma_<N> and voom_<N>.
#'
#' @param pivot_lazy polars LazyFrame. The annotated AREA-level pivot
#'   for one (MARKER, FIGURE, AREA, SUBAREA) key, with at least an
#'   AREA column and one column per sample.
#' @param sample_sheet data.frame. Per-sample metadata returned by
#'   sem_prepare_study_for_analysis(); must contain Sample_ID, the
#'   independent variable, and any covariates referenced.
#' @param family_test character, of the form 'limma_<degree>' or
#'   'voom_<degree>' (optional '_<partition>' suffix is ignored).
#' @param covariates character vector of covariate column names.
#' @param key data.frame row carrying MARKER/FIGURE/AREA/SUBAREA.
#' @param transformation_y character label passed through unchanged.
#' @param independent_variable character, single column name.
#' @param area_to_remove character. AREA values to skip (resume cache).
#' @param ... ignored (kept for caller-symmetry).
#'
#' @return data.frame with one row per kept area, same schema
#'   assoc_apply_stat_model_batch() returns. NULL on degenerate input.
#'
#' @keywords internal
#' @noRd
assoc_apply_stat_model_batch_lazy <- function(pivot_lazy,
                                         sample_sheet,
                                         family_test,
                                         covariates = NULL,
                                         key,
                                         transformation_y,
                                         independent_variable,
                                         area_to_remove = character(0),
                                         ...) {

  # Parser: <engine>_<degree>[_<partition>]
  parts  <- unlist(strsplit(as.character(family_test), "_"))
  engine <- parts[1]
  if (!engine %in% c("limma", "voom") || length(parts) < 2L) {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_apply_stat_model_batch_lazy: malformed family_test='",
              family_test, "'")
    return(NULL)
  }
  degree <- suppressWarnings(as.integer(parts[2]))
  if (is.na(degree) || degree < 1L) {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_apply_stat_model_batch_lazy: invalid degree in '",
              family_test, "'")
    return(NULL)
  }
  partition_percentage <- if (length(parts) >= 3L)
    suppressWarnings(as.numeric(parts[3])) else 1
  if (is.na(partition_percentage)) partition_percentage <- 1

  # Drop genomic-coord and probe-metadata columns LAZILY - no
  # materialisation yet.
  schema_names <- names(pivot_lazy$collect_schema())
  drop_cols <- intersect(schema_names,
                          c("CHR", "START", "END", "PROBE",
                            "K27", "K450", "K850"))
  if (length(drop_cols) > 0L) {
    pivot_lazy <- pivot_lazy$drop(drop_cols)
  }

  # AI-043 resume: filter out areas already in the on-disk CSV.
  # AI-061+ (2026-06-09): no more gsub("-","_") normalisation on either
  # side - names stay pass-through from the upstream annotation. The
  # downstream CSV and the pivot AREA column carry identical raw names,
  # so $is_in() matches exactly.
  # NB: $is_in() must receive a polars Expression / Series, NOT a bare R
  # character vector - otherwise polars 1.x parses each string as a column
  # reference and fails with "Column(s) not found: '<first value>' not found".
  # Wrap via pl$lit()$implode() so the values are treated as a literal set.
  if (length(area_to_remove) > 0L) {
    pivot_lazy <- pivot_lazy$filter(
      !polars::pl$col("AREA")$is_in(
        polars::pl$lit(area_to_remove)$implode()
      )
    )
  }

  # AI-044 / AI-061 (2026-06-09): apply transformation_y + universal
  # degenerate-burden filter LAZILY before materialisation. Before this
  # change the lazy path silently skipped io_data_preparation() entirely,
  # so any `transformation_y` ≠ "none" on a limma_/voom_ inference_detail
  # produced a CSV with UN-transformed values (silent bug), and rows with
  # var(Y) == 0 made it through to lmFit producing NaN t-stats. See
  # `io_data_preparation_lazy()` for the polars-native equivalent of the
  # R-side `io_data_preparation()` Y-side transformations + AI-044 filter.
  # Unification of the two paths is tracked as AI-097 in the backlog.
  schema_pre <- names(pivot_lazy$collect_schema())
  sample_cols_pre <- setdiff(schema_pre, c("AREA", "PROBE", "CHR", "START", "END",
                                            "K27", "K450", "K850"))
  pivot_lazy <- io_data_preparation_lazy(
    pivot_lazy        = pivot_lazy,
    sample_cols       = sample_cols_pre,
    transformation_y  = transformation_y,
    apply_degenerate_filter = TRUE,
    key               = key,
    family_test       = family_test
  )

  # AI-061+ (2026-06-09): SEPARATE LAZY PREP FROM FIT.
  # We need n_genes + sample_cols BEFORE materialising y_mat so the
  # memory gate can decide monolithic vs chunked. Both pieces are
  # discoverable lazily - schema gives us sample columns, $select($len)
  # gives us a row count without ever pulling values into R.
  schema_post <- names(pivot_lazy$collect_schema())
  sample_cols_all <- setdiff(schema_post, c("AREA", "PROBE", "CHR", "START", "END",
                                              "K27", "K450", "K850"))
  n_genes <- as.integer(as.data.frame(
    pivot_lazy$select(polars::pl$len()$alias("n"))$collect())$n[1])
  if (is.na(n_genes) || n_genes == 0L) {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_apply_stat_model_batch_lazy: nothing left after resume",
              " filter for ", key$MARKER, "/", key$FIGURE, "/",
              key$AREA, "/", key$SUBAREA, ".")
    return(NULL)
  }

  # Align sample_sheet to the lazy pivot's sample columns. Drop samples
  # without a matching row OR with NA in IV/covariates. This is the
  # same logic that used to operate on the materialised y_mat - pulled
  # forward so we can compute the design BEFORE the memory gate.
  ss <- sample_sheet[match(sample_cols_all,
                            as.character(sample_sheet$Sample_ID)), , drop = FALSE]
  use_iv <- c(independent_variable, covariates)
  use_iv <- use_iv[nzchar(use_iv) & use_iv %in% colnames(ss)]
  keep <- !is.na(ss$Sample_ID) &
          stats::complete.cases(ss[, use_iv, drop = FALSE])
  if (sum(keep) < (degree + length(covariates) + 2L)) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_apply_stat_model_batch_lazy: too few complete samples (",
              sum(keep), " < ", degree + length(covariates) + 2L, ")")
    return(NULL)
  }
  ss <- ss[keep, , drop = FALSE]
  sample_cols_kept <- sample_cols_all[keep]

  # Design matrix: poly(IV, degree, raw=TRUE) + covariates. Name the
  # polynomial columns the same way assoc_model_polynomial's
  # I(...) formula winds up after core_name_cleaning - 'I_<IV>_<deg>' - so
  # the CSV columns landed by build_pname/build_ename match the
  # polynomial CSV schema bit-for-bit.
  iv_vec <- as.numeric(ss[, independent_variable])
  poly_mat <- stats::poly(iv_vec, degree, raw = TRUE)
  colnames(poly_mat) <- paste0("I_", independent_variable, "_", seq_len(degree))

  cov_used <- character(0)
  if (length(covariates) > 0L && any(nzchar(covariates))) {
    cov_used <- intersect(covariates[nzchar(covariates)], colnames(ss))
  }
  if (length(cov_used) > 0L) {
    cov_mat <- as.matrix(ss[, cov_used, drop = FALSE])
    design  <- cbind(`(Intercept)` = 1, poly_mat, cov_mat)
  } else {
    design <- cbind(`(Intercept)` = 1, poly_mat)
  }

  # ---- MEMORY GATE & DISPATCH ----
  # Decides if a monolithic lmFit fits in budget. Chunked path activates
  # when the monolithic y_mat would push past
  # `total_RAM × SEMSEEKER_BULK_MODEL_MEM_FRACTION` (default 0.6).
  ssEnv_local      <- tryCatch(core_get_session_info(), error = function(e) NULL)
  tech_is_longread <- !is.null(ssEnv_local$tech) &&
                       ssEnv_local$tech %in% c("WGBS", "LONGREAD")
  memgate <- .core_bulk_model_memory_gate(
    n_probes  = n_genes,
    n_samples = length(sample_cols_kept),
    n_coef    = ncol(design)
  )
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " assoc_apply_stat_model_batch_lazy ", family_test,
            " [", key$MARKER, "/", key$FIGURE, "/", key$AREA, "/",
            key$SUBAREA, "] gate decision=", memgate$decision,
            " (mono=", round(memgate$mono_peak_GB, 1),
            " GB | chunk=", round(memgate$chunk_peak_GB, 1),
            " GB | avail=", round(memgate$available_GB, 1), " GB)")

  fit <- switch(memgate$decision,
    "monolithic" = assoc_lmfit_monolithic_lazy(
      pivot_lazy        = pivot_lazy,
      sample_cols_kept  = sample_cols_kept,
      design            = design,
      engine            = engine,
      key               = key,
      family_test       = family_test
    ),
    "chunked"    = assoc_lmfit_chunked_by_chr(
      pivot_lazy        = pivot_lazy,
      sample_cols_kept  = sample_cols_kept,
      design            = design,
      engine            = engine,
      key               = key,
      family_test       = family_test,
      probe_features    = if (!tech_is_longread)
                            tryCatch(anno_probe_features_get("PROBE"),
                                      error = function(e) NULL) else NULL,
      tech_is_longread  = tech_is_longread
    ),
    "abort"      = {
      stop(sprintf(
        "assoc_apply_stat_model_batch_lazy: even the biggest chunk's lmFit exceeds budget on [%s/%s/%s/%s]. needed=%.1f GB (chunk peak), avail=%.1f GB. Raise SEMSEEKER_BULK_MODEL_MEM_FRACTION or move to a bigger machine.",
        key$MARKER, key$FIGURE, key$AREA, key$SUBAREA,
        memgate$chunk_peak_GB, memgate$available_GB
      ))
    }
  )
  if (is.null(fit)) return(NULL)

  r_model_label <- if (engine == "voom") "limma::voom+lmFit+eBayes"
                   else                  "limma::lmFit+eBayes"
  area_cols <- rownames(fit$coefficients)
  fit <- limma::eBayes(fit)

  # Build the result data.frame in the same schema assoc_apply_stat_model_batch
  # returns - same column names, so the FDR + selector machinery in
  # assoc_analysis_save_results() picks them up unchanged.
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

  first_poly_pcol <- pnames[2L]  # coef 1 methods::is (Intercept), coef 2 is poly_1
  if (!is.null(first_poly_pcol) && first_poly_pcol %in% colnames(result_temp)) {
    result_temp$PVALUE <- result_temp[[first_poly_pcol]]
  }

  colnames(result_temp) <- core_name_cleaning(colnames(result_temp))

  # Release the heavy locals BEFORE the function returns. The MArrayLM
  # carries several gene-sized matrices; without explicit cleanup the
  # next batch piles a fresh fit on top of the previous one (R's GC
  # is lazy), driving the process into Jetsam OOM by the second
  # batch (limma_2 SIGNAL@PROBE, 2026-06-05). y_mat is now released
  # inside `assoc_lmfit_monolithic_lazy()` / `assoc_lmfit_chunked_by_chr()`.
  rm(fit, design)
  if (exists("poly_mat", inherits = FALSE)) rm(poly_mat)
  if (exists("cov_mat",  inherits = FALSE)) rm(cov_mat)
  gc(verbose = FALSE)

  result_temp
}
