assoc_execute_model <- function(family_test, tempDataFrame, sig.formula, burdenValue, independent_variable, transformation_y="", plot=FALSE,
  samples_sql_condition="", key=c()){

  area <- as.character(key$AREA)
  subarea <- as.character(key$SUBAREA)
  marker <- as.character(key$MARKER)
  figure <- as.character(key$FIGURE)

  samples_sql_condition <- as.character(samples_sql_condition)
  transformation_y <- as.character(transformation_y)
  family_test <- as.character(family_test)
  # model_result <- data.frame()
  #
  # tryCatch({
  if( family_test=="binomial" | family_test=="poisson" | family_test=="gaussian")
    model_result <- assoc_glm_model(family_test, tempDataFrame, sig.formula, transformation_y, plot , samples_sql_condition=samples_sql_condition,key)

  if(family_test=="wilcoxon" | family_test=="t.test" | family_test=="pearson" | family_test=="kendall" | family_test=="spearman" | family_test=="jsd"
    | family_test=="chisq.test" | family_test=="fisher.test" | family_test=="kruskal.test" | family_test=="bartlett.test")
    model_result <- assoc_test_model(family_test =  family_test, tempDataFrame =  tempDataFrame, sig.formula =  sig.formula,
      burdenValue = burdenValue, independent_variable = independent_variable, transformation_y =  transformation_y, plot =  plot,
      samples_sql_condition = samples_sql_condition, key)

  if (grepl("wilcoxon.paired",family_test) | grepl("t.test.paired",family_test ))
    model_result <- assoc_test_model_paired(family_test, tempDataFrame, sig.formula,burdenValue,independent_variable, transformation_y, plot, samples_sql_condition=samples_sql_condition,key)

  if(grepl("mean-permutation",family_test))
    model_result <- assoc_mean_permutation(family_test, sig.formula, tempDataFrame, independent_variable,plot, samples_sql_condition=samples_sql_condition,key)

  if(grepl("quantile-permutation",family_test))
    model_result <- assoc_quantile_permutation_model(family_test, sig.formula, tempDataFrame, independent_variable, samples_sql_condition=samples_sql_condition,key)

  if(grepl("spearman-permutation",family_test))
    model_result <- assoc_spearman_permutation(family_test, sig.formula, tempDataFrame, independent_variable, samples_sql_condition=samples_sql_condition,key)

  if(grepl("polynomial",family_test))
    model_result <- assoc_model_polynomial(family_test, tempDataFrame, sig.formula, transformation_y, plot, samples_sql_condition=samples_sql_condition,key)

  # AI-040 Fase 1: limma_<degree>[_<partition>] backend (Suggests:limma).
  # Dispatcher = guard point (AI-038): fail fast with install hint if the
  # opt-in dependency is missing, before constructing any design matrix.
  # NOTE: in normal pipeline runs limma_<N> is intercepted upstream by
  # assoc_apply_stat_model() and routed through the batch path. This per-area
  # branch is kept for direct assoc_execute_model() calls (unit tests, custom
  # downstream callers).
  if (grepl("^limma_", family_test)) {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("family_test='", family_test,
           "' requires the 'limma' package. Install it with:\n",
           "  BiocManager::install('limma')")
    }
    model_result <- assoc_model_limma(
      family_test, tempDataFrame, sig.formula, transformation_y, plot,
      samples_sql_condition = samples_sql_condition, key)
  }

  # AI-040 Fase 2: voom_<degree> is BATCH-ONLY. Per-area voom can't
  # estimate its mean-variance trend (lowess on 1 point), so any per-area
  # entry point is an architectural mistake. Refuse explicitly instead
  # of returning a degenerate result.
  if (grepl("^voom_", family_test)) {
    stop("family_test='", family_test,
         "' is batch-only: voom estimates its mean-variance trend across ",
         "many genomic areas at once. Call assoc_apply_stat_model() instead of ",
         "assoc_execute_model() - assoc_apply_stat_model dispatches batch families ",
         "to assoc_apply_stat_model_batch() automatically.")
  }

  if(family_test=="multinomial")
    model_result <- assoc_multinomial_model(family_test, tempDataFrame, sig.formula, transformation_y, plot, samples_sql_condition=samples_sql_condition,key)

  if(grepl("log10_",family_test) | grepl("pow10_",family_test) | grepl("log_",family_test) | grepl("exp",family_test))
    model_result <- assoc_model_nls(family_test, tempDataFrame, sig.formula, transformation_y, plot, samples_sql_condition=samples_sql_condition,key)

  if(grepl("mediation-quantreg",family_test))
    model_result <- assoc_mediation_lqm(family_test, tempDataFrame, sig.formula, transformation_y, plot, samples_sql_condition=samples_sql_condition,key)

  if(grepl("mediation-ridge",family_test))
    model_result <- assoc_mediation_ridge_model(family_test, tempDataFrame, sig.formula, transformation_y, plot, samples_sql_condition=samples_sql_condition,key)

  if(grepl("mediation-linear",family_test))
    model_result <- assoc_mediation_linear_model(family_test, tempDataFrame, sig.formula, transformation_y, plot, samples_sql_condition=samples_sql_condition,key)

  if (grepl("quantreg-permutation", family_test))
  {
    # Determine the null device for the current platform
    null_device <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
    # Redirect output to the null device
    sink(null_device)
    model_result <- assoc_quantreg_permutation_model(family_test, sig.formula, tempDataFrame, independent_variable, transformation_y, plot , samples_sql_condition=samples_sql_condition,key)
    sink()
  }
  if (grepl("quantreg", family_test) & !grepl("quantreg-permutation", family_test)  & !grepl("mediation-quantreg", family_test))
    model_result <- assoc_quantreg_model(family_test, sig.formula, tempDataFrame, independent_variable, transformation_y, plot,
      samples_sql_condition=samples_sql_condition,key)

  #   return (model_result)
  # }, error = function(e) {
  #   return (model_result)
  # })

  #
  model_result <- as.data.frame(model_result)
  return (model_result)

}
