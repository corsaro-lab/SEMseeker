#' GLM association model
#'
#' @param family_test regression model family to apply (e.g. "gaussian", "binomial", "poisson")
#' @param tempDataFrame data frame to use for the model
#' @param sig.formula formula to apply the model
#' @param transformation_y transformation applied to the dependent variable (for labelling plots)
#' @param plot logical; if TRUE, generate and save a scatter/fit plot
#' @param samples_sql_condition SQL condition string used to filter samples (used for plot file naming)
#' @param key named list with AREA, SUBAREA, MARKER and FIGURE identifiers for this test
#'
#' @return A list (as a data.frame row) with GLM fit results including AIC,
#'   per-coefficient p-values and estimates, model performance metrics, and
#'   optionally residual diagnostics.
#'
assoc_glm_model <- function(family_test, tempDataFrame, sig.formula, transformation_y, plot, samples_sql_condition=samples_sql_condition, key)
{

  # The plotting branches below read ssEnv$color_palette. The lookup was
  # missing, so every plot = TRUE call died on an unbound ssEnv.
  ssEnv <- core_get_session_info()

  area <- as.character(key$AREA)
  subarea <- as.character(key$SUBAREA)
  marker <- as.character(key$MARKER)
  figure <- as.character(key$FIGURE)

  result_glm  <- stats::glm( sig.formula, family = as.character(family_test), data = as.data.frame(tempDataFrame))

  # AI-044 (2026-06-09): emit AIC_VALUE + BIC_VALUE in metrics_properties.rda
  # canonical uppercase form. Legacy lowercase aic_value retained for
  # backward-compat with downstream consumers that read it by exact name.
  res <- data.frame(
    "aic_value" = result_glm$aic,
    "AIC_VALUE" = as.numeric(stats::AIC(result_glm)),
    "BIC_VALUE" = as.numeric(stats::BIC(result_glm))
  )

  dep_var <- assoc_sig_formula_vars(sig.formula)
  dependent_variable <- dep_var$dependent_variable
  independent_variable <- dep_var$independent_variable

  coefficients <- summary(result_glm)$coeff
  # extract the p-value
  for (i in seq_len(nrow(coefficients))) {
    # i <- 1
    p_value <- coefficients[i,4]
    row_name <- rownames(coefficients)[i]
    pval_name <- core_name_cleaning(paste0(row_name,"_PVALUE",sep=""))
    p_value <- data.frame(p_value)
    colnames(p_value) <- pval_name

    coef_estimate <- data.frame(coef_estimate = coefficients[i,1])
    colnames(coef_estimate) <- core_name_cleaning(paste0(row_name,"_ESTIMATE",sep=""))

    res <- cbind(res, p_value, coef_estimate)
  }

  # faking predictions until next improvement
  training.samples <- caret::createDataPartition(tempDataFrame[, dependent_variable], p = 1, list = FALSE)
  train.data <- tempDataFrame[training.samples, ]
  test.data <- tempDataFrame[-training.samples, ]
  predictions <- stats::predict(result_glm)
  res <- cbind(res, assoc_model_performance(predictions, train.data[,dependent_variable], c(),c()))

  if(family_test=="gaussian")
  {
    ggp <- NULL
    residuals <-  result_glm$resid
    #calculate shapiro of working residuals
    res$residuals_normality <- as.logical(if(length(residuals)>3 & length(unique(residuals))>3) (stats::shapiro.test(residuals)$p.value > 0.05) else NA)
    # calculate residuals homoschedadasticity
    res$residuals_homoschedasticity <- as.logical(lmtest::bptest(result_glm)$p.value > 0.05)

    if(plot)
    {
      #
      ggp <- ggplot2::ggplot(tempDataFrame, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
        ggplot2::geom_point( color = ssEnv$color_palette[1] ) +
        ggplot2::geom_smooth(method = "lm", se = TRUE, level = 0.95, color = ssEnv$color_palette_darker[3]) +
        ggplot2::ggtitle("") +
        ggplot2::xlab(independent_variable) +
        ggplot2::ylab(dependent_variable)

      # data_to_save <- cbind(train.data, predicted = apply(train.data[,independent_variable],1, function(x) (start_a * exp(start_b * x))))
      # colnames(data_to_save) <- c("Independent_Variable","Dependent_Variable","Predicted")
    }

  }

  if(family_test=="binomial")
  {
    ggp <- NULL
    if(plot)
    {
      # ggp <- plot_box_plot(dataFrameToPlot, independent_variable,burdenValue, transformation_y, family_test)
      #
      ggp <- ggplot2::ggplot(train.data, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
        ggplot2::geom_point( color = ssEnv$color_palette[1] ) +
        ggplot2::geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE, level = 0.95, color = ssEnv$color_palette_darker[3]) +
        ggplot2::ggtitle("") +
        ggplot2::xlab(independent_variable) +
        ggplot2::ylab(dependent_variable)

      # data_to_save <- cbind(train.data, predicted = apply(train.data[,independent_variable],1, function(x) (start_a * exp(start_b * x))))
      # colnames(data_to_save) <- c("Independent_Variable","Dependent_Variable","Predicted")
    }

  }

  if (!is.null(ggp) & plot==TRUE){

    chartFolder <- io_dir_check_and_create(ssEnv$result_folderChart,c("FITTED_MODEL", core_name_cleaning(samples_sql_condition)))
    filename  <-  io_file_path_build(chartFolder,
      c(as.character(family_test), independent_variable,"Vs",as.character(transformation_y), dependent_variable, covariates, key$COMBINED),
      ssEnv$plot_format)

    ggplot2::ggsave(
      filename,
      plot = ggp,
      scale = 1,
      width = 9,
      height = 9,
      units = c("in"),
      dpi = as.numeric(ssEnv$plot_resolution_ppi)
    )

  }

  #  calculate r-squared
  # res$r_squared <- summary(result_glm)$r.squared

  res$r_model <- "stats::glm"
  return (res)
}
