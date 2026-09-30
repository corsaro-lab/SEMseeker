assoc_multinomial_model <- function (family_test, tempDataFrame, sig.formula , transformation_y, plot, samples_sql_condition=samples_sql_condition, key)
{

  area <- as.character(key$AREA)
  subarea <- as.character(key$SUBAREA)
  marker <- as.character(key$MARKER)
  figure <- as.character(key$FIGURE)

  ssEnv <- core_get_session_info()

  # The family string is split so that a trailing "_predictor" suffix can swap
  # response and predictor, the same shape assoc_relu_model() uses. Without the
  # split, multinomial_params held the whole string, length() was always 1, and
  # the swap below could never fire. It still cannot fire today for a different
  # reason: assoc_validate_family() accepts the bare "multinomial" and no
  # suffixed form, so the branch is there for when it does.
  multinomial_params <- unlist(strsplit(as.character(family_test), "_"))
  res <-  data.frame("r_model","multinomial")

  tempDataFrame <- as.data.frame(tempDataFrame)
  dep_var <- assoc_sig_formula_vars(sig.formula)
  dependent_variable <- dep_var$dependent_variable
  independent_variable <- dep_var$independent_variable
  covariates <- dep_var$covariates

  if(length(multinomial_params)==4)
    if(multinomial_params[4]=="predictor")
    {
      dependent_variable <- dep_var$independent_variable
      independent_variable <- dep_var$dependent_variable
    }

  if (nrow(tempDataFrame) == 0)
    return(res)

  tempDataFrame[, dependent_variable] <- as.factor(tempDataFrame[, dependent_variable])
  multinomial_model_result <- nnet::multinom( sig.formula, data = tempDataFrame)

  # check id multinomial_model_result is null
  if(is.null(multinomial_model_result))
    return(res)


  # Wald test on the multinom fit: summary() carries both the coefficients and
  # their standard errors. The lines that READ model_summary were here, the line
  # that builds it was not, so this function could never return a result. The
  # coef() call above it was overwritten on the next line and is gone.
  model_summary <- summary(multinomial_model_result)
  coefficients <- model_summary$coefficients
  std_errors <- model_summary$standard.errors

  # Compute z-scores and p-values
  z_scores <- coefficients / std_errors
  p_values <- 2 * (1 - pnorm(abs(z_scores)))

  # Reshape to long format
  results_df <- data.frame(
    classe = rep(rownames(coefficients), times = ncol(coefficients)),
    variabile = rep(colnames(coefficients), each = nrow(coefficients)),
    coefficiente = as.vector(coefficients),
    errore_std = as.vector(std_errors),
    z_score = as.vector(z_scores),
    p_value = as.vector(p_values)
  )



  res$pvalue <- 0
  # for each degree extract the p-value
  for (i in seq_len(nrow(coefficients))) {
    # i <- 1
    p_value <- coefficients[i,4]
    row_name <- rownames(coefficients)[i]
    pval_name <- core_name_cleaning(paste0("pvalue_",row_name))
    pval_name <- core_name_cleaning(gsub("_STATS_POLY_EVAL_PARSE_TEXT_EQ","",pval_name))
    pval_name <- core_name_cleaning(gsub("_RAW_EQ_TRUE","",pval_name))
    pval_name <- core_name_cleaning(gsub("INDEPENDENT_VARIABLE",independent_variable,pval_name))
    p_value <- data.frame(p_value)
    colnames(p_value) <- pval_name
    res <- cbind(res, p_value)
    if (grepl(independent_variable, pval_name))
      res$pvalue <- max(c(res$pvalue, p_value[1,1]), na.rm = TRUE)
  }

  # remove rowname from res
  rownames(res) <- NULL
  if(plot)
  {
    chartFolder <- io_dir_check_and_create(ssEnv$result_folderChart,c("FITTED_MODEL", core_name_cleaning(samples_sql_condition)))

    if(is.null(covariates) || length(covariates)  ==  0)
      file_suffix <- ""
    else
    {
      long_covariates <- length(covariates) > 2
      # split each covariates by _
      if (long_covariates)
      {
        covariates <- unlist(t(strsplit( gsub(" ","",covariates),split  =  "_", fixed  =  TRUE)))
        covariates <- unique(covariates)
      }
      covariates <- paste(covariates, collapse = "_")
    }

    filename  <-  io_file_path_build(chartFolder,
      c(as.character(family_test), independent_variable,"Vs",as.character(transformation_y), dependent_variable, covariates, area, subarea),
      ssEnv$plot_format)

    # Predict the values for the plot
    # The model is fitted on the whole tempDataFrame: this function never split
    # train from test. The block below was transplanted from the relu model,
    # which does split, and arrived without partition_percentage, train.data,
    # test.data, predictions_test or degree. Plot what is actually here.
    plot_data <- tempDataFrame
    plot_data$predicted <- stats::predict(multinomial_model_result, newdata = plot_data)

    if(length(covariates)>0)
    {
      # Plot the data and the multinomial fit
      ggp <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
        ggplot2::geom_point(color = ssEnv$color_palette[1]) +
        ggplot2::geom_line(ggplot2::aes_string(y = "predicted"), color = ssEnv$color_palette_darker[2]) +
        ggplot2::xlab(independent_variable) +
        ggplot2::ylab(dependent_variable) +
        ggplot2::ggtitle("")
      # formula_string <- paste0(dependent_variable, " ~ poly(", independent_variable, ", ", degree, ", raw = TRUE)",
      #   ifelse(length(covariates) > 0, paste0(" + ", paste(covariates, collapse = " + ")), ""))
      #
      # ggp <- ggplot(train.data, aes(x = independent_variable, y = dependent_variable)) +
      #   ggplot2::geom_point( color = ssEnv$color_palette[1] ) +
      #   stat_smooth(method = lm,
      #     formula = formula,
      #     color = ssEnv$color_palette_darker[2]) +
      #   ggplot2::theme_minimal()
    }
    else
    {
      # do a plot with train.data, test.data and predictions with 3 different colors 1 color for train.data, 1 color for test.data and 1 color for predictions
      ggp <- ggplot2::ggplot(plot_data, ggplot2::aes(eval(parse(text=independent_variable)), eval(parse(text=dependent_variable))) ) +
        ggplot2::geom_point( color = ssEnv$color_palette[1] ) +
        ggplot2::geom_line(ggplot2::aes_string(y = "predicted"), color = ssEnv$color_palette_darker[2]) +
        ggplot2::xlab(independent_variable) +
        ggplot2::ylab(dependent_variable)
    }



    ggplot2::ggsave(
      filename,
      plot = ggp,
      scale = 1,
      width = 8,
      height = 8,
      units = c("in"),
      dpi = as.numeric(ssEnv$plot_resolution_ppi)
    )

    # data_to_save <- cbind(train.data, predicted = apply(train.data[,independent_variable],1, function(x) (poly(x, degree, raw = TRUE))))
    # colnames(data_to_save) <- c("Independent_Variable","Dependent_Variable","Predicted")

  }
  # # do a plot with train.data, test.data and predictions with 3 different colors 1 color for train.data, 1 color for test.data and 1 color for predictions
  # ggplot2::ggplot(plot_data, ggplot2::aes(eval(parse(text=independent_variable)), eval(parse(text=dependent_variable))) ) +
  #   ggplot2::geom_point( color = ssEnv$color_palette[1] ) +
  #   ggplot2::stat_smooth(method = lm, formula = y ~ stats::poly(x, degree, raw = TRUE))
  #
  # # do a plot with train.data, test.data and predictions with 3 different colors 1 color for train.data, 1 color for test.data and 1 color for predictions
  # ggplot2::ggplot(plot_data, ggplot2::aes(eval(parse(text=independent_variable)), eval(parse(text=dependent_variable))) ) +
  #   ggplot2::geom_point( color = ssEnv$color_palette[1] ) + ggplot2::stat_smooth(method = lm, formula = y ~ poly(x, degree, raw = TRUE)) +
  #   ggplot2::geom_point(data = test.data, ggplot2::aes(y = predictions), color = ssEnv$color_palette[2])
  #
  # # do a plot with train.data, test.data and predictions with 3 different colors 1 color for train.data, 1 color for test.data and 1 color for predictions
  # ggplot2::ggplot(plot_data, ggplot2::aes(eval(parse(text=independent_variable)), eval(parse(text=dependent_variable))) ) +
  #   ggplot2::geom_point( color = ssEnv$color_palette[1] ) + ggplot2::stat_smooth(method = lm, formula = y ~ poly(x, degree, raw = TRUE)) +
  #   ggplot2::geom_point(data = test.data, ggplot2::aes(y = predictions), color = ssEnv$color_palette[2]) +
  #   ggplot2::geom_point(data = data.frame(train.data,multinomial_model_result$residuals) , ggplot2::aes(y = multinomial_model_result$residuals), color = "cyan")
  #
  return (res)
}

# multinomial_model_result("multinomial_4_0.8", sample_sheet, "DELTARQ_HYPO ~ Age")
