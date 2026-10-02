assoc_model_nls <- function (family_test, tempDataFrame, sig.formula, transformation_y, plot, samples_sql_condition=samples_sql_condition,key)
{
  area <- as.character(key$AREA)
  subarea <- as.character(key$SUBAREA)
  marker <- as.character(key$MARKER)
  figure <- as.character(key$FIGURE)

  ssEnv <- core_get_session_info()
  nls_mdl_params <- unlist(strsplit(as.character(family_test),"_"))
  partition_percentage <- as.numeric(nls_mdl_params[2])

  res<- data.frame("PARTITION_PERC" = as.numeric(partition_percentage))

  tempDataFrame <- as.data.frame(tempDataFrame)

  dep_var <- assoc_sig_formula_vars(sig.formula)
  dependent_variable <- dep_var$dependent_variable
  independent_variable <- dep_var$independent_variable
  covariates <- dep_var$covariates
  applied_function <- nls_mdl_params[1]

  if(length(nls_mdl_params)==3)
    if(nls_mdl_params[3]=="predictor")
    {
      dependent_variable <- dep_var$independent_variable
      independent_variable <- dep_var$dependent_variable
    }


  if (grepl("^log",applied_function))
  {
    # Bias correction to avoid applied_function  of 0 or negative numbers
    bias_independent_variable <- min(tempDataFrame[, independent_variable])
    bias_covariates <- min(tempDataFrame[, covariates])
    bias <- min(bias_covariates,bias_independent_variable, 0)
    if (bias <=0)
    {
      tempDataFrame[, independent_variable] <- tempDataFrame[, independent_variable] - ceiling(bias) + 1
      tempDataFrame[, covariates] <- tempDataFrame[, covariates] - ceiling(bias) + 1
    }
  }

  if (grepl("^exp",applied_function))
  {
    # Bias correction to avoid applied_function  of 0 or negative numbers
    bias_dependent_variable <- min(tempDataFrame[, dependent_variable])
    bias <- min(bias_dependent_variable, 0)
    if (bias <0)
      tempDataFrame[, dependent_variable] <- tempDataFrame[, dependent_variable] - ceiling(bias) + 1
  }

  # Split data into training and test sets
  training.samples <- caret::createDataPartition(tempDataFrame[, dependent_variable], p = partition_percentage, list = FALSE)
  train.data <- tempDataFrame[training.samples, ]
  test.data <- tempDataFrame[-training.samples, ]


  cov_temp <- covariates
  ### changes
  for (cov in cov_temp) {
    # apply ", applied_function , "  on independent variable and covariates
    cov_temp[cov_temp==cov] <- paste0("", applied_function , " (", cov, ")")
  }
  ### changes

  # Build dynamic formula for lm (to get starting values)
  if (length(covariates) > 0)
    lm_formula_str <- paste0(dependent_variable, " ~ ", applied_function , " (", independent_variable, ") + ", paste(cov_temp, collapse = " + "))
  else
    lm_formula_str <- paste0(dependent_variable, " ~ ", applied_function , " (", independent_variable, ")")
  lm_formula <- as.formula(lm_formula_str)
  lm_fit <- stats::lm(lm_formula, data = train.data)

  # Extract starting values for a and b
  start_list <- coef(lm_fit)

  names(start_list) <- letters[seq_len(length(c(independent_variable,covariates))+1)]
  cov_letters <- names(start_list)[-1][-1]

  # Build dynamic formula for nls
  if(length(covariates) == 0)
    rhs <- paste0("a + b * ", applied_function , " (", independent_variable, ")")
  else
    rhs <- paste0("a + b * ", applied_function , " (", independent_variable, ")",paste(paste0(" + ", cov_letters, " * (", covariates, ")"), collapse = ""))
    # rhs <- paste0("a + b * ", applied_function , " (", independent_variable, ")",paste(paste0(" + ", cov_letters, " * ", applied_function , " (", covariates, ")"), collapse = ""))
  nls_formula <- as.formula(paste0(dependent_variable, " ~ ", rhs))

  # Fit the NLS model
  nls_model_model_result <- tryCatch({
    stats::nls(
      formula = nls_formula,
      start = start_list,
      data = train.data,
      control = stats::nls.control(warnOnly = TRUE, maxiter = 1000, minFactor = 1e-8)
    )
  }, error = function(e) {
    # Handle error
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), "Error fitting the model: ", e$message)
    stop()
  })
  # check if nls_model_model_result is a model
  if (!inherits(nls_model_model_result, "nls"))
    return(res)

  predictions <- stats::predict(nls_model_model_result)
  if(partition_percentage < 1)
  {
    # Make predictions
    predictions_test <- stats::predict(nls_model_model_result, newdata = test.data)
    res <- cbind(res, assoc_model_performance(predictions, train.data[,dependent_variable], predictions_test, test.data[,dependent_variable]))
  }
  else
    res <- cbind(res, assoc_model_performance(predictions, train.data[,dependent_variable], c(),c()))

  # If model fitting was successful
  tryCatch({
    summary_result <- summary(nls_model_model_result)
    coefficients <- summary_result$coefficients
  }, error = function(e) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), "Error extracting coefficients from the model: ", e$message)
    return(res)
  })

  if(!exists("summary_result"))
    return(res)

  if (!is.null(summary_result) && ncol(summary_result$coefficients) >= 4) {
    # safe to extract p-values
  } else {
    warning("P-values could not be computed due to singular matrix or failed model summary.")
  }

  rownames(coefficients) <- core_name_cleaning(c("INTERCEPT", independent_variable, covariates))

  # for a and b extract the p-value
  for (i in seq_len(nrow(coefficients))) {
    # i <- 1
    p_value <- coefficients[i,4]
    row_name <- rownames(coefficients)[i]
    pval_name <- core_name_cleaning(paste0(row_name,"_PVALUE",sep=""))
    p_value <- data.frame(p_value)
    colnames(p_value) <- pval_name

    nls_mdl_coef_estimate <- data.frame(nls_mdl_a_estimate = coefficients[i,1])
    colnames(nls_mdl_coef_estimate) <- core_name_cleaning(paste0(row_name,"_ESTIMATE",sep=""))

    res <- cbind(res, p_value, nls_mdl_coef_estimate)
  }

  # remove rowname from res
  rownames(res) <- NULL
  # res$sig.fromula <- as.character(sig.formula)

  if(plot & !any(is.na(predictions)))
  {

    # predicted <- stats::predict(nls_model_model_result, train.data, type = "response")
    # train.data$predicted <- predicted
    # train.data$fitted <- nls_model_model_result$m$fitted

    chartFolder <- io_dir_check_and_create(ssEnv$result_folderChart,c("FITTED_MODEL", core_name_cleaning(samples_sql_condition)))
    filename  <-  io_file_path_build(chartFolder,
      c(as.character(family_test), independent_variable,"Vs",as.character(transformation_y), dependent_variable, covariates, key$COMBINED),
      ssEnv$plot_format)

    ggp <- ggplot2::ggplot(train.data, ggplot2::aes(x =!!ggplot2::sym(independent_variable), y =!!ggplot2::sym(dependent_variable))) +
      ggplot2::geom_point(color = ssEnv$color_palette[1]) +
      ggplot2::geom_smooth(method = "lm", formula = as.formula(paste0("y ~ ", applied_function , " (x)")), se = TRUE, color = ssEnv$color_palette_darker[2]) + # The fitted line will be an exponential curve
      ggplot2::labs(
        title = "", # Puoi lasciare vuoto o aggiungere un titolo
        x = independent_variable,
        y = dependent_variable
      )

    if(partition_percentage < 1)
    {
      ggp <- ggp + ggplot2::geom_point(data = test.data, ggplot2::aes(y = predictions_test, x = eval(parse(text=independent_variable))), color = ssEnv$color_palette_darker[3])+
        ggplot2::geom_point(data = test.data, ggplot2::aes(y = dependent_variable, x = eval(parse(text=independent_variable))), color = ssEnv$ssEnv$color_palette[2])+
        ggplot2::xlab(independent_variable) +
        ggplot2::ylab(dependent_variable)
    }

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
  return (res)
}

