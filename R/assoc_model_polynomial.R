assoc_model_polynomial <- function (family_test, tempDataFrame, sig.formula , transformation_y, plot, samples_sql_condition=samples_sql_condition, key)
{

  area <- as.character(key$AREA)
  subarea <- as.character(key$SUBAREA)
  marker <- as.character(key$MARKER)
  figure <- as.character(key$FIGURE)

  ssEnv <- core_get_session_info()

  # polynomial_degree_partition-partition_percentage
  polynomial_params <- unlist(strsplit(as.character(family_test),"_"))

  if(length(polynomial_params)!=3)
  {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " The polynomial model must have 3 parameters: polynomial followed by degree and partition percentage eg: polynomial_2_1 " )
  }

  degree <- as.numeric(polynomial_params[2])
  res <- data.frame("PL_DEGREE"= degree)
  partition_percentage <- as.numeric(polynomial_params[3])
  res$PL_PERC <- partition_percentage
  # res$sig.fromula <- as.character(sig.formula)

  if(length(polynomial_params)==4)
    res$r_model <- paste0("polynomial",polynomial_params[4], sep="")
  else
    res$r_model <- "polynomial"


  tempDataFrame <- as.data.frame(tempDataFrame)
  dep_var <- assoc_sig_formula_vars(sig.formula)
  dependent_variable <- dep_var$dependent_variable
  independent_variable <- dep_var$independent_variable
  covariates <- dep_var$covariates

  if(length(polynomial_params)==4)
    if(polynomial_params[4]=="predictor")
    {
      dependent_variable <- dep_var$independent_variable
      independent_variable <- dep_var$dependent_variable
    }

  if (nrow(tempDataFrame) == 0)
    return(res)

  # Split the data into training and test set
  training.samples <- tempDataFrame[, dependent_variable] %>% caret::createDataPartition(p = partition_percentage, list = FALSE)

  train.data  <- tempDataFrame[training.samples, ]
  test.data <- tempDataFrame[-training.samples, ]

  # One formula for both cases. The no-covariate branch used to write its own
  # model, stats::poly(eval(parse(text = ...)), degree, raw = TRUE), which spans
  # the same column space as I(x^1) + ... + I(x^degree) - a raw poly IS the
  # monomials - so the fit was identical and only the coefficient names differed.
  # They differed badly; see the loop below.
  formula <- assoc_polynomial_formula_build(dependent_variable, independent_variable,
                                            degree, covariates)
  polynomial_model_result <- stats::lm(formula, data = train.data,
                                       na.action = stats::na.exclude)



  # check id polynomial_model_result is null
  if(is.null(polynomial_model_result))
    return(res)

  predictions <- stats::predict(polynomial_model_result)
  if(partition_percentage < 1)
  {
    # Make predictions
    predictions_test <- stats::predict(polynomial_model_result, newdata = test.data)
    res <- cbind(res, assoc_model_performance(predictions, train.data[,dependent_variable], predictions_test, test.data[,dependent_variable]))
  }
  else
    res <- cbind(res, assoc_model_performance(predictions, train.data[,dependent_variable], c(),c()))

  # Coefficients and Confidence Intervals
  coefficients <- coef(summary(polynomial_model_result))
  # conf_int <- confint(polynomial_model_result)

  # One cleaning for the three columns of a term, so a p-value, its estimate and
  # its standard error carry the same name. They did not: the p-value ran through
  # three gsub calls and the estimate through none, so the estimate came out as
  # STATS_POLY_EVAL_PARSE_TEXT_EQ_INDEPENDENT_VARIABLE_DEGREE_RAW_EQ_TRUE_1,
  # carrying the literal string INDEPENDENT_VARIABLE instead of the name of the
  # variable - two runs on two different predictors produced the same column
  # name - while its own p-value was in a column named for the variable. The
  # first of those gsub calls never fired at all: its pattern begins with an
  # underscore and the text it was meant to strip begins the name. All three are
  # gone with the model that needed them, because I(x^k) cleans to I_<VAR>_<k>
  # by itself.
  #
  # STD_ERROR is new. The value was always in coefficients[, 2] and was never
  # written, which is what left a cross-study pooling of these coefficients with
  # nothing to weight them by.
  coefficient_columns <- c(ESTIMATE = 1L, STD_ERROR = 2L, PVALUE = 4L)

  for (i in seq_len(nrow(coefficients))) {
    term_name <- core_name_cleaning(rownames(coefficients)[i])

    for (what in names(coefficient_columns)) {
      value <- data.frame(coefficients[i, coefficient_columns[[what]]])
      colnames(value) <- paste0(term_name, "_", what)
      res <- cbind(res, value)
    }
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
      c(as.character(family_test), independent_variable,"Vs",as.character(transformation_y), dependent_variable, covariates, key$COMBINED),
      ssEnv$plot_format)

    # Predict the values for the plot
    train.data$predicted <- predict(polynomial_model_result, newdata = train.data)

    if(length(covariates)>0)
    {
      # Plot the data and the polynomial fit
      ggp <- ggplot2::ggplot(train.data, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
        ggplot2::geom_point(color = ssEnv$color_palette[1]) +
        ggplot2::stat_smooth(method = lm, formula = y ~ poly(x, degree, raw = TRUE), color = ssEnv$color_palette_darker[3]) +
        ggplot2::xlab(independent_variable) +
        ggplot2::ylab(dependent_variable) +
        ggplot2::ggtitle("")
    }
    else
    {
      # do a plot with train.data, test.data and predictions with 3 different colors 1 color for train.data, 1 color for test.data and 1 color for predictions
      ggp <- ggplot2::ggplot(train.data, ggplot2::aes(eval(parse(text=independent_variable)), eval(parse(text=dependent_variable))) ) +
        ggplot2::geom_point( color = ssEnv$color_palette[1] ) +
        ggplot2::stat_smooth(method = lm, formula = y ~ poly(x, degree, raw = TRUE), color = ssEnv$color_palette_darker[3]) +
        ggplot2::xlab(independent_variable) +
        ggplot2::ylab(dependent_variable)
    }

    if(partition_percentage < 1)
      ggp <- ggp + ggplot2::geom_line(ggplot2::aes_string(y = "predicted"), color = ssEnv$color_palette_darker[2])


    if (partition_percentage < 1)
      ggp <- ggp + ggplot2::geom_point(data = test.data, ggplot2::aes(y = predictions_test, x = eval(parse(text=independent_variable))), color = ssEnv$color_palette_darker[3]) +
      ggplot2::xlab(independent_variable) +
      ggplot2::ylab(dependent_variable)


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
  # ggplot2::ggplot(train.data, ggplot2::aes(eval(parse(text=independent_variable)), eval(parse(text=dependent_variable))) ) +
  #   ggplot2::geom_point( color = ssEnv$color_palette[1] ) +
  #   ggplot2::stat_smooth(method = lm, formula = y ~ stats::poly(x, degree, raw = TRUE))
  #
  # # do a plot with train.data, test.data and predictions with 3 different colors 1 color for train.data, 1 color for test.data and 1 color for predictions
  # ggplot2::ggplot(train.data, ggplot2::aes(eval(parse(text=independent_variable)), eval(parse(text=dependent_variable))) ) +
  #   ggplot2::geom_point( color = ssEnv$color_palette[1] ) + ggplot2::stat_smooth(method = lm, formula = y ~ poly(x, degree, raw = TRUE)) +
  #   ggplot2::geom_point(data = test.data, ggplot2::aes(y = predictions), color = ssEnv$color_palette[2])
  #
  # # do a plot with train.data, test.data and predictions with 3 different colors 1 color for train.data, 1 color for test.data and 1 color for predictions
  # ggplot2::ggplot(train.data, ggplot2::aes(eval(parse(text=independent_variable)), eval(parse(text=dependent_variable))) ) +
  #   ggplot2::geom_point( color = ssEnv$color_palette[1] ) + ggplot2::stat_smooth(method = lm, formula = y ~ poly(x, degree, raw = TRUE)) +
  #   ggplot2::geom_point(data = test.data, ggplot2::aes(y = predictions), color = ssEnv$color_palette[2]) +
  #   ggplot2::geom_point(data = data.frame(train.data,polynomial_model_result$residuals) , ggplot2::aes(y = polynomial_model_result$residuals), color = "cyan")
  #
  return (res)
}

# polynomial_model_result("polynomial_4_0.8", sample_sheet, "DELTARQ_HYPO ~ Age")
