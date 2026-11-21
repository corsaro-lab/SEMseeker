# Function to calculate mediation effects using Ridge Regression
assoc_mediation_ridge_model <-  function(family_test,tempDataFrame, sig.formula, transformation_y, plot, samples_sql_condition=samples_sql_condition, key)
{
  # function(tempDataFrame, treatment, mediator, outcome,covariates, lambda = 0.1, permutations = 1000) {


  area <- as.character(key$AREA)
  subarea <- as.character(key$SUBAREA)
  marker <- as.character(key$MARKER)
  figure <- as.character(key$FIGURE)

  ssEnv <- core_get_session_info()

  ridge_params <- unlist(strsplit(as.character(family_test),"_"))
  lambda <- as.numeric(ridge_params[2])
  res <- data.frame("lambda" = lambda)
  permutations_test <- as.numeric(ridge_params[3])
  permutations <- as.numeric(ridge_params[4])

  # decompose the formula, to move the mediator
  # the expected formula is "mediator ~ outcome + treatment + covariates"
  vars <- assoc_sig_formula_vars(sig.formula)
  mediator <- vars[[1]] # burden_Value
  outcome <- vars[[2]] # dependent variable
  treatment <- vars[[3]][1]

  # core_log_event("DEBUG: ",  format(Sys.time(), "%a %b %d %X %Y"), " - Starting mediation analysis for: ", mediator, " and ", outcome)

  res$mediator <- mediator
  res$outcome <- outcome
  res$treatment <- treatment

  if(length(vars[[3]])>1)
  {
    covariates <- vars[[3]][-c(1)]
    # transform covariates in a vector of string
    covariates <- vapply(covariates, as.character, character(1))
  } else
    covariates <- c()


  # Prepare data matrices for glmnet::glmnet
  if(length(covariates) != 0)
    x <- as.matrix(tempDataFrame[, c(treatment, covariates)])
  else
    x <- as.matrix(tempDataFrame[, c(treatment)])
  z <- as.matrix(tempDataFrame[, mediator])
  y <- as.matrix(tempDataFrame[, outcome])

  # Ensure the matrices are correctly dimensioned
  x <- matrix(x, ncol = (1 + length(covariates)))
  z <- matrix(z, ncol = 1)
  y <- matrix(y, ncol = 1)

  perms <- sort(unique(c(permutations_test, permutations)))

  tryCatch(
    {
      # Cross-validation to select the best lambda
      cv_model_mediator <- glmnet::cv.glmnet(x, z, alpha=0)
      # Visualizza il miglior valore di lambda
      best_lambda_mediator <- cv_model_mediator$lambda.min
      res$best_lambda_mediator <- best_lambda_mediator

      # Cross-validation to select the best lambda
      cv_model_outcome <- glmnet::cv.glmnet(cbind(x, z), y, alpha=0)
      # Visualizza il miglior valore di lambda
      best_lambda_outcome<- cv_model_outcome$lambda.min
      res$best_lambda_outcome <- best_lambda_outcome

      # Fit the mediator model using Ridge Regression
      mediator_model <- glmnet::glmnet(x, z, alpha = 0, lambda = best_lambda_mediator)
      treatment_on_mediator_effect <- coef(mediator_model)[2]

      predictions <- predict(mediator_model, newx = x, s = "lambda.min")
      expected_values <- z
      fitted_values <- predictions
      res_mdl <- assoc_model_performance(fitted_values, expected_values, c(), c())
      res <- cbind(res, res_mdl)

      # Fit the outcome model using Ridge Regression
      outcome_model <- glmnet::glmnet(cbind(x, z), y, alpha = 0, lambda = best_lambda_outcome)
      mediator_on_outcome_effect <- coef(outcome_model)[3]
      direct_effect <- coef(outcome_model)[2]

      # Calculate the indirect effect (indirect_effect)
      indirect_effect <- treatment_on_mediator_effect * mediator_on_outcome_effect

      # Calculate the total effect (c = c' + indirect_effect)
      total_effect <- direct_effect + indirect_effect
    },
    error = function(e) {
      core_log_event("DEBUG: ",  format(Sys.time(), "%a %b %d %X %Y"), " - Error in mediation analysis for: ", mediator, " and ", outcome)
      core_log_event("DEBUG: ",  format(Sys.time(), "%a %b %d %X %Y"), " - Error message: ", e$message)
    },
    finally = {
      if (!exists("mediator_on_outcome_effect") | !exists("treatment_on_mediator_effect") |
          !exists("direct_effect") | !exists("indirect_effect") | !exists("total_effect"))
        return(res)
      # Save the results
      res$mediator_on_outcome_effect <- mediator_on_outcome_effect
      res$treatment_on_mediator_effect <- treatment_on_mediator_effect
      res$direct_effect <- direct_effect
      res$indirect_effect <- indirect_effect
      res$total_effect <- total_effect
    }
  )

  # Function to calculate mediation effects for bootstrapped samples
  boot_mediation_ridge <- function(tempDataFrame, indices) {
    tryCatch(
      {
        boot_data <- tempDataFrame[indices, ]
        if(length(covariates) != 0)
          x_boot <- as.matrix(boot_data[, c(treatment, covariates)])
        else
          x_boot <- as.matrix(boot_data[, c(treatment)])
        z_boot <- as.matrix(boot_data[, mediator])
        y_boot <- as.matrix(boot_data[, outcome])

        x_boot <- matrix(x_boot, ncol = (1 + length(covariates)))
        z_boot <- matrix(z_boot, ncol = 1)
        y_boot <- matrix(y_boot, ncol = 1)

        mediator_model <- glmnet::glmnet(x_boot, z_boot, alpha = 0, lambda = best_lambda_mediator)
        outcome_model <- glmnet::glmnet(cbind(x_boot, z_boot), y_boot, alpha = 0, lambda = best_lambda_outcome)
        treatment_on_mediator_effect <- coef(mediator_model)[2]
        mediator_on_outcome_effect <- coef(outcome_model)[3]
        direct_effect <- coef(outcome_model)[2]
        indirect_effect <- treatment_on_mediator_effect * mediator_on_outcome_effect
        total_effect <- direct_effect + indirect_effect
      },
      error = function(e) {
        core_log_event("DEBUG: ",  format(Sys.time(), "%a %b %d %X %Y"), " - Error in bootstrapping for: ", mediator, " and ", outcome)
        core_log_event("DEBUG: ",  format(Sys.time(), "%a %b %d %X %Y"), " - Error message: ", e$message)
        return(c(NA, NA, NA, NA, NA))
      }, finally = {
        if (!exists("mediator_on_outcome_effect") | !exists("treatment_on_mediator_effect") |
            !exists("direct_effect") | !exists("indirect_effect") | !exists("total_effect"))
          return(NULL)
        return(c(treatment_on_mediator_effect, mediator_on_outcome_effect, direct_effect, indirect_effect, total_effect))
      })
  }

  for (i in seq_along(perms)) {

    permutations <- perms[i]

    # Perform bootstrapping
    boot_results <- boot::boot(data = tempDataFrame, statistic = boot_mediation_ridge, R = permutations)

    if(is.null(boot_results))
      return(res)

    # Calculate p-values
    p_value <- function(observed, bootstrapped) {
      (sum(abs(bootstrapped) >= abs(observed)) + 1) / (length(bootstrapped) + 1)
    }

    p_values <- vapply(seq_len(4), function(i) p_value(c(treatment_on_mediator_effect, mediator_on_outcome_effect, direct_effect, indirect_effect)[i], boot_results$t[, i]), numeric(1))
    p_values <- as.data.frame(matrix(p_values, nrow = 1))
    colnames(p_values) <- c("Treatment_on_Mediator_Pvalue", "Mediator_on_Outcome_Pvalue", "Direct_Effect_Pvalue", "Indirect_Effect_Pvalue")

    res$pvalue <- max(p_values)
    res <- cbind(res, p_values)
    if(!all(p_values < as.numeric(ssEnv$alpha)))
      break

  }
  return(res)
}

# # Example usage
# set.seed(123)
# df <- data.frame(
#   treatment = rbinom(100, 1, 0.5),
#   mediator = rnorm(100),
#   outcome = rnorm(100),
#   covariates = rnorm(100)
# )
#
# # Calculate mediation effects with Ridge Regression and bootstrapping
# mediation_result <- calculate_mediation_effect_ridge(df, treatment = "treatment", mediator = "mediator", outcome = "outcome", covariates = "covariates")
#
# # Print the results
# print(mediation_result)

