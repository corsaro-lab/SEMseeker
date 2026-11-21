assoc_mediation_lqm <- function(family_test,tempDataFrame, sig.formula, transformation_y, plot, samples_sql_condition=samples_sql_condition, key)
{
  # mediation_lqm_model
  ssEnv <- core_get_session_info()
  area <- as.character(key$AREA)
  subarea <- as.character(key$SUBAREA)
  marker <- as.character(key$MARKER)
  figure <- as.character(key$FIGURE)

  quantreg_params <- unlist(strsplit(as.character(family_test),"_"))
  tau <- as.numeric(quantreg_params[2])
  res <- data.frame("tau" = tau)
  permutations_test <- as.numeric(quantreg_params[3])
  permutations <- as.numeric(quantreg_params[4])

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

  df_mediator <- data.frame(tempDataFrame)
  df_mediator$mediator <- df_mediator[,mediator]
  df_mediator$treatment <- df_mediator[,treatment]
  if(length(covariates) != 0)
    df_mediator <- df_mediator[,c("treatment", "mediator",covariates)]
  else
    df_mediator <- df_mediator[,c("treatment", "mediator")]

  df_outcome <- data.frame(tempDataFrame)
  df_outcome$outcome <- df_outcome[,outcome]
  df_outcome$mediator <- df_outcome[,mediator]
  df_outcome$treatment <- df_outcome[,treatment]
  if(length(covariates) != 0)
    df_outcome <- df_outcome[,c("outcome", "mediator", "treatment",covariates)]
  else
    df_outcome <- df_outcome[,c("outcome", "mediator", "treatment")]

  df_direct <- data.frame(tempDataFrame)
  df_direct$outcome <- df_direct[,outcome]
  df_direct$treatment <- df_direct[,treatment]
  if(length(covariates) != 0)
    df_direct <- df_direct[,c("outcome", "treatment",covariates)]
  else
    df_direct <- df_direct[,c("outcome", "treatment")]

  perms <- sort(unique(c(permutations_test, permutations)))


  lqm_control <- list(loop_tol_ll = 1e-5, loop_max_iter = 10000, verbose = FALSE )


  if(assoc_calculate_collinearity_score(df_outcome[, -which(colnames(df_outcome) == "outcome")]) > 0.5)
  {
    core_log_event("DEBUG: ",  format(Sys.time(), "%a %b %d %X %Y"), " - Collinearity detected in the mediator model")
    res$collinearity_mediator <- TRUE
    return(res)
  } else
    res$collinearity_mediator <- FALSE

  for ( p in seq_len(2))
  {
    permutations <- perms[p]
    res$permutations <- permutations
    tryCatch(
      expr ={
        res <- data.frame("permutations"=permutations)

        b1_fun <- function(data,i){
          # model.x <- lm("outcome ~ . ", data= data[i,] , na.action =stats::na.omit)
          # b1 <- model.x$coefficients[2]
          model.x <- lqmm::lqm("outcome ~ . ", tau=tau,  data= data[i,] , na.action =stats::na.omit, control = lqm_control)
          b1 <- model.x$theta[2]
          return(b1)
        }

        b2_fun <- function(data,i){
          # model.m <- lm("mediator ~ . ", data= data[i,] , na.action =stats::na.omit)
          # b2 <- model.m$coefficients[2]
          model.m <- lqmm::lqm("mediator ~ . ", tau=tau, data= data[i,], na.action = stats::na.omit, control = lqm_control)
          b2 <- model.m$theta[2]
          return(b2)
        }

        b3_fun <- function(data,i){
          # model.y <- lm("outcome ~ . ", data= data[i,] , na.action =stats::na.omit)
          # b3 <- model.y$coefficients[2]
          model.y <- lqmm::lqm(formula = as.formula("outcome ~ . "), tau=tau,  data= data[i,] , na.action = stats::na.omit, control = lqm_control)
          b3 <- model.y$theta[3]
          return(b3)
        }

        b1_res <- boot::boot(data=df_direct, statistic = b1_fun, R=permutations, parallel = "no", ncpus=1)
        b3_res <- boot::boot(data=df_outcome, statistic = b3_fun, R=permutations, parallel = "no", ncpus=1)

        mediation_fun <- function(data, i){
          model.y <- lqmm::lqm(formula = as.formula("outcome ~ ."), tau=tau,  data= data[i,] , na.action = stats::na.omit, control = lqm_control)
          b3 <- model.y$theta[3]
          # model.x <- lm("outcome ~ treatment + mediator ", data= data[i,] , na.action =stats::na.omit)
          # b3 <- model.x$coefficients[2]
          # remove mediator from the data
          data[i,] <- data[i,][, -which(colnames(data[i,]) == "mediator")]
          model.x <- lqmm::lqm("outcome ~ . ", tau=tau,  data= data[i,] , na.action =stats::na.omit, control = lqm_control)
          b1 <- model.x$theta[2]
          # model.m <- lm("outcome ~ treatment ", data= data[i,] , na.action =stats::na.omit)
          # b1 <- model.m$coefficients[2]
          mediation_res <- b1 - b3
          return(mediation_res)
        }

        prop_effect_fun <- function(data, i){
          model.y <- lqmm::lqm(formula = as.formula("outcome ~ ."), tau=tau,  data= data[i,] , na.action = stats::na.omit, control = lqm_control)
          b3 <- model.y$theta[3]
          # model.x <- lm("outcome ~ treatment + mediator ", data= data[i,] , na.action =stats::na.omit)
          # b3 <- model.x$coefficients[2]
          model.m <- lqmm::lqm("mediator ~ . ", tau=tau, data= data[i,], na.action = stats::na.omit, control = lqm_control)
          b2 <- model.m$theta[2]
          # model.m <- lm("mediator ~ treatment ", data= data[i,] , na.action =stats::na.omit)
          # b2 <- model.m$coefficients[2]
          mediation_res <- b3 * b2
          return(mediation_res)
        }


        # calculate the mediation effect
        mediation_res <- boot::boot(data=df_outcome, statistic = mediation_fun, R=permutations, parallel = "no", ncpus=1)
        boot.res <- boot::boot.ci(mediation_res, type="bca")
        res$mediation_effect_ci_lower <- boot.res$bca[4]
        res$mediation_effect_ci_upper <- boot.res$bca[5]
        res$mediation_effect_significative <- !(boot.res$bca[4]<0 & boot.res$bca[5]>0)

        # Compute the proportion mediated
        prop_res <- boot::boot(data=df_outcome, statistic = prop_effect_fun, R=permutations, parallel = "no", ncpus=1)
        boot.res <- boot::boot.ci(prop_res, type="bca")
        res$prop_effect_ci_lower <- boot.res$bca[4]
        res$prop_effect_ci_upper <- boot.res$bca[5]
        res$prop_effect_significative <- !(boot.res$bca[4]<0 & boot.res$bca[5]>0)
        # calculate pvalue


        # Compute b1: Total Effect
        boot.res <- boot::boot.ci(b1_res, type="bca")
        res$total_effect_ci_lower <- boot.res$bca[4]
        res$total_effect_ci_upper <- boot.res$bca[5]
        res$total_effect_significative <- !(boot.res$bca[4]<0 & boot.res$bca[5]>0)

        # Compute b3: Direct Effect
        boot.res <- boot::boot.ci(b3_res, type="bca")
        res$direct_effect_ci_lower <- boot.res$bca[4]
        res$direct_effect_ci_upper <- boot.res$bca[5]
        res$direct_effect_significative <- !(boot.res$bca[4]<0 & boot.res$bca[5]>0)

        res$pvalue <- as.numeric(!(res$mediation_effect_significative & res$prop_effect_significative & res$total_effect_significative & res$direct_effect_significative))
        res$tau <- tau
        if(res$pvalue == 1)
          break
      },
      finally = {

      }
    )
  }
  res$family_test <- "mediation-quantreg"
  return(res)
}
