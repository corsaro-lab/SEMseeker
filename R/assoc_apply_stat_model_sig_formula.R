assoc_apply_stat_model_sig_formula <- function (family_test, burdenValue, independent_variable, covariates)
{

  if(grepl("wilcoxon.paired",family_test) | grepl("t.test.paired",family_test))
  {
    predictors <- independent_variable
    sig.formula <- stats::as.formula(paste0(burdenValue,"~", predictors, sep=""))
  }

  #
  if(family_test=="wilcoxon" | family_test=="t.test" | family_test =="jsd" | family_test=="chisq.test" | family_test=="fisher.test" | family_test=="kruskal.test")
  {
    predictors <- independent_variable
    sig.formula <- stats::as.formula(paste0(burdenValue,"~", predictors, sep=""))
  }

  if(grepl("mean-permutation",family_test) | grepl("quantile-permutation",family_test) | grepl("spearman-permutation",family_test) )
  {
    predictors <- independent_variable
    sig.formula <- stats::as.formula(paste0(burdenValue,"~", predictors, sep=""))
  }

  if( family_test=="pearson" | family_test=="kendall" | family_test=="spearman")
  {
    predictors <- independent_variable
    sig.formula <- stats::as.formula(paste0(burdenValue,"~", predictors, sep=""))
  }

  if ( family_test=="multinomial" | family_test=="binomial")
  {
    # inversion of roles for variable
    if(is.null(covariates) || length(covariates)==0)
    {
      predictors <- burdenValue
    } else
    {
      predictors <- paste0(paste0(c(burdenValue, covariates),collapse="+", sep=""))
    }
    sig.formula <- stats::as.formula(paste0(independent_variable,"~", predictors, sep=""))
  }

  if(family_test=="gaussian" | family_test=="poisson" | grepl("quantreg", family_test) | grepl("polynomial",family_test) | grepl("exp",family_test)
    | grepl("pow10",family_test) | grepl("log10",family_test) | grepl("log",family_test) | grepl("mediation-ridge",family_test) )
  {
    if(is.null(covariates) || length(covariates)==0)
      predictors <- independent_variable
    else
      predictors <- paste0(paste0(c(independent_variable, covariates),collapse="+", sep=""))
    sig.formula <- stats::as.formula(paste0(burdenValue,"~", predictors, sep=""))
  }

  return (sig.formula)
}
