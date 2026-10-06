assoc_is_family_dicotomic <- function(family_test) {

  if(family_test=="gaussian" | family_test=="spearman" | family_test=="pearson" |
      family_test=="kendall" | grepl("quantreg-permutation", family_test) | grepl("quantreg", family_test)
    | family_test=="poisson"  | grepl("mean-permutation", family_test)  |
      grepl("polynomial", family_test) | grepl("exp", family_test)  | grepl("pow10", family_test) | grepl("log10", family_test) | grepl("log", family_test) | family_test=="spearman-permutation" | family_test=="quantile-permutation"
    | grepl("mediation-quantreg", family_test) | grepl("mediation-linear", family_test) | grepl("mediation-ridge", family_test)
    | grepl("^(limma|voom)_", family_test))  # regression families
    return(FALSE)

  if(family_test=="wilcoxon" | family_test=="t.test" | grepl("wilcoxon.paired",family_test) | grepl("t.test.paired",family_test) |  family_test =="jsd"
    | family_test=="kruskal.test" | family_test=="chisq.test" | family_test=="fisher.test" | family_test=="binomial" | family_test=="binomial_bulk" | family_test=="multinomial")
    return(TRUE)

  core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " family_test not recognized ")
  stop()
}
