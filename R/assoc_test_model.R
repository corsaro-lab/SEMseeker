#' Statistical test model dispatcher
#'
#' @param family_test which family test to apply (e.g. "wilcoxon", "t.test", "kruskal.test", "pearson", "spearman")
#' @param tempDataFrame data frame to use with the test
#' @param sig.formula formula to apply
#' @param burdenValue name of the burden (dependent) column in tempDataFrame
#' @param independent_variable name of the independent variable column
#' @param transformation_y transformation applied to the dependent variable (for labelling plots)
#' @param plot logical; if TRUE, generate and save diagnostic box/scatter plots
#' @param samples_sql_condition SQL condition string used to filter samples (used for plot file naming)
#' @param key named list with AREA, SUBAREA, MARKER and FIGURE identifiers for this test
#'
#' @return A list (as a data.frame row) with test results including the p-value,
#'   test statistic, effect size, power, and model identifier; the exact fields
#'   depend on the chosen \code{family_test}.
#'
#' @section What the rank comparison reports, and the one number behind three names:
#' For \code{"wilcoxon"} three columns describe the same comparison and it is
#' worth knowing which is which, because two of them were once the same number.
#'
#' Vargha-Delaney A is the probability that a value drawn from the first group
#' exceeds one drawn from the second. It is therefore the area under the curve,
#' and it is also the Mann-Whitney statistic divided by the product of the group
#' sizes: \code{A}, \code{AUC} and \code{W/(n1 n2)} are one quantity under
#' three names, which was measured and not assumed.
#'
#' \describe{
#'   \item{\code{C_STATISTIC_AUC}}{That quantity, under the name the metric
#'     registry already carried for it, with \code{_CI_LOWER} and
#'     \code{_CI_UPPER} from \code{assoc_auc_confidence_interval()}.}
#'   \item{\code{RANK_BISERIAL_CORRELATION}}{\code{2A - 1}. Until 0.99.5 this
#'     column held \code{W/(n1 n2)}, that is A itself, so it reported the AUC
#'     under the name of a different statistic. The two are not a rescaling a
#'     reader can undo: on a worked example A was 0.2667 while \code{2A - 1} was
#'     -0.4667, the opposite sign.}
#'   \item{\code{effect_size_estimate}, \code{effect_size_magnitude}}{A again,
#'     and the qualitative label \code{effsize} attaches to it. The estimate
#'     duplicates \code{C_STATISTIC_AUC} on purpose, because consumers read it.}
#' }
#'
#' The \strong{sign} of \code{RANK_BISERIAL_CORRELATION} depends on which group
#' is first, that is on the order of the levels of the independent variable.
#' Since \code{inference_details$independent_variable_order} exists that order
#' is the request's to declare, so the sign is decided by the question rather
#' than by the alphabetical accident of the labels.
#'
#' @section The power is post-hoc and standardised:
#' \code{pwr::pwr.t2n.test()} takes Cohen's \code{d}, and both branches used to
#' hand it something else. The conversion now happens in
#' \code{assoc_power_from_d()}, in one place, and the two failures it ends are
#' worth stating because the column looks like evidence:
#'
#' \itemize{
#'   \item \code{"wilcoxon"} passed A, bounded in \code{[0, 1]}. At
#'     \code{A = 0.5}, the exact null, that reads as a medium effect and the
#'     reported power was 0.532 where the truth is the significance level. Most
#'     positions carry no effect, so the column sat above a half across the
#'     genome. It is now \code{d = sqrt(2) * qnorm(A)}.
#'   \item \code{"t.test"} passed the raw difference of the means, which carries
#'     the units of the data. The same comparison read on the beta and the
#'     M-value scale has \code{d = -1.50} either way and true power 1.000; the
#'     reported power was 0.053 and 0.171. A large effect read as no power, and
#'     the number moved with the scale, which is what standardising exists to
#'     prevent. It is now \code{effsize::cohen.d()}.
#' }
#'
#' Post-hoc power computed from the observed effect is a description of the
#' comparison that was run, not an argument that it was adequately sized. That
#' question is prospective and \code{assoc_statistical_power()} answers it.
#'
assoc_test_model <- function (family_test, tempDataFrame, sig.formula,burdenValue,independent_variable , transformation_y, plot , samples_sql_condition=samples_sql_condition, key)
{

  area <- as.character(key$AREA)
  subarea <- as.character(key$SUBAREA)
  marker <- as.character(key$MARKER)
  figure <- as.character(key$FIGURE)

  ssEnv <- core_get_session_info()
  res <- data.frame(pvalue=NA)
  if(family_test=="chisq.test")
  {
    # if (plot)
    #   plot_box_plot(tempDataFrame, independent_variable,burdenValue, transformation_y, family_test,samples_sql_condition, key)

    tempDataFrame <- as.data.frame(tempDataFrame)
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    dependent_variable <- dep_var[[2]]
    independent_variable <- dep_var[[3]] # sample_group
    sample1 <- round(tempDataFrame[,dependent_variable],3)
    sample2 <- tempDataFrame[,independent_variable]
    # create a contingency table
    contingency_table <- table( dependent_variable= sample1, independent_variable= sample2)

    result_chisq <- suppressWarnings(stats::chisq.test(as.matrix(contingency_table)))
    res$pvalue <- result_chisq$p.value
    res$r_model <- "stats_chisq.test"
    res$statistic_parameter <- result_chisq$statistic
    degrees_of_freedom <- result_chisq$parameter
    effect_size <- sqrt(result_chisq$statistic/nrow(tempDataFrame))
    res$effect_size <- effect_size
    power_result <- pwr::pwr.chisq.test(w = effect_size, N = nrow(tempDataFrame) , df = degrees_of_freedom, sig.level = as.numeric(ssEnv$alpha), power = )
    res$power <- power_result$power
  }

  if(family_test=="bartlett.test")
  {
    if (plot)
      plot_box_plot(tempDataFrame, independent_variable,burdenValue, transformation_y, family_test,samples_sql_condition, key)

    tempDataFrame <- as.data.frame(tempDataFrame)
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    dependent_variable <- dep_var[[2]]
    independent_variable <- dep_var[[3]] # sample_group

    bartlett_result <- suppressWarnings(stats::bartlett.test(x = tempDataFrame[, dependent_variable], g = tempDataFrame[, independent_variable]))
    res$pvalue <- bartlett_result$p.value
    res$r_model <- "bartlett.test"
    res$statistic_parameter <- bartlett_result$statistic
    degrees_of_freedom <- bartlett_result$parameter
    effect_size <- sqrt(bartlett_result$statistic/nrow(tempDataFrame))
    res$effect_size <- effect_size
    power_result <- pwr::pwr.chisq.test(w = effect_size, N = nrow(tempDataFrame) , df = degrees_of_freedom, sig.level = as.numeric(ssEnv$alpha), power = NULL)
    res$power <- power_result$power
  }

  if (family_test=="fisher.test")
  {
    if (plot)
      plot_box_plot(tempDataFrame, independent_variable,burdenValue, transformation_y, family_test,samples_sql_condition, key)

    tempDataFrame <- as.data.frame(tempDataFrame)
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    dependent_variable <- dep_var[[2]]
    independent_variable <- dep_var[[3]] # sample_group
    sample1 <- round(tempDataFrame[,dependent_variable],3)
    sample2 <- tempDataFrame[,independent_variable]
    # create a contingency table
    contingency_table <- table( dependent_variable= sample1, independent_variable= sample2)

    result_fisher <- suppressWarnings(stats::fisher.test(as.matrix(contingency_table)))
    res$pvalue <- result_fisher$p.value
    res$r_model <- "stats_fisher.test"
    res$statistic_parameter <- result_fisher$estimate
  }

  if (family_test=="kruskal.test")
  {
    if (plot)
      plot_box_plot(tempDataFrame, independent_variable,burdenValue, transformation_y, family_test, samples_sql_condition, key)

    tempDataFrame <- as.data.frame(tempDataFrame)
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    dependent_variable <- dep_var[[2]]
    independent_variable <- dep_var[[3]] # sample_group
    dependent_variable <- round(as.numeric(tempDataFrame[,dependent_variable]),3)
    group <- as.factor(tempDataFrame[,independent_variable])
    if(length(levels(group)) == 1)
    {
      res$pvalue <- NA
      res$r_model <- "stats_kruskal.test"
      res$kw_runk_sum <- 0
      return(res)
    }
    result_fisher <- suppressWarnings(stats::kruskal.test(x = dependent_variable, g = group))
    res$pvalue <- result_fisher$p.value
    res$r_model <- "stats_kruskal.test"
    res$kw_runk_sum <- result_fisher$statistic
    # pwr does not export pwr.kruskal.test; power calculation skipped for Kruskal-Wallis
    res$power <- NA_real_

    # pairwise wilcoxon test
    tryCatch({
      kw_result <- stats::pairwise.wilcox.test(dependent_variable, group)
    }, error = function(e) {
      return(res)
    })
    kruska_wallis_summary <- kw_result$p.value

    # for each group combination extract the p-value
    for (i in seq_len(nrow(kruska_wallis_summary))) {
      # i <-1
      for (j in seq_len(ncol(kruska_wallis_summary))) {
        # j <- 1
        p_value <- kruska_wallis_summary[i,j][1]
        row <- as.character(rownames(kruska_wallis_summary)[i])
        col <- as.character(colnames(kruska_wallis_summary)[j])
        if(i!=j)
          next
        pval_name <- paste0("PVALUE_KW_",as.character(col),"_",as.character(row),sep="")
        p_value <- data.frame(p_value)
        colnames(p_value) <- pval_name
        if (exists("res"))
          res <- cbind(res, p_value)
        else
          res <- data.frame(p_value)
      }
    }

    # remove rowname from res
    rownames(res) <- NULL

  }

  if(family_test=="jsd")
  {
    # Sample observations for the first and second sample
    tempDataFrame <- as.data.frame(tempDataFrame)
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    SPLIT <- split(round(tempDataFrame[,dep_var[[2]]],3), tempDataFrame[,dep_var[[3]]])
    sample1 <- SPLIT[[1]]
    sample2 <- SPLIT[[2]]

    # Combine the unique elements from both samples to create a common event space
    common_events <- unique(c(sample1, sample2))

    # Create adjusted frequency tables for both samples
    frequency_table1_adjusted <- tabulate(match(sample1, common_events), nbins = length(common_events))
    frequency_table2_adjusted <- tabulate(match(sample2, common_events), nbins = length(common_events))

    # Convert adjusted frequencies to probabilities
    probability_distribution1_adjusted <- frequency_table1_adjusted / sum(frequency_table1_adjusted)
    probability_distribution2_adjusted <- frequency_table2_adjusted / sum(frequency_table2_adjusted)

    # Calculate the Jensen-Shannon distance
    res$statistic_parameter <- suppressMessages(suppressWarnings(philentropy::JSD(rbind(probability_distribution1_adjusted, probability_distribution2_adjusted))))
    res$r_model <- "philentropy.JSD"
  }

  if(family_test=="wilcoxon")
  {

    if(length(levels(tempDataFrame[,independent_variable])) !=2)
    {
      res$pvalue <- NA
      res$r_model <- "stats_wilcox.test"
      res$Wilcox_Value <- 0
      return(res)
    }

    if (plot)
      plot_box_plot(tempDataFrame, independent_variable,burdenValue, transformation_y, family_test,samples_sql_condition, key)

    result_w  <- suppressWarnings(stats::wilcox.test(formula= sig.formula, data = as.data.frame(tempDataFrame), exact=TRUE))
    res$pvalue <- result_w$p.value

    #
    res[1,"Wilcox_Value"] <- result_w$statistic

    res$r_model <- "stats_wilcox.test"
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    SPLIT <- split(tempDataFrame[,dep_var[[2]]], tempDataFrame[,dep_var[[3]]])

    ## calculate the JSD
    sample1 <- SPLIT[[1]]
    sample2 <- SPLIT[[2]]

    # Combine the unique elements from both samples to create a common event space
    common_events <- unique(c(sample1, sample2))

    # Create adjusted frequency tables for both samples
    frequency_table1_adjusted <- tabulate(match(sample1, common_events), nbins = length(common_events))
    frequency_table2_adjusted <- tabulate(match(sample2, common_events), nbins = length(common_events))

    # Convert adjusted frequencies to probabilities
    probability_distribution1_adjusted <- frequency_table1_adjusted / sum(frequency_table1_adjusted)
    probability_distribution2_adjusted <- frequency_table2_adjusted / sum(frequency_table2_adjusted)

    # Calculate the Jensen-Shannon distance
    res$jsd <- suppressMessages(suppressWarnings(philentropy::JSD(rbind(probability_distribution1_adjusted, probability_distribution2_adjusted))))

    # Three names for one number, and two of the three were wrong about it.
    #
    # Vargha-Delaney A IS the area under the curve: the probability that a value
    # drawn from the first group exceeds one drawn from the second. W/(n1*n2) is
    # the same quantity again - measured, identical to VD.A - so
    # RANK_BISERIAL_CORRELATION held the AUC under the name of a different
    # statistic. The rank-biserial correlation is 2A - 1: on a worked example
    # A = 0.2667 while 2A - 1 = -0.4667, a different number with the opposite
    # sign, so a reader of that column read an association the statistic it is
    # named after calls negative.
    #
    # Its sign depends on which group is SPLIT[[1]], that is on the order of the
    # levels, which is why it is worth saying now: since independent_variable_order
    # exists the order is the request's to declare, so the sign of the effect is
    # decided by the question rather than by the alphabet.
    n1 <- length(SPLIT[[1]])
    n2 <- length(SPLIT[[2]])

    es_res <- effsize::VD.A(SPLIT[[1]], SPLIT[[2]])
    auc <- as.numeric(es_res$estimate)
    res$effect_size_estimate <- es_res$estimate
    res$effect_size_magnitude <- es_res$magnitude

    # The AUC under the name the metric registry already carried for it, which
    # until now had no producer, with the interval it was reported without.
    auc_ci <- assoc_auc_confidence_interval(auc, n1, n2, as.numeric(ssEnv$alpha))
    res$C_STATISTIC_AUC          <- auc
    res$C_STATISTIC_AUC_CI_LOWER <- auc_ci$lower
    res$C_STATISTIC_AUC_CI_UPPER <- auc_ci$upper

    res$RANK_BISERIAL_CORRELATION <- 2 * auc - 1

    # pwr wants Cohen's d and was handed A, which is bounded in [0, 1]: at
    # A = 0.5, the exact null, that reads as a medium effect and the reported
    # power was 0.53 where the truth is the significance level.
    res$power <- assoc_power_from_d(sqrt(2) * stats::qnorm(auc), n1, n2,
                                    as.numeric(ssEnv$alpha))
  }


  if(family_test=="t.test")
  {
    if (plot)
      plot_box_plot(tempDataFrame, independent_variable,burdenValue, transformation_y, family_test,samples_sql_condition, key)

    result_w  <-stats::t.test(formula= sig.formula, data = as.data.frame(tempDataFrame))
    res$pvalue <- result_w$p.value
    res$r_model <- "stats_t.test"
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    SPLIT <- split(tempDataFrame[,dep_var[[2]]], tempDataFrame[,dep_var[[3]]])
    res$statistic_parameter <- mean(SPLIT[[1]]) - mean(SPLIT[[2]])

    # The raw difference of the means carries the units of the data, and pwr
    # wants it standardised. Measured on the same comparison read on two scales,
    # Cohen's d = -1.50 either way and the true power 1.000: the reported power
    # was 0.053 on the beta scale and 0.171 on the M-value scale. A large effect
    # read as no power, and the number moved with the scale, which is the one
    # thing standardising exists to prevent.
    cohen_d <- tryCatch(
      as.numeric(effsize::cohen.d(SPLIT[[1]], SPLIT[[2]],
                                  pooled = TRUE, paired = FALSE, na.rm = TRUE)$estimate),
      error = function(e) NA_real_)
    res$power <- assoc_power_from_d(cohen_d, length(SPLIT[[1]]), length(SPLIT[[2]]),
                                    as.numeric(ssEnv$alpha))
  }

  if( family_test=="pearson" | family_test=="kendall" | family_test=="spearman")
  {
    tempDataFrame <- tempDataFrame[complete.cases(tempDataFrame),]
    result_cor <- stats::cor.test(as.numeric(tempDataFrame[,burdenValue]), as.numeric(tempDataFrame[,independent_variable]), method = as.character(family_test))
    res$pvalue <- result_cor$p.value
    res$r_model <- "stats_cor.test"
    statistic_parameter <- result_cor$estimate
    res$rho <- statistic_parameter
    power_result <- pwr::pwr.r.test(n = nrow(tempDataFrame) , r = statistic_parameter , sig.level = as.numeric(ssEnv$alpha) , power = NULL)
    res$power <- power_result$power

    if (plot)
    {
      chartFolder <- io_dir_check_and_create(ssEnv$result_folderChart,c("CORRELATION",core_name_cleaning(as.character(samples_sql_condition))))
      # plot a scatter plot of the burden value vs the independent variable and a linear regression line
      filename  <-  io_file_path_build(chartFolder,toupper(c(family_test,as.character(transformation_y), independent_variable,"Vs", burdenValue,area, subarea)),ssEnv$plot_format)
      grDevices::png(filename, width = 9, height = 9, units="in", res = as.numeric(ssEnv$plot_resolution_ppi))
      plot(tempDataFrame[,burdenValue], tempDataFrame[,independent_variable], xlab = burdenValue, ylab = independent_variable, main = paste("Correlation between", burdenValue, "and", independent_variable), pch = 19)
      abline(lm(tempDataFrame[,independent_variable] ~ tempDataFrame[,burdenValue]), col = "blue")
      dev.off()
    }
  }

  return (res)

}
