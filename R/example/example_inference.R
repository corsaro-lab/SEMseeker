result_folder <- file.path(getwd(),"/tmp/GSE186766")

# independent variable: the sample_sheet used by semseeker in the upstream process
# applicable regressions model: gaussian, poisson, binomial, quantreg
# applicable statistics test: wilcoxon, stats::t.test,
# appliocable correlation tests: pearson, kendall, spearman

# MUTATIONS_* ~ tcdd_mother + exam_age

# transformation_y to dependent variable: mutations and lesions: scale, log, log2, log10, exp, none, , quantile_n eg: quantile_4
# depth analysis
# 1: sample level
# 2: type level (gene, DMR, cpgisland) (includes 1)
# 3: genomic area: gene, body, gene tss1550, gene whole, gene tss200,  (includes 1 and 2)
# filter_p_value report after adjusting saves only significative nominal p-value
# "BODY","TSS1500","5UTR","TSS200","1STEXON","3UTR","EXNBND","WHOLE"
# "ISLAND","N_SHORE","S_SHORE","N_SHELF","S_SHELF","OPENSEA","WHOLE"

inference_details <- expand.grid("independent_variable"= c("tcdd_mother"),
                                "covariates"=c("exam_age+tcdd_father","breast_feeding"),
                                "family_test"=c("gaussian","wilcoxon"),
                                "transformation_y"="scale",
                                "scope"="SAMPLE",
                                "filter_p_value" = FALSE)


association_analysis(inference_details = inference_details)
