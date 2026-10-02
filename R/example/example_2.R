# install.packages("devtools")
# require("devtools")
# install_github("corsaro-lab/SEMseeker")
# require("semseeker")
# packageVersion('semseeker')
#
# devtools::install_github("hrbrmstr/dtupdate")
# require(dtupdate)
# github_update()
#
# workingFolder <- file.path("~/Downloads/GSE139307")
# dir.create(workingFolder)
#
# sample_sheet <- build_data_set_from_geo("GSE139307",workingFolder, 0)
#
# # ChAMP need the sample name variable as first column
# # so let's move the Sample_ID as first column
# sample_sheet <- data.frame("Sample_ID"=sample_sheet$Sample_ID, sample_sheet[, colnames(sample_sheet)!="Sample_ID"])
#
# write.table(
#   sample_sheet,
#   paste(workingFolder, "/", "final_samplesheet.csv", sep = ""),
#   row.names = FALSE,
#   sep=",",
#   quote = FALSE
# )
#
# require(ChAMP)
# idat_folder <- workingFolder
# result_folder = file.path( workingFolder,"/result/")
#
# myLoadN <- champ.load(directory = idat_folder,
#                       method = "minfi",
#                       methValue="B",
#                       autoimpute=TRUE,
#                       filterDetP=TRUE,
#                       ProbeCutoff=0,
#                       SampleCutoff=0.1,
#                       detPcut=0.01,
#                       filterBeads=TRUE,
#                       beadCutoff=0.05,
#                       filterNoCG=TRUE,
#                       filterSNPs=TRUE,
#                       population=NULL,
#                       filterMultiHit=TRUE,
#                       filterXY=TRUE,
#                       force=TRUE,
#                       arraytype="450K")
#
# # normalize with ChAMP
# normalizedData<-champ.norm(signal=myLoadN$signal,
#                            rgSet=myLoadN$rgSet,
#                            mset=myLoadN$mset,
#                            resultsDir= result_folder,
#                            method="SWAN",
#                            plotBMIQ=FALSE,
#                            arraytype="450K",
#                            cores= detectCores(all.tests = FALSE, logical = TRUE) - 1
# )
#
# sample_sheet <- utils::read.csv2("~/Downloads/GSE139307/final_samplesheet.csv")
# normalizedData <- readRDS("~/Downloads/GSE139307/normalizeddata.rds")
#
#
# # we need the Sample_Group, the sample sheet of the example has the column pathology we can use to distinguish Cases vs Controls
# sample_sheet$Sample_Group <- sample_sheet$dioxin.group.ch1
# sample_sheet[sample_sheet$Sample_Group=="CONTROL", "Sample_Group"] <- "Control"
# sample_sheet[sample_sheet$Sample_Group!="Control", "Sample_Group"] <- "Case"
#
# # we have not the reference Sample_Group so we reuse the Control sample group as reference
# reference <- subset(sample_sheet,Sample_Group=="Control")
# reference$Sample_Group <- "Reference"
# sample_sheet <- rbind(sample_sheet, reference)
#
#
#
# semseeker(input         = normalizedData,
#           sample_sheet  = sample_sheet,
#           result_folder = file.path(workingFolder, "/semseeker_result/"))
