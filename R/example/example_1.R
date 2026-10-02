# #
# # Example of semseeker application with a GEO data set
# #
# # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186766
#
# install.packages("devtools")
# require("devtools")
# install_github("corsaro-lab/SEMseeker", force=TRUE)
# require("semseeker")
# require(ChAMP)
#
#
# workingFolder <- file.path(getwd(),"/tmp/")
# dir.create(workingFolder)
#
# workingFolder <- file.path(getwd(),"/tmp/GSE186766")
# dir.create(workingFolder)
#
# sample_sheet <- semseeker::build_data_set_from_geo("GSE186766",workingFolder, 0)
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
#                       arraytype="EPIC")
#
# # normalize with ChAMP
# normalizedData<-champ.norm(signal=myLoadN$signal,
#                     rgSet=myLoadN$rgSet,
#                     mset=myLoadN$mset,
#                     resultsDir= result_folder,
#                     method="SWAN",
#                     plotBMIQ=FALSE,
#                     arraytype="EPIC",
#                     cores= detectCores(all.tests = FALSE, logical = TRUE) - 1
# )
#
#
# # we need the Sample_Group, the sample sheet of the example has the column pathology we can use to distinguish Cases vs Controls
# sample_sheet$Sample_Group <- sample_sheet$pathology.ch1
# sample_sheet[sample_sheet$Sample_Group=="AN", "Sample_Group"] <- "Case"
# sample_sheet[sample_sheet$Sample_Group=="control", "Sample_Group"] <- "Control"
#
# # we have not the reference Sample_Group so we reuse the Control sample group as reference
# reference <- subset(sample_sheet,Sample_Group=="Control")
# reference$Sample_Group <- "Reference"
# sample_sheet <- rbind(sample_sheet, reference)
#
#
# semseeker(input         = normalizedData,
#           sample_sheet  = sample_sheet,
#           result_folder = file.path(workingFolder, "/semseeker_result/"))
#
