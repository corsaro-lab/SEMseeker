
sem_sample_group_check <- function(sample_sheet, signal_data)
{

  ssEnv <- core_get_session_info()

  # count rows of sample sheet and column of methylation
  # n_sample <- nrow(sample_sheet)
  # if ( nrow(unique(sample_sheet))!=n_sample)
  #   stop("")

  # check Sample_Group==Control has different Sample_ID of Sample_Group==Reference
  # check Sample_Group==Case has different Sample_ID of  Sample_Group==Reference

  sample_sheet <- as.data.frame(sample_sheet)
  sample_sheet <- sample_sheet[,!(colnames(sample_sheet) %in% core_name_cleaning(c("Probes_Count", ssEnv$keys$pasted)))]

  result <- NULL

  if((sum(sample_sheet$Sample_Group=="Case", na.rm=TRUE)<=3)
    || (sum(sample_sheet$Sample_Group=="Control", na.rm=TRUE)<=3)
    || (sum(sample_sheet$Sample_Group=="Reference", na.rm=TRUE)<=3))
    return(paste0("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Sample groups Case, Control and Reference must be present in the sample sheet, not enough samples each group. "))

    needColumns <- c("Sample_ID", "Sample_Group")
  missedColumns <- needColumns[!(needColumns %in% colnames(sample_sheet))]

  if (length(missedColumns) > 0) {
    result <- paste(result,  " Lost following columns ", missedColumns," ",format(Sys.time(), "%a %b %d %X %Y"), "Especting a column with name Sample_Group and possible values Reference,  Control and Case")
  }

  if(sum(is.na(sample_sheet$Sample_ID)))
  {
    result <- paste(result,  (" some samples have Sample_ID with NA !"))
  }

  if(sum(is.na(sample_sheet$Sample_Group)))
  {
    result <- paste(result,  (" some samples have Sample_Group with NA !"))
  }
  else
  {
    # Exist at least 3 samples per Sample_Group==Reference
    tabsample <- table(sample_sheet$Sample_Group)
    if(tabsample ["Reference"]<2)
    {
      result <- paste(result,  (" Required at least two samples for the Reference Sample Group"))
    }
  }


  if (sum(colnames(signal_data) %in% sample_sheet$Sample_ID)!=ncol(signal_data))
  {
    result <- paste(result, "\n", "Have a look: the methylation data first columns:", utils::head(colnames(signal_data), 4), "... some Sample_ID:", utils::head(sample_sheet$Sample_ID, 4), "\n" )
    lost_column <- colnames(signal_data)[!(colnames(signal_data) %in% sample_sheet$Sample_ID)]
    result <- paste(result, "\n","Lost column:", lost_column, "\n" )
    result <- paste(result, "\n", (" The methylation data has not the column's names as expected from the sample sheet in column Sample_ID!"))
  }
  populationGroups <- as.factor(sample_sheet$Sample_Group)
  # expected R as reference, S as Study, C as Control

  # reference population
  sample_sheet$Sample_Group <- R.utils::toCamelCase(tolower(sample_sheet$Sample_Group), capitalize=TRUE)
  sample_sheet$Sample_Group <- as.factor(sample_sheet$Sample_Group)
  matchedPopulation <- sort(as.character(levels(sample_sheet$Sample_Group))) %in% sort(as.vector(ssEnv$keys_sample_groups[,1]))
  if (is.element(FALSE, matchedPopulation)) {
    result <- paste(result,  " The Sample_Group should contain only: Reference, Control, Case" )
    result <- paste(result,  " The Sample_Group contains: ", paste(sort(as.character(levels(sample_sheet$Sample_Group))), collapse=", ") )
  }

  refence_group <- sample_sheet$Sample_ID[sample_sheet$Sample_Group=="Reference"]
  other_group <- sample_sheet$Sample_ID[sample_sheet$Sample_Group!="Reference"]
  if (sum(refence_group %in% other_group)==length(refence_group))
  {
    ssEnv$keys_sample_groups <-  data.frame("SAMPLE_GROUP"=c("Control","Case"))
  }

  refence_group <- sample_sheet$Sample_ID[sample_sheet$Sample_Group=="Reference"]
  other_group <- sample_sheet$Sample_ID[sample_sheet$Sample_Group!="Reference"]
  if (sum(refence_group %in% other_group)==length(refence_group))
  {
    ssEnv$keys_sample_groups <-  data.frame("SAMPLE_GROUP"=c("Control","Case"))
  }

  return(result)
}
