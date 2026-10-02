#' Preliminary exploratory analysys of the data
#'
#' @param categorical_variables vector of variables or variables (with plus sign) to create exploratory pivot
#' @param numerical_variables vector of variables to have summary
#' @param sample_sheet path to the sample sheet
#' @param signal_data  path to the signal data
#' @param max_missed_sample_sheet max number of missing values in each sample sheet column
#' @param max_missed_signal_data max number of missing values in each signal data column and row
#' @param result_folder path to the result folder
#' @param ... other parameters to define options for semseeker
#' @param sample_id_column name of the column with sample ids
#' @param exploration_phase numeric value to preserve history of exploratory analysis
#' @param sample_sheet_mapping rules to rename and remove columns
#' @param values_mapping rules to recode values and remove samples (source missed leave blank in the mapping file)
#'
#' @param delete_keyword character. Value in \code{values_mapping} that marks
#'   samples for removal (default \code{"REMOVE"}).
#' @param mapping_folder character. Optional path to folder containing mapping
#'   files; overrides \code{sample_sheet_mapping} / \code{values_mapping} if
#'   provided (default \code{NULL}).
#' @param removal_folder character. Optional path to folder containing removal
#'   rule files (default \code{NULL}).
#' @param removal_rules character vector. Inline removal rules applied before
#'   mapping (default \code{c()}).
#' @return Invisibly \code{NULL}. Cleaned sample sheet and signal data are
#'   written to the result folder together with exploratory summary reports.
#' @examples
#' result_dir <- tempdir()
#' \dontrun{
#' util_exploratory_analysis(
#'   categorical_variables = c("Sample_Group", "Sex"),
#'   numerical_variables   = c("Age"),
#'   sample_sheet          = sample_sheet,
#'   signal_data           = beta_matrix
#' )
#' }
#' @export
#'
util_exploratory_analysis <- function(categorical_variables,numerical_variables, sample_sheet,signal_data, max_missed_sample_sheet = 0.3,
  max_missed_signal_data = 0.10,sample_id_column="Sample_ID", delete_keyword="REMOVE",
  exploration_phase="0",
  mapping_folder = NULL,
  sample_sheet_mapping = NULL,
  values_mapping = c(),
  removal_folder = NULL,
  removal_rules=c(),
  result_folder,  ... ) {

  core_init_env( result_folder= result_folder, ...)

  ssEnv <- core_get_session_info()
  core_log_event("BANNER:", format(Sys.time(), "%a %b %d %X %Y"), " SEMseeker will explore your data for project \n in ", ssEnv$result_folderData)
  folder_path <- io_dir_check_and_create(ssEnv$result_folderData,paste0("Exploratory_", exploration_phase))

  # remove Sample_Group and unique
  sample_sheet <- io_source_data_get(sample_sheet)
  sample_sheet <- sample_sheet[,!grepl("Sample_Group", colnames(sample_sheet))]
  sample_sheet <- unique(sample_sheet)

  step <- 0

  # rename columns as defined in the sample_sheet_mapping
  if (!is.null(sample_sheet_mapping)) {
    step <- step + 1
    map <- io_source_data_get(paste0(mapping_folder,sample_sheet_mapping))
    colun_to_remove <- map[map$Mapping == delete_keyword, "Variable"]
    sample_sheet <- sample_sheet[,!(colnames(sample_sheet) %in% colun_to_remove)]
    map <- map[map$Mapping != delete_keyword,]
    # rename columns
    for (i in seq_len(nrow(map)))
    {
      names_col <- colnames(sample_sheet)
      names_col[names_col == map$Variable[i]] <- map$Mapping[i]
      colnames(sample_sheet) <- names_col
    }
    file_path <- io_file_path_build(folder_path,c(step, "cleaned_sample_sheet") ,"csv")
    write.csv2(sample_sheet, file_path, row.names = FALSE)
  }

  # maps values in columns as defined in the values_mapping
  if(length(values_mapping)>0)
  {
    step <- step + 1
    for (i in seq_along(values_mapping))
    {
      map <- io_source_data_get(paste0(mapping_folder,values_mapping[i]))
      from <- colnames(map)[1]
      to <- colnames(map)[2]
      if (from==to)
      {
        core_log_event("ERROR: ",format(Sys.time(), "%a %b %d %X %Y"),"The source and target columns are the same")
        stop("The source and target columns are the same")
      }
      # random_string <- stringi::stri_rand_strings(1,10)
      random_string <- "Missing_Data"
      sample_sheet[is.na(sample_sheet[,from]),from] <- random_string
      sample_sheet[sample_sheet[,from]=="NA",from] <- random_string
      map[is.na(map[,from]),from] <- random_string
      map[is.na(map[,to]),to] <- random_string
      map[map[,from]=="NA",from] <- random_string
      map[map[,to]=="NA",to] <- random_string
      sample_sheet <- merge(sample_sheet,map,by.x=from,by.y=from,all.x=TRUE)
      # sample_sheet[sample_sheet[,to] == random_string,from] <- NA
      # remove where to is REMOVE
      sample_sheet <- sample_sheet[sample_sheet[,to] != delete_keyword,]
      # remove column of from
      sample_sheet <- sample_sheet[,!(colnames(sample_sheet)==from)]
    }
    file_path <- io_file_path_build(folder_path,c(step, "cleaned_sample_sheet") ,"csv")
    write.csv2(sample_sheet, file_path, row.names = FALSE)
  }

  # remove samples as defined in the removal_rules file
  if(length(removal_rules)>0)
    for (variable in removal_rules) {
      variable_splitted <- strsplit(variable,"\\+")
      variable_1 <- variable_splitted[[1]][1]
      variable_2 <- variable_splitted[[1]][2]
      file_path <- io_file_path_build(removal_folder, c(variable,"pivot"),"xlsx")
      rules <- io_source_data_get(file_path)
      # get first column
      first_var <- data.frame("variable_1"= rules[,variable_1])
      # other vars
      removal_details <- data.frame()
      for (i in 2:ncol(rules))
      {
        values <- data.frame( "values"=rules[,i])
        values$variable_2 <- colnames(rules)[i]
        tmp <- cbind(first_var,values)
        removal_details <- plyr::rbind.fill(removal_details, tmp)
      }
      removal_details <- subset(removal_details, removal_details$values == delete_keyword)
      for (r in seq_len(nrow(removal_details)))
      {      # remove from sample sheet
        selector <- (sample_sheet[,variable_1] %in% removal_details[r,1]) & (sample_sheet[,variable_2] %in% removal_details[r,3])
        sample_sheet <- subset(sample_sheet, !selector)
      }
    }

  if(ncol(sample_sheet) < 2)
  {
    core_log_event("ERROR: ",format(Sys.time(), "%a %b %d %X %Y")," Not enough columns found in the sample sheet")
    stop("Not enough column found in the sample sheet")
  }

  tryCatch({
    sample_sheet <- sample_sheet[!is.na(sample_sheet[,sample_id_column]),]
  }, error = function(e) {
    core_log_event("ERROR: ",format(Sys.time(), "%a %b %d %X %Y")," Sample ID column not found in the sample sheet")
    stop("Sample ID column not found in the sample sheet")
  })
  signal_data_original <- signal_data

  step <- step + 1
  # Run the function on your dataframe
  sample_sheet_report <- util_describe_dataframe(sample_sheet)
  # order desc by Missing_Values_Percent
  sample_sheet_report <- sample_sheet_report[order(-sample_sheet_report$Missing_Values_Percent),]
  file_path <- io_file_path_build(folder_path,c(step,"sample_sheet_report"),"csv")
  sample_sheet_report[sample_sheet_report=="0 (0%)"] <- " "
  write.csv2(sample_sheet_report, file_path, row.names = FALSE)
  io_save_latex_table(sample_sheet_report, file_path, paste0("Descriptive statistics of the sample sheet ", collapse = " "))

  # remove columns with all NA
  sample_sheet[, colSums(!is.na(sample_sheet)) > max_missed_sample_sheet * nrow(sample_sheet)]
  # replace NA with "NA" where column is not numeric
  # sample_sheet <- sample_sheet %>% dplyr::mutate_if(~!is.numeric(.), ~tidyr::replace_na(., "NA"))

  # check if the sample sheet is empty
  if(length(categorical_variables)>0)
    for (variable in categorical_variables) {

      # variable <- c("sample_type","simplified_disease_stage")
      # variable <- c("age")
      # create a pivot table using sample_type items as column and simplified_disease_stage items as row

      if(grepl("\\+",variable))
        variable <- strsplit(variable,split="\\+")[[1]]
      if (length(variable)>2)
      {
        core_log_event("ERROR: ",format(Sys.time(), "%a %b %d %X %Y")," max two variables managed!")
        next
      }
      if(length(variable)==2)
      {
        if(variable[1] %in% colnames(sample_sheet) && variable[2] %in% colnames(sample_sheet))
        {
          # replace NA with "NA" where column is variable[1]
          sample_sheet[,variable[1]] <- tidyr::replace_na(sample_sheet[,variable[1]], "NA")
          sample_sheet[,variable[2]] <- tidyr::replace_na(sample_sheet[,variable[2]], "NA")
          pivot <- as.data.frame.matrix(table(sample_sheet[,variable[1]],sample_sheet[,variable[2]]))
          first_col <-variable[1]
        }
        else
        {
          core_log_event("ERROR: ",format(Sys.time(), "%a %b %d %X %Y")," variable ",variable[1]," or ",variable[2]," not found in the sample sheet")
          next
        }
      }
      else
      {
        if(variable %in% colnames(sample_sheet))
        {

          sample_sheet[,variable] <- tidyr::replace_na(sample_sheet[,variable], "NA")
          # count grouping by variable
          pivot <- sample_sheet %>%
            dplyr::group_by(.data[[variable]]) %>%
            dplyr::summarise(count = dplyr::n(), .groups = "drop")

          pivot <- as.data.frame(pivot)
          pivot <- data.frame("stat"=pivot$count, row.names = pivot[,variable])
          # rownames(pivot) <- pivot[,variable]
          # # remove variable column
          # pivot <- pivot[,-1]
          # # pivot <- as.data.frame.matrix(matrix(table(sample_sheet[,variable])), col.names=variable, row.names = names(table(sample_sheet[,variable])))
          # # rename V1 as stat
          # colnames(pivot)[colnames(pivot)=="count"] <- "stat"
          first_col <-variable
        }
        else
        {
          core_log_event("ERROR: ",format(Sys.time(), "%a %b %d %X %Y")," variable ",variable," not found in the sample sheet")
          next
        }
      }
      if (nrow(pivot)==0)
        next
      # sum total
      pivot_total <- sum(pivot)
      pivot_order <- order(rowSums(pivot),decreasing = TRUE)
      if(ncol(pivot)==1)
      {
        pivot$TEMP <- rownames(pivot)
        pivot <- as.data.frame(pivot[pivot_order,])
        pivot <- data.frame("stat"=pivot$stat, row.names = pivot$TEMP)
      }
      else
      {

        pivot <- as.data.frame(pivot[pivot_order,])
      }
      if(ncol(pivot)>1)
        pivot_rowsums <- data.frame("Total"=paste0(rowSums(pivot), " (",round(rowSums(pivot)/pivot_total*100,2),"%)"))

      if(ncol(pivot)==1)
        pivot_columsums <- colSums(pivot)
      else
      {
        pivot_columsums <- t(data.frame(colSums(pivot)))
        colnames(pivot_columsums) <- colnames(pivot)
      }
      # calculate relative frequency
      pivot_rel_freq <- round(pivot/pivot_total*100,2)
      # build  a table with the total and relative frequency as percentage in the same cell
      for (i in seq_len(nrow(pivot_rel_freq))) {
        for (j in seq_len(ncol(pivot_rel_freq))) {
          pivot_rel_freq[i,j] <- paste0(pivot[i,j]," (",pivot_rel_freq[i,j],"%)")
        }
      }
      pivot_rel_freq[,first_col] <- rownames(pivot_rel_freq)
      if(ncol(pivot)==1)
        colnames(pivot_rel_freq)[colnames(pivot_rel_freq)=="pivot[pivot_order, ]"] <- "stat"
      pivot_rel_freq <- pivot_rel_freq[,c(ncol(pivot_rel_freq),seq_len(ncol(pivot_rel_freq)-1))]
      if(ncol(pivot)>1)
        pivot_rel_freq <- cbind(pivot_rel_freq, pivot_rowsums)
      if(ncol(pivot)>1)
      {
        total_by_row <- sum(pivot_columsums)[1]
        pivot_columsums_names <- colnames(pivot_columsums)
        pivot_columsums  <- as.data.frame(t(as.data.frame(paste0(pivot_columsums," (",round(pivot_columsums/pivot_total*100,2),"%)"))))
        colnames(pivot_columsums) <- pivot_columsums_names
        total_row <- as.data.frame(cbind("Col1"="Total",pivot_columsums,"Total"=total_by_row))
      }
      else
        total_row <-data.frame("Col1"="Total","stat"=sum(pivot_columsums)[1])
      # rename Col1
      colnames(total_row)[which(colnames(total_row)=="Col1")] <- first_col
      pivot_rel_freq <- plyr::rbind.fill(pivot_rel_freq, total_row)
      file_path <- io_file_path_build(folder_path, c(step, variable,"pivot"),"csv")
      # replace 0 (0%) with null string
      pivot_rel_freq[pivot_rel_freq=="0 (0%)"] <- " "
      write.csv2(pivot_rel_freq, file_path, row.names = FALSE)
      io_save_latex_table(pivot_rel_freq, file_path, paste0("Descriptive statistics for ",variable,sep=" ", collapse = " "))
    }

  if(length(numerical_variables)>0)
  {
    num_res_data <- data.frame()
    for (variable in numerical_variables) {
      if (variable %in% colnames(sample_sheet) == FALSE)
      {
        core_log_event("ERROR: ",format(Sys.time(), "%a %b %d %X %Y")," variable ",variable," not found in the sample sheet")
        next
      }
      # calculate mean and standard deviation
      mean <- paste(round(mean(sample_sheet[,variable], na.rm = TRUE),2) ," \u00b1",round(sd(sample_sheet[,variable], na.rm = TRUE),2))
      # calculate median
      median <- round(stats::median(sample_sheet[,variable], na.rm = TRUE))
      # calculate min and max
      min <- round(min(sample_sheet[,variable], na.rm = TRUE),2)
      max <- round(max(sample_sheet[,variable], na.rm = TRUE),2)
      # calculate quantiles
      quantiles <- round(quantile(sample_sheet[,variable], c(0.25,0.5,0.75), na.rm = TRUE),2)
      # calculate number of missing values
      missing_values <- paste(round(sum(is.na(sample_sheet[,variable]))/length(sample_sheet[,variable])*100,2),"%")
      # build a table with the descriptive statistics
      res_data <- data.frame(Statistic=c("Mean \u00b1 Std","Median","Min","Max","1st Quartile","2nd Quartile","3rd Quartile","Missing Values"),
        Value=c(mean,median,min,max,quantiles,missing_values))
      res_data <- as.data.frame(t(res_data))
      colnames(res_data) <- res_data[1,]
      res_data <- res_data[-1,]
      res_data <- data.frame("Variable"=variable, res_data)
      num_res_data <- plyr::rbind.fill(num_res_data, res_data)
    }
    if(nrow(num_res_data)==0)
      return()
    file_path <- io_file_path_build(folder_path,c(step, "numeric_variable_pivot"),"csv")
    num_res_data[num_res_data=="0 (0%)"] <- " "
    write.csv2(num_res_data, file_path, row.names = FALSE)
    io_save_latex_table(num_res_data, file_path, "Descriptive statistics for numerioc variables" )
  }

  signal_data <- io_source_data_get(signal_data, TRUE)
  # AI-224: normalise BOTH sides. Cleaning only the signal columns and then
  # comparing them with the raw sample sheet identifiers is what produced the
  # empty intersections below (and an all-but-PROBE "cleaned" parquet) for
  # identifiers such as "C3L-00001-06".
  .normalized <- core_normalize_sample_ids(sample_sheet, signal_data,
                                           sample_id_column = sample_id_column)
  sample_sheet <- .normalized$sample_sheet
  signal_data  <- .normalized$signal_data
  rm(.normalized)

  # signal_data <- signal_data[1:1000,]
  gc()

  ssEnv <- core_get_meth_tech(signal_data)
  signal_data <- is.na(signal_data)
  # rm(signal_data)
  nrow_ex_ante <- nrow(signal_data)
  n_missed_per_row <- as.vector(rowSums(signal_data))
  n_na <- paste0(round((100*(sum(n_missed_per_row)/nrow_ex_ante)/ncol(signal_data)),2), "%", sep="")
  row_selector <- (n_missed_per_row > max_missed_signal_data * ncol(signal_data))
  nrow_to_delete <- sum(row_selector)

  starting_col <- ncol(signal_data)
  n_missed_per_col <- colSums(signal_data[!row_selector,])
  ncol_to_delete <- sum(as.vector(n_missed_per_col > max_missed_signal_data * nrow(signal_data)))

  remaining_rows <- ((nrow_ex_ante - nrow_to_delete)/nrow_ex_ante)/starting_col
  remaining_cols <-  (starting_col - ncol_to_delete)
  remaining_values <- paste0(round(100*remaining_rows*remaining_cols,2), "%")

  samples_in_sample_sheet_not_in_signal_data <- nrow(sample_sheet[!(c("PROBE",sample_sheet[, sample_id_column]) %in% colnames(signal_data)),])
  samples_in_signal_data_not_in_sample_sheet <- sum(!(colnames(signal_data) %in% c("PROBE",sample_sheet[, sample_id_column])))

  signal_data_report <- data.frame(
    "Missing values"=n_na,
    "Number of Positions"=nrow_ex_ante,
    "Number of Samples"=ncol(signal_data),
    "Number of Positions to delete"= paste0(nrow_to_delete, " (", round(100*nrow_to_delete/nrow_ex_ante,2), "%)"),
    "Number of Samples to delete"= paste0(ncol_to_delete, " (", round(100*ncol_to_delete/ncol(signal_data),2), "%)"),
    "Remaining Positions"= paste0(nrow_ex_ante - nrow_to_delete, " (", round(100*(nrow_ex_ante - nrow_to_delete)/nrow_ex_ante,2), "%)"),
    "Remaining Samples"= paste0(remaining_cols, " (", round(100*remaining_cols/starting_col,2), "%)"),
    "Remaining values"= remaining_values,
    "Samples in sample sheet not in signal data"=samples_in_sample_sheet_not_in_signal_data,
    "Samples in signal data not in sample sheet"=samples_in_signal_data_not_in_sample_sheet,
    "Tech"=ssEnv$tech,
    "Is Beta"=ssEnv$beta)

  file_path <- io_file_path_build(folder_path, c(step,"signal_data_report"),"csv")
  signal_data_report[signal_data_report=="0 (0%)"] <- " "
  write.csv2(signal_data_report, file_path, row.names = FALSE)
  io_save_latex_table(signal_data_report, file_path, "Signal data summary" )

  step <- step + 1
  rm(signal_data)
  gc()
  signal_data <- io_source_data_get(signal_data_original)
  # AI-224: same normalisation as above - sample_sheet identifiers were already
  # normalised, so only the freshly reloaded matrix needs it here.
  signal_data <- core_normalize_sample_ids(sample_sheet = NULL,
                                           signal_data = signal_data,
                                           sample_id_column = sample_id_column)$signal_data
  # signal_data <- signal_data[1:1000,]

  # remove rows with too many missing values
  signal_data <- signal_data[!row_selector,]
  # remove columns with too many missing values
  signal_data <- signal_data[,n_missed_per_col <= max_missed_signal_data * nrow(signal_data)]
  signal_data$PROBE <- rownames(signal_data)
  signal_data <- as.data.frame(signal_data[, colnames(signal_data) %in% c("PROBE",sample_sheet[, sample_id_column])])
  signal_data <- signal_data[,c("PROBE",colnames(signal_data)[seq_len(ncol(signal_data)-1)])]
  # save cleaned signal data
  file_path <- io_file_path_build(folder_path, c(step,"cleaned_signal_data"),"parquet")
  polars::as_polars_df(as.data.frame(signal_data))$write_parquet(file_path)

  # save cleaned sample_sheet
  sample_sheet <- sample_sheet[sample_sheet[,sample_id_column] %in% colnames(signal_data),]
  file_path <- io_file_path_build(folder_path, c(step,"cleaned_sample_sheet") ,"csv")
  write.csv2(sample_sheet, file_path, row.names = FALSE)


}
