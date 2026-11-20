io_source_data_get <- function(source_data, check_is_numeric=FALSE){

  # check id source data is a path or a dataframe
  if (!is.character(source_data))
    return(as.data.frame(source_data))

  if (!file.exists(source_data))
  {
    message("File ", source_data, " not found")
    return(NULL)
  }
  # check if ile extension is csv
  if ( grepl("\\.csv$", source_data))
    source <- withr::with_envvar(
      c(VROOM_CONNECTION_SIZE = "5000000"),
      as.data.frame(readr::read_csv2(source_data, show_col_types = FALSE))
    )

  if ( grepl("\\.xls$", source_data) |  grepl("\\.xlsx$", source_data))
  {
    source <- as.data.frame(readxl::read_excel(source_data,col_names = TRUE))
  }


  # check if the file is gz or rds
  if ( grepl("\\.gz$", source_data))
  {
    source <- as.data.frame(readr::read_delim(source_data, show_col_types = FALSE))
  }
  if ( grepl("\\.rds$", source_data) ||  grepl("\\.RDS$", source_data))
    source <- as.data.frame(readRDS(source_data))

  if ( grepl("\\.parquet$", source_data))
    source <- io_read_parquet_frame(source_data)

  # if exists a column PROBE, set it as rownames
  if ("PROBE" %in% colnames(source))
  {
    rownames(source) <- source[,"PROBE"]
    # remove colname PROBE
    source <- source[, -which(colnames(source) == "PROBE")]
  }

  if (!is.data.frame(source))
  {
    core_log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y") , " source_data ", source_data, " is missed !")
    stop()
  }

  # Normalize CHR to bare internal format if present
  if ("CHR" %in% colnames(source))
    source$CHR <- anno_normalize_chr(source$CHR, "internal")

  # check all values are numeric
  if (check_is_numeric & !all(vapply(source, is.numeric, logical(1))))
  {
    core_log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y") , " source_data ", source_data, " contains non-numeric values !")
    stop()
  }

  return(as.data.frame(source))
}
