sem_analyze_single_sample_both <- function( sample_detail, marker) {

  ssEnv <- core_get_session_info()
  start_time_single_sample <- Sys.time()
  # message(sample_detail$Sample_ID, " ", "SingleSampleBoth Sample analysis warmingUP ", Sys.time())
  result <- ""
  result <- result[-1]

  data_to_save <- data.frame("CHR"="", "START"="", "END"="")
  data_to_save <- data_to_save[-1,]
  marker <- as.character(marker)
  figures <- c("HYPER","HYPO")

  for( i in seq_along(figures))
  {
    figure <- figures[i]
    folder_to_save <- io_dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),paste0(marker,"_", figure, sep = "")))
    fileName <- io_file_path_build(folder_to_save,c(sample_detail$Sample_ID,marker,figure),"bed")
    if(file.exists(fileName))
    {
      data_to_saveTemp <- utils::read.table(fileName, sep="\t", col.names =c("CHR", "START", "END") )
      data_to_saveTemp$CHR <- anno_normalize_chr(data_to_saveTemp$CHR, "internal")
      data_to_save <- rbind(data_to_save, data_to_saveTemp )
    }
  }

  figure <- "BOTH"
  folder_to_save <- io_dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),paste0(marker,"_", figure, sep = "")))
  fileName <- io_file_path_build(folder_to_save,c(sample_detail$Sample_ID,marker,figure),"bed")

  io_dump_sample_as_bed_file(
    data_to_dump = data_to_save,
    fileName = fileName
  )

  result <- data.frame(var_name = if (!is.null(data_to_save)) nrow(data_to_save) else 0)

  end_time_single_sample <- Sys.time()
  time_taken <- end_time_single_sample - start_time_single_sample
  # message(sample_detail$Sample_ID, " ", "Completed sample ", time_taken)

  return(result)
}
