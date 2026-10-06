io_bed_file_name <- function(sample_id, sample_group, marker, figure,
                          skip_dir_create = FALSE)
{
  ssEnv <- core_get_session_info()
  if (is.na(sample_group) || sample_group=="")
  {
    stop("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " io_bed_file_name: sample_group is empty")
  }

  # when the caller has already ensured the destination dir exists
  # (e.g. sem_analyze_population creates all combos once at the top of the
  # function), skip the per-call file.exists+dir.create check. On a 4000-
  # sample population with 4 marker/figure combos this saves ~16k stat
  # syscalls per pipeline pass.
  if (skip_dir_create) {
    folder_to_save <- file.path(ssEnv$result_folderData,
                                as.character(sample_group),
                                paste0(marker, "_", figure))
  } else {
    folder_to_save <- io_dir_check_and_create(ssEnv$result_folderData,
                                           c(as.character(sample_group),
                                             paste0(marker, "_", figure)))
  }
  bed_ext <- unique(ssEnv$keys_markers_figures_default[ ssEnv$keys_markers_figures_default$MARKER==marker & ssEnv$keys_markers_figures_default$FIGURE==figure, "EXT"])
  if(length(bed_ext)==0){
    stop("ERROR: io_bed_file_name: bed_ext is empty")
  }
  file_name <- io_file_path_build(folder_to_save,c(sample_id,marker,figure),bed_ext, add_gz=TRUE)
  return(file_name)
}
