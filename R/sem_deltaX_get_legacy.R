#' Legacy R-loop implementation of derived-marker computation.
#'
#' Superseded by \code{\link{sem_deltaX_get}} (polars-native, ~14x faster on
#' ewas_data_hub). Kept for reference and possible A/B testing - call sites
#' have been migrated to the new \code{sem_deltaX_get()} (sem_core,
#' core_recover). This legacy function still writes per-sample bedgraph dumps
#' which the new path avoids.
#'
#' @return Invisible \code{NULL}. Side effect: writes derived pivot parquets
#'   and per-sample bedgraph files.
#'
#' @importFrom doRNG %dorng%
#'
#' @keywords internal
#' @noRd
sem_deltaX_get_legacy <- function()
{
  ssEnv <- core_get_session_info()
  sample_sheet <- sem_study_summary_get()
  keys <- ssEnv$keys_markers_figures
  # remove figure colunn
  keys <- keys[,-which(colnames(keys)=="FIGURE")]
  keys <- keys[,-which(colnames(keys)=="COMBINED")]
  keys <- unique(keys)

  area_position <-"POSITION"
  subarea_position <- "WHOLE"

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = seq_len(nrow(keys)))
  else
    progress_bar <- ""

  variables_to_export <- c("ssEnv", "io_dir_check_and_create", "subarea",
    "progress_bar","progression_index", "progression", "progressor_uuid",
    "owner_session_uuid", "trace","anno_probe_features_get", "localKeys",
    "io_file_path_build","%>%","core_get_session_info","core_log_event")
  sample_sheet_res <- data.frame()
  keys <- subset(keys,!is.na(SOURCE))
  keys <- subset(keys,!is.na(Q))
  keys <- subset(keys,(Q!=1))
  if(nrow(keys)==0)
    return()

  # sample_sheet_res <-  foreach::foreach(k =1:nrow(keys), .combine= cbind , .export = variables_to_export) %dorng%
  for ( k in seq_len(nrow(keys)))
  {
    sample_sheet_temp <- as.data.frame(sample_sheet[,c("Sample_ID","Sample_Group")])

    key <- keys[k,]
    source_marker <- key$SOURCE
    marker <- key$MARKER

    pivot_file_nameparquet_dest_hyper <- io_pivot_file_name_parquet(marker,"HYPER","POSITION","WHOLE")
    pivot_file_nameparquet_dest_hypo <- io_pivot_file_name_parquet(marker,"HYPO","POSITION","WHOLE")
    if(file.exists(pivot_file_nameparquet_dest_hyper) & file.exists(pivot_file_nameparquet_dest_hypo))
      next

    pivot_file_nameparquet <- io_pivot_file_name_parquet(source_marker,"HYPER",area_position,subarea_position)
    if (!file.exists(pivot_file_nameparquet)) {
      core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
        " [sem_deltaX_get] Source parquet not found, skipping marker=", marker,
        " source=", source_marker, " file=", pivot_file_nameparquet)
      next
    }
    pivot_hyper <- polars::pl$scan_parquet(pivot_file_nameparquet)
    positions_hyper <- as.data.frame(pivot_hyper$select(c("CHR","START","END"))$collect())
    pivot_hyper <- pivot_hyper$drop(c("CHR","START","END"))
    vector_shaped_hyper <- as.vector(as.matrix(as.data.frame(pivot_hyper$collect())))
    dim_pivot_hyper <- dim(as.data.frame(pivot_hyper$collect()))
    colname_pivot_hyper <- colnames(pivot_hyper)
    rm(pivot_hyper)
    vector_shaped_hyper[vector_shaped_hyper==0] <- NA
    length_hyper <- length(vector_shaped_hyper)

    pivot_file_nameparquet <- io_pivot_file_name_parquet(source_marker,"HYPO",area_position,subarea_position)
    if (!file.exists(pivot_file_nameparquet)) {
      core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
        " [sem_deltaX_get] Source parquet not found, skipping marker=", marker,
        " source=", source_marker, " figure=HYPO file=", pivot_file_nameparquet)
      next
    }
    pivot_hypo <- polars::pl$scan_parquet(pivot_file_nameparquet)
    positions_hypo <- as.data.frame(pivot_hypo$select(c("CHR","START","END"))$collect())
    pivot_hypo <- pivot_hypo$drop(c("CHR","START","END"))
    vector_shaped_hypo <- as.vector(as.matrix(as.data.frame(pivot_hypo$collect())))
    dim_pivot_hypo <- dim(as.data.frame(pivot_hypo$collect()))
    colname_pivot_hypo <- colnames(pivot_hypo)
    rm(pivot_hypo)
    length_hypo <- length(vector_shaped_hypo)
    vector_shaped_hypo[vector_shaped_hypo==0] <- NA

    vector_shaped <- c(vector_shaped_hyper,vector_shaped_hypo)
    rm(vector_shaped_hyper,vector_shaped_hypo)
    positions <- rbind(positions_hyper,positions_hypo)

    ####
    if(endsWith(marker,"P"))
      breaks_q <- as.numeric(key$Q)
    else if(endsWith(marker,"Q"))
      breaks_q <- unique(quantile(x = vector_shaped,probs=0:as.numeric(key$Q)/as.numeric(key$Q), na.rm =TRUE))

    vector_both <-  cut(vector_shaped, breaks=breaks_q, labels=FALSE)

    # cut(df$Age, quantile(df$Age, 0:4/4))
    # vector_both <- as.numeric(dplyr::ntile(x=vector_shaped , n= as.numeric(key$Q)))

    rm(vector_shaped)
    sem_save_figure(colname_pivot_hyper, dim_pivot_hyper,vector_both[seq_len(length_hyper)],positions_hyper,sample_sheet_temp,area_position,subarea_position,marker,"HYPER")
    sem_save_figure(colname_pivot_hypo, dim_pivot_hypo,tail(vector_both, length_hypo), positions_hypo,sample_sheet_temp,area_position,subarea_position,marker,"HYPO")

    if(ssEnv$showprogress)
      progress_bar(sprintf("Got: %s", stringr::str_pad(marker, 10, pad = " ")))

  }
}


sem_save_figure <- function(colname_pivot,dim_pivot,vector_figure,positions,sample_sheet,area, subarea,marker,figure)
{

  ssEnv <- core_get_session_info()
  vector_figure <- matrix(vector_figure, nrow=dim_pivot[1], ncol=dim_pivot[2])
  vector_figure <- as.data.frame(vector_figure)
  colnames(vector_figure) <- colname_pivot
  vector_figure <- cbind(positions,vector_figure)
  pivot_file_nameparquet <- io_pivot_file_name_parquet(marker,figure,area,subarea)
  vector_figure <- polars::as_polars_df(as.data.frame(vector_figure))
  vector_figure <- vector_figure$sort(c("CHR", "START"), descending = FALSE)
  vector_figure$write_parquet(pivot_file_nameparquet)

  vector_figure <- as.data.frame(vector_figure)

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 4:ncol(vector_figure))
  else
    progress_bar <- ""

  for(c in 4:(ncol(vector_figure)))
  {
    sample_id <- colnames(vector_figure)[c]
    bed_fname <- io_bed_file_name(sample_id,sample_sheet[sample_sheet$Sample_ID==sample_id,"Sample_Group"],marker,figure)
    if(ssEnv$showprogress)
      progress_bar(sprintf("Saving bed file %s of %s, %s.",sample_id,marker, figure))
    if(file.exists(bed_fname))
      next
    vector_figure_sample <- vector_figure[,c]
    sample_group <- sample_sheet[sample_sheet$Sample_ID==sample_id,"Sample_Group"]
    result <- cbind(positions,vector_figure_sample)
    # drop rows with NA
    result <- result[complete.cases(result),]
    result <- result[result[,4]>0,]
    bed_filename <- io_bed_file_name(sample_id,sample_group,marker,figure)
    io_dump_sample_as_bed_file(result,bed_filename)
  }
}
