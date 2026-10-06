#' @param observed_probes character. Probe identifiers present in the input.
#' @param keys data.frame of AREA / SUBAREA (optionally COMBINED) to report on.
#'   Defaults to the areas selected for the run. The mandatory coverage gate passes the FULL
#'   annotation set here: `areas` defaults to POSITION only
#'   (util_keys_create.R:83) and this function skips POSITION and PROBE, so on
#'   a default run the session keys are empty and no chart was ever produced.
#'   The coverage picture must not depend on which areas the run analyses.
#' @keywords internal
#' @noRd
sem_coverage_analysis <- function(observed_probes, keys = NULL)
{
  # probe_features <- PROBES_CHR_CHR
  # area <- c("CHR")
  # probes_prefix <- "PROBES_CHR_"
  ssEnv <- core_get_session_info()
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Started Coverage analysis.")
  if (is.null(keys))
    keys <- ssEnv$keys_areas_subareas
  if (!is.null(keys) && nrow(keys) > 0 && !("COMBINED" %in% colnames(keys)))
    keys$COMBINED <- apply(keys[, c("AREA", "SUBAREA")], 1, function(x)
      core_name_cleaning(gsub(" ", "", paste0(x[x != ""], collapse = "_"))))
  keys <- keys[keys$AREA != "POSITION",]
  keys <- keys[keys$AREA != "PROBE",]
  if(plyr::empty(keys))
  {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " No keys found for coverage analysis.")
    return()
  }
  # local accumulators. The previous exists()-based logic resolved
  # through globalenv(), so an object of the same name left in an interactive
  # session changed the result - same defect class fixed in
  # sem_study_summary_total() by the fix for empty per-sample burden columns. Now that coverage runs on EVERY SEM
  # analysis this function must be deterministic.
  cov_result <- NULL
  tot_result <- NULL

  for ( k in seq_len(nrow(keys)))
  {
    # k <- 16
    subarea <- as.character(keys[k,"SUBAREA"])
    area <-  as.character(keys[k,"AREA"])
    area_subarea <- as.character(keys[k,"COMBINED"])

    if(!grepl("_",area_subarea))
      area_subarea <- paste(area_subarea,"_","WHOLE",sep="")

    # core_log_event(area,"\n")
    covered_count <- NULL

    probe_features <- anno_probe_features_get(area_subarea)
    if(plyr::empty(probe_features) | nrow(probe_features)==0)
      next

    total_count <- stats::aggregate(probe_features$PROBE, list(probe_features[, area_subarea]), FUN=length)
    colnames(total_count) <- c(as.character(area_subarea),"COUNT_TOTAL")

    probe_filtered <- probe_features[ probe_features$PROBE %in% observed_probes,]
    # if no probe is covered
    if (nrow(probe_filtered)>0)
    {
      covered_count <- stats::aggregate(probe_filtered$PROBE,list(probe_filtered[,area_subarea]), FUN=length)
      colnames(covered_count) <- c(as.character(area_subarea),"COUNT_COVERED")
      coverage <- merge(total_count, covered_count, by=area_subarea, all.x = TRUE)
    }
    else
    {
      coverage <- total_count
      if(nrow(total_count)>0)
        coverage$COUNT_COVERED <- 0
      # coverage <- merge(total_count, covered_count, by=area_subarea, all.x = TRUE)
    }

    coverage[is.na(coverage)] <- 0
    if(nrow(coverage)<2)
      next

    coverage <- coverage[ !is.na(coverage[,area_subarea]),]
    coverage <- coverage[ coverage[,area_subarea]!="",]

    coverage$PERC <- plyr::round_any(100 * as.numeric(coverage$COUNT_COVERED) / coverage$COUNT_TOTAL, 5)
    cov_stat <- stats::aggregate(coverage$PERC, list(coverage$PERC), FUN=length)
    colnames(cov_stat) <- c("COV_PERC","COUNT")

    cov_stat$AREA <- area
    cov_stat$SUBAREA <- subarea
    if(!is.null(cov_result))
      cov_result <- rbind(cov_result, cov_stat)
    else
      cov_result <- cov_stat

    colnames(total_count) <- c("AREA","COUNT_TOTAL")
    total_count$AREA <- area
    total_count$SUBAREA <- subarea

    if(!is.null(tot_result))
      tot_result <- rbind(total_count, tot_result)
    else
      tot_result <- total_count

  }

  # No area produced a usable table (e.g. tech without Illumina annotation):
  # nothing to plot. Previously this line dereferenced a possibly non-existent
  # tot_result and raised "object 'tot_result' not found".
  if(is.null(tot_result) || nrow(tot_result)<1 || is.null(cov_result))
    return()

  chartFolder <- io_dir_check_and_create(ssEnv$result_folderChart,"COVERAGE")
  filename  <-  io_file_path_build(chartFolder,c("EACH_AREA_COVERAGE_ANALYSIS"),ssEnv$plot_format)

  temp_cov_result <- subset(cov_result,cov_result$SUBAREA !="WHOLE")
  if(nrow(temp_cov_result)>2)
  {
    number_areas <- stats::aggregate(temp_cov_result$COUNT, list(temp_cov_result$AREA,temp_cov_result$SUBAREA), FUN=sum)
    colnames(number_areas) <- c("AREA","SUBAREA","TOTAL_AREAS")
    temp_cov_result <- merge(temp_cov_result, number_areas, by=c("AREA","SUBAREA"), all.x = TRUE)
    temp_cov_result$GENOMIC_AREA <- paste(temp_cov_result$AREA, temp_cov_result$SUBAREA, sep=" ")
    temp_cov_result$AREA_PERC <- round(100* temp_cov_result$COUNT / temp_cov_result$TOTAL_AREAS,2)

    ggp <- ggplot2::ggplot(temp_cov_result, ggplot2::aes(.data$GENOMIC_AREA, .data$COV_PERC)) +    # Create default ggplot2 heatmap
      ggplot2::geom_tile(ggplot2::aes(fill = .data$AREA_PERC)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggplot2::geom_text(ggplot2::aes(label = .data$AREA_PERC), size = 3) +
      ggplot2::scale_fill_gradient(low = "white", high = "#1b98e0") +
      ggplot2::labs(x = "Genomic Area", y = " Percentage of covered Probes") +
      ggplot2::labs(fill  = "Percentage\nover\neach\narea") +
      ggplot2::scale_y_continuous(breaks = seq(0,100, by=5))

    invisible
    (
      ggplot2::ggsave(
        filename,
        plot = ggp,
        # scale = 1,
        # width = 1240,
        # height = 1240,
        # units = c("px"),
        dpi = as.numeric(ssEnv$plot_resolution_ppi)
      )
    )
  }


  temp_cov_result <- subset(cov_result, cov_result$SUBAREA !="WHOLE" & cov_result$SUBAREA !="CHR" & cov_result$SUBAREA !="PROBE" & cov_result$AREA != "PROBE")
  if(nrow(temp_cov_result>2))
  {
    total_areas <- sum(number_areas$TOTAL_AREAS)
    temp_cov_result$GENOMIC_AREA <- paste(temp_cov_result$AREA, temp_cov_result$SUBAREA, sep=" ")
    temp_cov_result$AREA_PERC <- round(100* temp_cov_result$COUNT / total_areas,2)
    filename  <-  io_file_path_build(chartFolder,c("WHOLE_GENOME_COVERAGE_ANALYSIS"),ssEnv$plot_format)
    # temp_cov_result$COV_PERC <- sprintf('%02d', str_pad(temp_cov_result$COV_PERC, 3, pad = "0"))
    temp_cov_result$COV_PERC <- sprintf('%03d', temp_cov_result$COV_PERC)
    temp_cov_result$AREA_PERC <-  round(temp_cov_result$AREA_PERC,2)
    h_total <- temp_cov_result %>%
      dplyr::group_by(.data$GENOMIC_AREA) %>%
      dplyr::summarise(AREA_PERC = sum(.data$AREA_PERC)) %>%
      dplyr::mutate(COV_PERC = 'Total')
    v_total <- temp_cov_result %>%
      dplyr::group_by(.data$COV_PERC) %>%
      dplyr::summarise(AREA_PERC = sum(.data$AREA_PERC)) %>%
      dplyr::mutate(GENOMIC_AREA = 'Total')

    ggp <- ggplot2::ggplot(data = temp_cov_result, ggplot2::aes( x =.data$GENOMIC_AREA, y=.data$COV_PERC )) +    # Create default ggplot2 heatmap
      ggplot2::geom_tile(ggplot2::aes(fill = .data$AREA_PERC)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggplot2::geom_text(ggplot2::aes(label = .data$AREA_PERC), size = 3) +
      ggplot2::scale_fill_gradient(low = "white", high = "#1b98e0") +
      ggplot2::labs(x = "Genomic Area", y = " Percentage of covered Probes") +
      ggplot2::labs(fill  = "Percentage\nover\nwhole\ngenome", color="% Covered\narea\nover\nstudied") +
      ggplot2::geom_point(data = h_total, ggplot2::aes(color = .data$AREA_PERC), size = 8, shape = 19) +
      ggplot2::geom_point(data = v_total, ggplot2::aes(color = .data$AREA_PERC), size = 8, shape = 19) +
      ggplot2::scale_color_gradient2(low = "white",high = "grey",midpoint = 0) +
      ggplot2::geom_text(data = h_total, size = 3, ggplot2::aes(label = .data$AREA_PERC)) +
      ggplot2::geom_text(data = v_total, size = 3, ggplot2::aes(label = .data$AREA_PERC))

    invisible(
      ggplot2::ggsave(
        filename,
        plot = ggp,
        # scale = 1,
        # width = 1240,
        # height = 1240,
        # units = c("px"),
        dpi = as.numeric(ssEnv$plot_resolution_ppi)
      )
    )

  }

  #show total count of probe_features per area
  # temp_cov_result <- subset(cov_result, SUBAREA !="WHOLE" & SUBAREA !="CHR" & SUBAREA !="PROBE")
  # temp_cov_result$GENOMIC_AREA <- paste(temp_cov_result$AREA, temp_cov_result$SUBAREA, sep=" ")
  # temp_cov_result$AREA_PERC <- round(100* temp_cov_result$COUNT / total_areas,2)
  # filename = paste0( chartFolder,"/","WHOLE_GENOME_COVERAGE_ANALYSIS_ABSOLUTE_COUNT.png",sep="")
  # h_total <- temp_cov_result %>%
  #   dplyr::group_by(.data$GENOMIC_AREA) %>%
  #   dplyr::summarise(AREA_COUNT = sum(.data$COUNT_COVERED)) %>%
  #   dplyr::mutate(COUNT_TOTAL = 'Total')
  # v_total <- temp_cov_result %>%
  #   dplyr::group_by(.data$COV_PERC) %>%
  #   dplyr::summarise(AREA_COUNT = sum(.data$COUNT_COVERED)) %>%
  #   dplyr::mutate(GENOMIC_AREA = 'Total')
  #
  # ggp <- ggplot2::ggplot(data = temp_cov_result, ggplot2::aes( x =.data$GENOMIC_AREA, y=.data$COUNT_COVERED )) +    # Create default ggplot2 heatmap
  #   ggplot2::geom_tile(ggplot2::aes(fill = .data$AREA_COUNT)) +
  #   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #   ggplot2::geom_text(ggplot2::aes(label = .data$AREA_COUNT), size=3) +
  #   ggplot2::scale_fill_gradient(low = "white", high = "#1b98e0") +
  #   ggplot2::labs(x = "Genomic Area", y = " Absolute count of covered Probes") +
  #   ggplot2::labs(fill  = "Absolute Count\nover\nwhole\ngenome", color="% Absolute Count \narea\nover\nstudied") +
  #   ggplot2::geom_point(data = h_total, ggplot2::aes(color = .data$AREA_COUNT), size = 10, shape = 19) +
  #   ggplot2::geom_point(data = v_total, ggplot2::aes(color = .data$AREA_COUNT), size = 10, shape = 19) +
  #   ggplot2::scale_color_gradient2(low = "white",high = "grey",midpoint = 0) +
  #   ggplot2::geom_text(data = h_total, size = 3, ggplot2::aes(label = .data$AREA_COUNT)) +
  #   ggplot2::geom_text(data = v_total, size = 3, ggplot2::aes(label = .data$AREA_COUNT))
  #
  # ggplot2::ggsave(
  #   filename,
  #   plot = ggp,
  #   scale = 1,
  #   width = 1240,
  #   height = 1240,
  #   units = c("px"),
  #   dpi = as.numeric(ssEnv$plot_resolution_ppi)
  # )
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Coverage analysis executed." )

  # cov_result <- reshape2::dcast(data = cov_result, AREA + SUBAREA ~ COV_PERC, value.var = "COUNT", sum)
  # number_areas <- rowSums(cov_result[,3:ncol(cov_result)])
  # cov_result[,3:ncol(cov_result)] <- round(100*cov_result[,3:ncol(cov_result)]/number_areas,0)
  # chartFolder <- io_dir_check_and_create(ssEnv$result_folderChart,"COVERAGE")
  # tt <- as.data.frame(cov_result[,3:ncol(cov_result)])
  # rownames(tt) <- paste(cov_result$AREA,cov_result$SUBAREA, sep=" ")
  # colors <- grDevices::hsv(0.560, seq(0,1,length.out = 100) , 1)
  # if (!plyr::empty(tt))
  #   if(nrow(tt) > 2 & ncol(tt) > 2)
  #   {
  #     filename = paste0( chartFolder,"/","COVERAGE_ANALYSIS.png",sep="")
  #     grDevices::png(file= filename, width=1240, height = 1240, pointsize = 15, res = ssEnv$plot_resolution)
  #     stats::heatmap(
  #       x =  as.matrix(tt),
  #       col = colors,
  #       scale = "column",
  #       Colv = NA,
  #       Rowv = NA,
  #       # RowSideColors = as.vector( as.character(cov_result$COV_PERC)),
  #       margins = c(15, 15),
  #       main = "Coverage analysis of probe_features for each genomic area"
  #     )
  #     graphics::legend(x="right", legend=c("min 0%", "mean 50%", "max 100%"), fill=colors[c(1,50,100)])
  #     grDevices::dev.off()
  #   }
  #
  #
  # tot_result <- subset(tot_result, SUBAREA!="CHR" & SUBAREA != "WHOLE")
  # # tot_result <- aggregate(tot_result$SUBAREA, list(tot_result$COUNT_TOTAL), FUN=length)
  # tot_result <- reshape2::dcast(data = tot_result, AREA + SUBAREA ~ COUNT_TOTAL, value.var = "COUNT_TOTAL", length)
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Completed Coverage analysis.")

  return(cov_result)
}
