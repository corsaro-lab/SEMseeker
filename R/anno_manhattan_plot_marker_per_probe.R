#' Title anno_manhattan_plot_marker_per_probe
#'
#' @param max_sample max number of samples to plot
#' @param min_sample min number of samples to plot
#' @param min_signal_probe min signal value of the probe to be plotted
#' @param label_font_size size of the labels
#' @param hyper_color color to assign to hypermethylated probes
#' @param hypo_color color to assign to hypomethylated probes
#' @param non_outlier_color color to assign to probes that are not outliers
#' @param limit_label_color color to assign to the labels of the limit lines
#' @param limit_line_color color to assign to the limit lines
#' @param limit_line_color_median color to assign to the median limit line
#' @param reference_samples_color color to assign to the reference samples
#' @param result_folder foder where the results are stored
#' @param probe_name_max  cg name of the probe to represent tyh probe with maximum burden
#' @param probe_name_min cg name of the probe to represent tyh probe with miniumal burden
#' @param show_labels show labels in the plot
#' @param parallel_strategy strategy to use for parallelization
#' @param ... other parameter
#' @param maxResources percentage of max system's resource to use
#'
#' @return Invisibly \code{NULL}. The function saves PNG plot files to disk as a
#'   side effect.
anno_manhattan_plot_marker_per_probe <- function(probe_name_max = "cg11680158", probe_name_min = "cg11680158", max_sample=0, min_sample=0 , min_signal_probe=0, label_font_size=3,
  hyper_color = "blue", hypo_color = "orange",  non_outlier_color = "grey",  limit_label_color = "blue",
  limit_line_color = "red", limit_line_color_median = "black",
  reference_samples_color = "cyan", show_labels = FALSE,
  result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, ...)
{

  ssEnv <- core_init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE,
    figures=c("HYPER","HYPO"), areas="PROBE", ...)

  core_log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), " SEMseeker will generate images for significative probes \n")

  study_summary <-   sem_study_summary_get()
  anno_annotate_position_pivots()

  # ssEnv <- core_init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, tech = tech)
  chart_folder <- io_dir_check_and_create(ssEnv$result_folderChart, "MARKER_PER_PROBE")

  sample_sheet <- utils::read.csv2(io_file_path_build(ssEnv$result_folderData, "1_sample_sheet_original","csv"), sep=";", header = TRUE, stringsAsFactors = FALSE)
  references <- sample_sheet[sample_sheet$Sample_Group == "Reference", "Sample_ID"]

  localKeys <- unique(ssEnv$keys_markers_figures$MARKER)
  tempKeys <- localKeys

  # find the probe with the highet signal values
  high_signal_probes <- 0
  median_signal_probes <- 0
  # low_signal_probes <- 0.1 * high_signal_probes
  low_signal_probes <- 0
  probes_stat_fname <- io_file_path_build( ssEnv$result_folderData, "PROBES_STAT","csv")

  if(((probe_name_max != "cg11680158") | (probe_name_min != "cg11680158")) & file.exists(probes_stat_fname))
    unlink(probes_stat_fname)

  if (file.exists(probes_stat_fname))
    probes_stat <- utils::read.csv2(probes_stat_fname, header = TRUE, stringsAsFactors = FALSE)
  else
  {
    probes_stat <- expand.grid(ssEnv$keys_markers_figures$MARKER, c("MEDIAN", "MAX", "MIN","Q1","Q3","COUNT"))
    probes_stat$POSITION <- ""
    probes_stat$VALUE <- ""
    colnames(probes_stat) <- c("MARKER", "METRIC","POSITION","VALUE")
  }



  localKeys <- unique(localKeys[(localKeys %in% unique(probes_stat[(is.na(probes_stat$VALUE) | probes_stat$VALUE==""),"MARKER"]))])
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = seq_along(localKeys))
  else
    progress_bar <- ""

  if (length(localKeys) != 0)
    for(k in seq_along(localKeys)){
      marker <- as.character(localKeys[k])
      {
        if(ssEnv$showprogress)
          progress_bar(sprintf("Collecting stats for marker %s", marker))
        count_m <- io_get_pivot_both(marker)
        if (nrow(count_m)==0)
          tempKeys <- tempKeys[-k]
        else
        {

          core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Collecting probe stat for marker: ", marker, " \n")
          areas_name <- count_m$AREA
          count_m <- count_m[, -1]
          probes_stat[probes_stat$MARKER == marker  & probes_stat$METRIC== "COUNT","VALUE"] <- nrow(count_m) * ncol(count_m)
          # count per each row the number of no zero values
          markers_sum <- apply(count_m, 1, sum)
          #Q1
          markers_sum_q1 <- quantile(markers_sum, 0.25)
          q1_index <- which.min(abs(markers_sum - markers_sum_q1))
          probes_name_q1 <- areas_name[q1_index]
          if (length(probes_name_q1) > 1)
            probes_name_q1 <- probes_name_q1[1]
          probes_stat[(probes_stat$MARKER == marker) & probes_stat$METRIC== "Q1","POSITION"] <- probes_name_q1
          probes_stat[(probes_stat$MARKER == marker) & probes_stat$METRIC== "Q1","VALUE"] <- markers_sum_q1
          #Q3
          markers_sum_q3 <- quantile(markers_sum, 0.75)
          q3_index <- which.min(abs(markers_sum - markers_sum_q3))
          probes_name_q3 <- areas_name[q3_index]
          if (length(probes_name_q3) > 1)
            probes_name_q3 <- probes_name_q3[1]
          probes_stat[(probes_stat$MARKER == marker) & probes_stat$METRIC== "Q3","POSITION"] <- probes_name_q3
          probes_stat[(probes_stat$MARKER == marker) & probes_stat$METRIC== "Q3","VALUE"] <- markers_sum_q3
          #MAX
          if((probe_name_max != "cg11680158"))
          {
            max_index <- which.max(markers_sum)
            markers_sum_max <- markers_sum[max_index]
            probe_name_max <- areas_name[max_index]
          }
          else
          {
            markers_sum_max <- max(markers_sum)
            if (length(probe_name_max) > 1)
              probe_name_max <- probe_name_max[1]
          }
          probes_stat[(probes_stat$MARKER == marker) & probes_stat$METRIC== "MAX","POSITION"] <- probe_name_max
          probes_stat[(probes_stat$MARKER == marker) & probes_stat$METRIC== "MAX","VALUE"] <- markers_sum_max
          #MIN
          if((probe_name_min != "cg11680158"))
          {
            min_index <- which.min(markers_sum)
            markers_sum_min <- markers_sum[min_index]
            probe_name_min <- areas_name[min_index]
          }
          else
          {
            markers_sum_min <- min(markers_sum)
            if (length(probe_name_min) > 1)
              probe_name_min <- probe_name_min[1]
          }
          probes_stat[(probes_stat$MARKER == marker) & probes_stat$METRIC== "MIN","POSITION"] <- probe_name_min
          probes_stat[(probes_stat$MARKER == marker) & probes_stat$METRIC== "MIN","VALUE"] <- markers_sum_min
          #MEDIAN
          markers_median <- quantile(markers_sum, 0.5)
          median_index <- which.min(abs(markers_sum - markers_median))
          probe_name_median <- areas_name[median_index]
          if (length(probe_name_median) > 1)
            probe_name_median <- probe_name_median[1]
          probes_stat[(probes_stat$MARKER == marker) & probes_stat$METRIC== "MEDIAN","POSITION"] <- probe_name_median
          probes_stat[(probes_stat$MARKER == marker) & probes_stat$METRIC== "MEDIAN","VALUE"] <- markers_median
        }
      }
    }

  utils::write.csv2(probes_stat, probes_stat_fname, row.names = FALSE)
  # remove SIGNAL from tempKeys
  tempKeys <- tempKeys[!grepl("SIGNAL", tempKeys)]
  # AI-027: read via unified dispatcher; cached parquet is the normal
  # path (CASE 1). The fname variable is kept for the error message below.
  fname <- io_pivot_file_name_parquet("SIGNAL", io_signal_figure(), "PROBE","WHOLE")
  signal_pivot_lazy <- io_read_pivot("SIGNAL", io_signal_figure(), "PROBE", "WHOLE")
  if (is.null(signal_pivot_lazy))
    stop("SIGNAL PROBE pivot not found at ", fname,
         " - manhattan plot requires the SEM signal pivot.")
  signal_data <- as.data.frame(signal_pivot_lazy$collect())

  # threshold_data <- fst::read_fst(io_file_path_build( ssEnv$result_folderData,"1_signal_thresholds","fst"))
  threshold_data <- as.data.frame(polars::pl$read_parquet(io_file_path_build( ssEnv$result_folderData,"1_signal_thresholds","parquet")))

  # AI-106 (2026-06-09): no more colname sanitisation. Sample_ID are
  # passed through from the upstream pivot writers (io_signal_save,
  # sem_analyze_population_bulk). The `colnames(signal_data) %in%
  # colnames(marker_data)` match below relies on both pivots writing
  # Sample_ID in the SAME raw format.
  if(max_sample!=0)
    signal_data <- signal_data[,min_sample:max_sample]
  samples <- as.vector(colnames(signal_data[,!grepl("AREA", colnames(signal_data))]))

  signal_data$AREA <- toupper(signal_data$AREA)
  probes_stat <- utils::read.csv2(probes_stat_fname, header = TRUE, stringsAsFactors = FALSE)
  probes_stat <- probes_stat[!is.na(probes_stat$POSITION),]

  ########################################################################################################################################################################
  ############################ PLOT MARKERS FOR A GIVEN PROBE ################################################################################################
  ########################################################################################################################################################################

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = seq_len(length(tempKeys)*5))
  else
    progress_bar <- ""

  for(j in seq_along(tempKeys))
  {
    marker <- as.character(tempKeys[j])
    # check if the plot already has been saved


    # List files that match both strings
    files <- list.files(chart_folder, pattern = paste0(marker,".png"), full.names = TRUE)
    # Check if any matching file exists
    if (length(files) > 0)
      next

    marker_data <- io_get_pivot_both(marker)
    if(nrow(marker_data)==0)
      next
    if(any(!(colnames(signal_data) %in% colnames(marker_data))))
    {
      lost_columns <- colnames(signal_data)[!(colnames(signal_data) %in% colnames(marker_data))]
      #  add column with zeros to marker data
      for(l in seq_along(lost_columns))
        marker_data[[lost_columns[l]]] <- 0
    }
    # marker_data <- marker_data[,colnames(marker_data) %in% colnames(signal_data)]
    marker_data <- marker_data[,colnames(signal_data)]
    marker_data$AREA <- toupper(marker_data$AREA)

    for ( selection in c("min","max","median","q1","q3"))
    {
      # selection <- "max"
      core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), "Generating plot for ", marker, " and ", selection, " probe")
      probe_name <- probes_stat[probes_stat$METRIC==toupper(selection) & max(probes_stat[probes_stat$METRIC==toupper(selection),"VALUE"])==probes_stat[probes_stat$METRIC==toupper(selection),"VALUE"],"POSITION"]
      if(length(probe_name)>1)
        probe_name <- probe_name[1]
      probe_name <- toupper(probe_name)

      probe_signal <- as.vector(t(signal_data[signal_data$AREA==probe_name,!grepl("AREA", colnames(signal_data))]))
      probe_signal[is.na(probe_signal)] <- 0

      threshold_data$PROBE <- toupper(threshold_data$PROBE)
      probe_threshold <- threshold_data[threshold_data$PROBE==probe_name,]

      q1 <- probe_threshold$q1
      q3 <- probe_threshold$q3
      y_med <- probe_threshold$signal_median_values
      iqr <- probe_threshold$iqr
      y_sup <- probe_threshold$signal_superior_thresholds
      y_inf <- probe_threshold$signal_inferior_thresholds


      probe_signal <- as.numeric(probe_signal)
      outlier <-  as.vector(ifelse(probe_signal > y_sup, "Hyper", ifelse(probe_signal < y_inf, "Hypo" ,"Non-outlier")))
      reference_col_id <- which(colnames(signal_data[,-1]) %in% references)
      outlier[reference_col_id] <- "Reference"
      marker_segment_color <-  ifelse(probe_signal > y_sup, hyper_color, ifelse(probe_signal < y_inf, hypo_color ,"white"))

      probe_marker <- marker_data[marker_data$AREA==probe_name,!grepl("AREA", colnames(marker_data))]
      probe_marker[is.na(probe_marker)] <- 0
      probe_marker <- as.numeric(probe_marker)

      outlier[probe_marker==0] <- "Non-outlier"
      outlier[reference_col_id] <- "Reference"

      data_to_plot <- data.frame(samples = samples, signal_to_plot = as.numeric(probe_marker), outlier)
      marker_bar_end <- 0

      # pp <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = as.numeric(as.factor(samples)), y = signal_to_plot, colour = outlier)) +
      #   ggplot2::geom_segment(ggplot2::aes(xend = as.numeric(factor(samples)), yend = marker_bar_end), color = marker_segment_color, alpha = 0.5) +
      #   ggplot2::geom_point(alpha = 0.6) +
      #   ggplot2::scale_color_manual(values = c("Hyper" = hyper_color, "Hypo"= hypo_color, "Non-outlier" = non_outlier_color,"Reference" = reference_samples_color)) +
      #   # ggplot2::labs(x = "Sample", y = marker, title = paste0(marker, " found among all samples for probe ", probe_name)) +
      #   ggplot2::labs(x = "Sample", y = marker, title = "") +
      #   ggplot2::theme_classic() +
      #   ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = ggplot2::element_text(hjust = 0.5))
      # # ggplot2::scale_x_continuous(limits = c(1, NA)) # Adjust the limits to exclude 0

      pp <- ggplot2::ggplot(
        data_to_plot,
        ggplot2::aes(x = as.numeric(factor(samples)), y = signal_to_plot, colour = outlier)
      ) +
        ggplot2::geom_segment(
          ggplot2::aes(
            x = as.numeric(factor(samples)),
            xend = as.numeric(factor(samples)),
            y = marker_bar_end,
            yend = signal_to_plot
          ),
          color = marker_segment_color,
          alpha = 0.5
        ) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::scale_color_manual(values = c(
          "Hyper" = hyper_color,
          "Hypo" = hypo_color,
          "Non-outlier" = non_outlier_color,
          "Reference" = reference_samples_color
        )) +
        ggplot2::labs(x = "Sample", y = marker, title = "") +
        ggplot2::theme_classic() +
        ggplot2::theme(
          legend.position = "none",
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.title = ggplot2::element_text(hjust = 0.5)
        )

      plot_filename <- io_file_path_build(chart_folder,c(selection, probe_name,marker),ssEnv$plot_format)
      ggplot2::ggsave(filename =plot_filename,plot = pp,units = "in", width = 6, height = 2, dpi=as.numeric(ssEnv$plot_resolution_ppi))


      ########################################################################################################################################################################
      ############################ PLOT OUTLIERS FOR A GIVEN PROBE ################################################################################################
      ########################################################################################################################################################################

      data_to_plot <- data.frame(samples = samples, signal_to_plot = as.numeric(probe_signal), outlier)
      delta_label_value <-  ifelse(probe_signal > y_sup, paste("\u03B4 ",
        format(probe_signal - y_sup, scientific=TRUE , digits=2),sep=""),
        ifelse(probe_signal < y_inf,  paste("\u03B4:", format(y_inf - probe_signal,scientific=TRUE, digits=2),sep=""),""))
      delta_label_value[as.numeric(probe_marker)==0] <- ""

      y_shift <- 0
      delta_label_position <-  ifelse(probe_signal > y_sup, (probe_signal + y_shift), ifelse(probe_signal < y_inf, (probe_signal - y_shift) ,0))
      marker_bar_end <-  ifelse(probe_signal > y_sup, y_sup, ifelse(probe_signal < y_inf, y_inf ,probe_signal))
      # marker_segment_color <-  ifelse(probe_signal > y_sup, hyper_color, ifelse(probe_signal < y_inf, hypo_color ,"white"))
      limit_line_color_q <- "red"

      pp <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = as.numeric(as.factor(samples)), y = signal_to_plot, color = outlier)) +
        ggplot2::geom_segment(ggplot2::aes(xend = as.numeric(factor(samples)), yend = marker_bar_end), color = marker_segment_color, alpha = 0.1) +
        ggplot2::geom_hline(yintercept = y_med, linetype = "dashed", color = limit_line_color_median, size = 0.1) +
        ggplot2::geom_hline(yintercept = y_sup, linetype = "dashed", color = limit_line_color, size = 0.1) +
        ggplot2::geom_hline(yintercept = y_inf, linetype = "dashed", color = limit_line_color, size = 0.1) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::geom_hline(yintercept = q1, linetype = "dashed", color = limit_line_color_q, size = 0.1) +
        ggplot2::geom_hline(yintercept = q3, linetype = "dashed", color = limit_line_color_q, size = 0.1) +
        ggplot2::geom_text(ggplot2::aes(x = 0, y = q1, label = ifelse(show_labels, "Q1","")), color = limit_line_color_q, hjust = -0.2, vjust = -0.3) +
        ggplot2::geom_text(ggplot2::aes(x = 0, y = q3, label =  ifelse(show_labels, "Q3","")), color = limit_line_color_q, hjust = -0.2, vjust = -0.3) +
        ggplot2::geom_text(ggplot2::aes(x = 0, y = y_med, label = ifelse(show_labels, "Median","")), color = limit_line_color_median, hjust = -0.2, vjust = -0.3) +
        ggplot2::geom_text(ggplot2::aes(x = 0, y = y_sup, label = ifelse(show_labels, "Upper Limit","")), color = limit_label_color, hjust = -0.2, vjust = 1.5) +
        ggplot2::geom_text(ggplot2::aes(x = 0, y = y_inf, label = ifelse(show_labels, "Lower Limit","")), color = limit_label_color, hjust = -0.2, vjust = -0.3) +
        ggplot2::geom_text(
          ggplot2::aes( angle = 90, x = as.numeric(factor(samples)), y = delta_label_position, label = ifelse(delta_label_value,"")),
          color = marker_segment_color,
          check_overlap = TRUE,
          size = label_font_size,
          hjust = +1,
          vjust = -0.45) +
        ggplot2::scale_color_manual(values = c("Hyper" = hyper_color, "Hypo"= hypo_color, "Non-outlier" = non_outlier_color,"Reference" = reference_samples_color)) +
        ggplot2::scale_shape_manual(values = c("Hyper" = 24, "Hypo" = 24, "Non-outlier" = 19, "Reference" = 6)) +  # Change "Reference" to shape 4 (X mark)
        # ggplot2::labs(x = "Sample", y = "Signal", title = paste0("Signal outliers of all samples for probe ", probe_name)) +
        ggplot2::labs(x = "Sample", y = "Signal", title = "") +
        ggplot2::theme_classic() +
        ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = ggplot2::element_text(hjust = 0.5))
      # ggplot2::scale_x_continuous(limits = c(1, NA)) # Adjust the limits to exclude 0

      plot_filename <- io_file_path_build(chart_folder,c(selection, probe_name,"SIGNAL_OUTLIER"),ssEnv$plot_format)
      ggplot2::ggsave(filename = plot_filename,plot = pp,units = "in", width = 6, height = 4, dpi=as.numeric(ssEnv$plot_resolution_ppi))

      if(ssEnv$showprogress)
        progress_bar(sprintf("Plotting images for probe %s", probe_name))

    }
  }

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Generation of plot completed!")

}


