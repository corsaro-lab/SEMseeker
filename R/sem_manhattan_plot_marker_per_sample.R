#' Manhattan plot of SEM markers for a single sample
#'
#' @param sample_name character. Sample identifier as it appears in the
#'   sample sheet (default \code{"NAME"}).
#' @param probes_range integer vector. Probe index range to display
#'   (default \code{1000:2000}).
#' @param hyper_color character. Colour for hypermethylated probes
#'   (default \code{"blue"}).
#' @param hypo_color character. Colour for hypomethylated probes
#'   (default \code{"orange"}).
#' @param non_outlier_color character. Colour for non-outlier probes
#'   (default \code{"grey"}).
#' @param limit_label_color character. Colour for threshold labels.
#' @param result_folder character. Path to the SEMseeker result folder.
#' @param maxResources numeric. Maximum percentage of CPU cores to use
#'   (default 90).
#' @param parallel_strategy character. Parallelisation backend
#'   (default \code{"multisession"}).
#'   Asking for \code{"multicore"} is accepted and converted to
#'   \code{"multisession"}: it means fork(), which is unsafe with this
#'   package's native thread pool on every platform that offers it, and
#'   absent on Windows. The conversion is logged.
#' @param ... Additional arguments passed to \code{core_init_env()}.
#'
#' @return Invisibly \code{NULL}. A PNG plot is saved under
#'   \code{Charts/MARKER_PER_SAMPLE/} in \code{result_folder}.
#' @examples
#' result_dir <- tempdir()
#' \dontrun{
#' sem_manhattan_plot_marker_per_sample(
#'   sample_name  = "CASE_001",
#'   probes_range = 1:5000,
#'   result_folder = result_dir
#' )
#' }
#' @export
#'
sem_manhattan_plot_marker_per_sample <- function( sample_name = "NAME", probes_range = 1000:2000,
  hyper_color = "blue",   hypo_color = "orange",   non_outlier_color = "grey", limit_label_color = ssEnv$color_palette[1],
  result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, ...)
{


  limit_line_color <- ssEnv$color_palette[2]
  limit_line_color_q <- "cyan"
  limit_line_color_median <- "black"

  ssEnv <- core_init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE,
    figures=c("BOTH","HYPER","HYPO"), areas="PROBE", ...)
  # ssEnv <- core_init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE,
  #   figures=c("BOTH","HYPER","HYPO"), areas="PROBE", markers = "DELTAQ")
  chart_folder <- io_dir_check_and_create(ssEnv$result_folderChart, "MARKER_PER_SAMPLE")

  localKeys <- unique(ssEnv$keys_markers_figures$MARKER)
  tempKeys <- localKeys
  for(k in seq_along(localKeys)){
    marker <- as.character(localKeys[k])
    pivot_subfolder <- io_dir_check_and_create(result_folderPivot,marker)
    fname_both <- io_file_path_build( pivot_subfolder ,c(marker, "BOTH", "PROBE","WHOLE"),"csv", add_gz = TRUE)
    fname_hypo <- io_file_path_build( pivot_subfolder ,c(marker, "HYPO", "PROBE","WHOLE"),"csv", add_gz = TRUE)
    fname_hyper <- io_file_path_build( pivot_subfolder ,c(marker, "HYPER", "PROBE","WHOLE"),"csv", add_gz = TRUE)

    if (!file.exists(fname_both) | !file.exists(fname_hypo) | !file.exists(fname_hyper))
      tempKeys <- tempKeys[-k]
  }

  # remove SIGNAL from tempKeys
  tempKeys <- tempKeys[!grepl("SIGNAL", tempKeys)]

  for(j in seq_along(tempKeys))
  {

    marker <- as.character(tempKeys[j])
    pivot_subfolder <- io_dir_check_and_create(result_folderPivot,marker)
    fname_both <- io_file_path_build( pivot_subfolder ,c(marker, "BOTH", "PROBE","WHOLE"),"csv", add_gz = TRUE)
    if(!file.exists(fname_both))
      next

    data_name <- utils::read.csv2(fname_both, sep  =  ";", nrows = 1)
    marker_data_both <- utils::read.csv2(fname_both, sep  =  ";", skip = 1)
    # remove column X1
    marker_data_both <- marker_data_both[, !grepl("X1", colnames(marker_data_both))]
    # rename column SAMPLE_GROUP as SAMPLEID
    colnames(marker_data_both)[colnames(marker_data_both)=="SAMPLE_GROUP"] <- "SAMPLEID"
    colnames(marker_data_both) <- colnames(data_name)
    # marker_data_both <- utils::read.csv2(fname_both, sep  =  ";")
    # marker_data_both <- marker_data_both[-1,]

    # if sample_name is not in the columns, then skip
    if(!(sample_name %in% colnames(marker_data_both)))
    {
      # find the sample with highest number of probes
      sample_name <- colnames(marker_data_both)[which.max(colSums(marker_data_both[,2:ncol(marker_data_both)]>0))]
    }
    if(!(sample_name %in% colnames(marker_data_both)))
      next
    sample_marker_both <- marker_data_both[,sample_name]
    sample_marker_both <- as.numeric(as.vector(t(sample_marker_both)))
    probes <- marker_data_both$SAMPLEID

    fname_hypo <- io_file_path_build( pivot_subfolder ,c(marker, "HYPO", "PROBE","WHOLE"),"csv", add_gz = TRUE)
    if(!file.exists(fname_hypo))
      next


    marker_data_hypo <- utils::read.csv2(fname_hypo, sep  =  ";")
    marker_data_hypo <- marker_data_hypo[-1,]
    # add to marker data hypo missed rows from both
    missed <- setdiff(marker_data_both$SAMPLEID,marker_data_hypo$SAMPLEID)
    missed <- data.frame(SAMPLEID = missed, matrix(0, nrow = length(missed), ncol = ncol(marker_data_hypo)-1))
    colnames(missed) <- colnames(marker_data_hypo)
    marker_data_hypo <- rbind(marker_data_hypo,missed)
    # add to marker missed column from both
    missed <- setdiff(colnames(marker_data_both),colnames(marker_data_hypo))
    missed <- data.frame(SAMPLEID = probes, matrix(0, nrow = length(probes), ncol = length(missed)))
    marker_data_hypo <- merge(marker_data_hypo,missed, by ="SAMPLEID")

    # give the same order as in both
    marker_data_hypo <- marker_data_hypo[match(marker_data_both$SAMPLEID,marker_data_hypo$SAMPLEID),]

    sample_marker_hypo <- marker_data_hypo[,sample_name]
    sample_marker_hypo <- as.numeric(as.vector(t(sample_marker_hypo)))
    sample_marker_hypo <- ifelse(sample_marker_hypo!=0,"Hypo","")

    fname_hyper <- io_file_path_build( pivot_subfolder ,c(marker, "HYPER", "PROBE","WHOLE"),"csv", add_gz = TRUE)
    if(!file.exists(fname_hyper))
      next
    marker_data_hyper <- utils::read.csv2(fname_hyper, sep  =  ";")
    marker_data_hyper <- marker_data_hyper[-1,]
    # add to marker data hypo missed rows from both
    missed <- setdiff(marker_data_both$SAMPLEID,marker_data_hyper$SAMPLEID)
    missed <- data.frame(SAMPLEID = missed, matrix(0, nrow = length(missed), ncol = ncol(marker_data_hyper)-1))
    colnames(missed) <- colnames(marker_data_hyper)
    marker_data_hyper <- rbind(marker_data_hyper,missed)

    # add to marker missed column from both
    missed_samples <- setdiff(colnames(marker_data_both),colnames(marker_data_hyper))
    missed <- data.frame(SAMPLEID = probes, matrix(0, nrow = length(probes), ncol = length(missed_samples)))
    colnames(missed) <- c("SAMPLEID",missed_samples)
    marker_data_hyper <- merge(marker_data_hyper,missed, by ="SAMPLEID")

    # give the same order as in both
    marker_data_hyper <- marker_data_hyper[match(marker_data_both$SAMPLEID,marker_data_hyper$SAMPLEID),]

    sample_marker_hyper <- marker_data_hyper[,sample_name]
    sample_marker_hyper <- as.numeric(as.vector(t(sample_marker_hyper)))
    sample_marker_hyper <- ifelse(sample_marker_hyper!=0,"Hyper","")

    outlier <- paste(sample_marker_hyper,sample_marker_hypo, sep ="")
    outlier <- ifelse(outlier=="","Non-outlier",outlier)

    marker_segment_color <-  ifelse(outlier =="Hyper", hyper_color, ifelse(outlier=="Hypo", hypo_color ,"white"))

    sample_marker_both[outlier=="Non-outlier"] <- 0
    data_to_plot <- data.frame(probes = probes, signal_to_plot = as.numeric(sample_marker_both), outlier)
    marker_bar_end <- 0

    data_to_plot <- data_to_plot[probes_range,]
    marker_segment_color <- marker_segment_color[probes_range]
    outlier <- outlier[probes_range]
    probes <- probes[probes_range]

    # Assuming data_to_plot is your data frame and it's already filtered or adjusted
    # to exclude the factor level corresponding to x=0

    # calculate ratio of axis as
    ratio <- length(probes)/(max(data_to_plot$signal_to_plot) * 4)

    pp <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = as.numeric(as.factor(probes)), y = signal_to_plot, color = outlier)) +
      ggplot2::geom_segment(ggplot2::aes(xend = as.numeric(factor(probes)), yend = marker_bar_end), color = marker_segment_color, alpha = 0.5) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::scale_color_manual(values = c("Hyper" = hyper_color, "Hypo"= hypo_color, non_outlier_color)) +
      ggplot2::labs(x = "Probe", y = marker, title ="") +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::scale_x_continuous(limits = c(1, NA)) +
      ggplot2::coord_fixed(ratio = ratio)

    # pp <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = as.numeric(as.factor(probes)), y = signal_to_plot, color = outlier)) +
    #   ggplot2::geom_point(alpha = 0.6) +
    #   ggplot2::geom_segment(ggplot2::aes(xend = as.numeric(factor(probes)), yend = marker_bar_end), color = marker_segment_color, alpha = 0.5) +
    #   ggplot2::scale_color_manual(values = c("Hyper" = hyper_color, "Hypo"= hypo_color, non_outlier_color)) +
    #   # ggplot2::labs(x = "Probe", y = marker, title = paste0(marker," of all probes for sample ", sample_name)) +
    #   ggplot2::labs(x = "Probe", y = marker, title ="") +
    #   ggplot2::theme_classic() +
    #   ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = ggplot2::element_text(hjust = 0.5))

    plot_filename <- io_file_path_build(chart_folder,c(sample_name,marker),ssEnv$plot_format)
    ggplot2::ggsave(filename =plot_filename,plot = pp,
      # units = "in", width = 6, height = 2,
      dpi=as.numeric(ssEnv$plot_resolution))
  }

}
