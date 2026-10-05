#' Compare association results across subsamples of one study
#'
#' Same comparison as \code{\link{meta_association_overlaps_studies}} but
#' within a single study, across the subsamples produced by a replication run.
#' It joins the per-subsample inference results on the area taxonomy and writes,
#' per area, which subsamples called it significant.
#'
#' @param inference_details Inference specification rows identifying the
#'   requests whose results are compared.
#' @param alpha Significance threshold applied to \code{pvalue_column}.
#' @param adjust_per_area Adjust p-values within each area separately.
#' @param adjust_globally Adjust p-values across the whole result set.
#' @param pvalue_column Name of the p-value column to test against
#'   \code{alpha}. Defaults to the all-scope BH-adjusted column.
#' @param statistic_parameter Name of the effect-size column carried into the
#'   comparison table alongside the p-value.
#' @param adjustment_method Multiple-testing correction passed to
#'   \code{stats::p.adjust()}.
#' @param old_label,new_label Optional relabelling of the subsample names in
#'   the output table.
#' @param run_prefix Prefix prepended to the output file names, to keep the
#'   results of several runs side by side.
#' @param result_folder Folder holding the subsamples' results and receiving
#'   the comparison output.
#' @param ... Passed through to the session setup.
#'
#' @return Called for its side effect: CSV files written under the inference
#'   folder. Returns \code{NULL} invisibly, and early if no results are found.
# compare inference associations of different sub samples
#' @examples
#' # The overlap is computed across the subsamples of one study, so the details
#' # are several rows of the same shape, one per subsample run:
#' inference_details <- data.frame(
#'   independent_variable = "Sample_Group",
#'   family_test          = "wilcoxon",
#'   transformation_y     = "none",
#'   aggregation          = "MEAN",
#'   scope                = "INSTANCE"
#' )
#' \dontrun{
#' meta_association_overlaps_subsamples(
#'   inference_details   = inference_details,
#'   statistic_parameter = "PVALUE",
#'   result_folder       = tempdir()
#' )
#' }
#' @export
meta_association_overlaps_subsamples <- function(inference_details,alpha = 0.05, adjust_per_area = FALSE,
  adjust_globally = FALSE,pvalue_column="PVALUE_ADJ_ALL_BH",statistic_parameter, adjustment_method = "BH",
  old_label = NULL, new_label = NULL, run_prefix = "",
  result_folder, ...)
{
  pvalue_column <- core_name_cleaning(pvalue_column)
  if (nrow(inference_details) ==0)
  {
    core_log_event("ERROR: ",format(Sys.time(), "%a %b %d %X %Y")," No inference details found!")
    return(NULL)
  }
  ssEnv <- core_init_env( result_folder =  result_folder, start_fresh = FALSE, ...)

  count_family <- length(unique(inference_details$family_test))
  if(count_family>1)
  {
    core_log_event("ERROR: ",format(Sys.time(), "%a %b %d %X %Y")," Inference details have different family tests !")
    return(NULL)
  }

  family_test <- as.character(inference_details[1,"family_test"])
  independent_variable <- as.character(inference_details[1,"independent_variable"])

  core_log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), " SEMseeker will perform the inference detail association analysys for sub samples overlaps.")

  color_palette <- ssEnv$color_palette
  localKeys <- unique(ssEnv$keys_areas_subareas_markers_figures[,c("MARKER","AREA")])

  aggregated_results <- data.frame()
  for (a in 1:nrow(localKeys))
  {
    MARKER <- localKeys[a, "MARKER"]
    AREA <- localKeys[a, "AREA"]
    # for each SAMPLES_SQL_CONDITION in inference_details
    for (i in seq_len(nrow(inference_details)))
    {
      # get the inference detail
      temp_res <- assoc_results_get(inference_detail = inference_details[i,], marker = MARKER,
        area = AREA, scope = "INSTANCE",
        pvalue_column= pvalue_column,
        adjustment_method = adjustment_method, significance = NULL)
      if(nrow(temp_res) != 0)
      {
        temp_res$SAMPLES_SQL_CONDITION <- core_name_cleaning(inference_details[i,"samples_sql_condition"])
        aggregated_results <- plyr::rbind.fill(aggregated_results, temp_res)
      }
    }
  }

  # replace empty SAMPLES_SQL_CONDITION with "ALL"
  aggregated_results$SAMPLES_SQL_CONDITION[aggregated_results$SAMPLES_SQL_CONDITION == ""] <- "ALL"

  aggregated_results <- assoc_filter_sql(inference_details$association_results_sql_condition, aggregated_results)

  core_log_event("INFO: ",format(Sys.time(), "%a %b %d %X %Y")," inference files aggregated!")

  # create a file with max of pvalue column for the inference detail and mean of statistic parameter
  if(nrow(aggregated_results)>0)
    # write a file for each marker
    for (marker in unique(aggregated_results$MARKER))
    {
      tt <- subset(aggregated_results, MARKER == marker)
      # AI-255: per-instance artefacts, said by the taxonomy
      tt <- subset(tt, SCOPE == "INSTANCE")

      tt$KEY <- paste0(tt$AREA,"_",tt$SUBAREA,"_",tt$MARKER,"_",tt$FIGURE,"_",tt$AREA_OF_TEST)
      SPLIT <- split(tt$KEY, tt$SAMPLES_SQL_CONDITION)
      # get the common keys
      common_keys <- Reduce(intersect, SPLIT)
      # get the common keys
      tt <- subset(tt, KEY %in% common_keys)
      # remove KEY column
      tt$KEY <- NULL

      if(statistic_parameter!="")
      {
        tt <- tt[,c("MARKER","FIGURE","SCOPE","AREA","SUBAREA","AGGREGATION","AREA_OF_TEST",statistic_parameter, pvalue_column)]
        # get only the taxonomy key + AREA_OF_TEST common to SAMPLES_SQL_CONDITION
        tt <- tt %>%
          dplyr::group_by(.data$AREA, .data$SUBAREA, .data$MARKER, .data$FIGURE, .data$SCOPE, .data$AGGREGATION, .data$AREA_OF_TEST) %>%
          dplyr::summarise(
            alpha = max(get(pvalue_column), na.rm = TRUE),
            statistic_parameter = mean(get(statistic_parameter), na.rm = TRUE)
          ) %>%
          dplyr::ungroup()
        colnames(tt)[colnames(tt) == statistic_parameter] <- statistic_parameter
      }
      else
      {
        tt <- tt[,c("MARKER","FIGURE","SCOPE","AREA","SUBAREA","AGGREGATION","AREA_OF_TEST",pvalue_column)]
        # summarise grouping by "AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST" and calculate the max of the pvalues
        tt <- tt %>%
          dplyr::group_by(.data$AREA, .data$SUBAREA, .data$MARKER, .data$FIGURE, .data$SCOPE, .data$AGGREGATION, .data$AREA_OF_TEST) %>%
          dplyr::summarise(
            alpha = max(get(pvalue_column), na.rm = TRUE)
          ) %>%
          dplyr::ungroup()
      }
      # rename pvalue_column to the original name
      colnames(tt)[colnames(tt) == "alpha"] <- pvalue_column

      dest_folder <- io_dir_check_and_create(ssEnv$result_folderInference,c("ASSOCIATION_CROSS_SAMPLE_ANALYSIS",core_name_cleaning(unique(inference_details$association_results_sql_condition),"ALL"), core_name_cleaning(paste(family_test, pvalue_column, run_prefix))))
      inference_detail <- inference_details[1,]
      inference_detail$samples_sql_condition <- ""
      filename <- io_inference_file_name(inference_detail, marker, dest_folder ,prefix="" )

      utils::write.csv2(tt, filename, row.names = FALSE)
    }


  aggregated_results <- data.frame()
  for (signif in c(TRUE))
  {
    for (a in seq_len(nrow(localKeys)))
    {
      MARKER <- localKeys[a, "MARKER"]
      AREA <- localKeys[a, "AREA"]
      # for each SAMPLES_SQL_CONDITION in inference_details
      for (i in seq_len(nrow(inference_details)))
      {
        temp_res <- assoc_results_get(inference_detail = inference_details[i,], marker = MARKER,
          area = AREA, scope = "INSTANCE",
          pvalue_column= pvalue_column,
          adjustment_method = adjustment_method, significance = signif)
        if(nrow(temp_res) != 0)
          temp_res$SAMPLES_SQL_CONDITION <- core_name_cleaning(inference_details[i,"samples_sql_condition"])
        aggregated_results <- plyr::rbind.fill(aggregated_results, temp_res)
      }
    }

    # replace empty SAMPLES_SQL_CONDITION with "ALL"
    aggregated_results$SAMPLES_SQL_CONDITION[aggregated_results$SAMPLES_SQL_CONDITION == ""] <- "ALL"
    aggregated_results <- assoc_filter_sql(inference_details$association_results_sql_condition, aggregated_results)

    if(nrow(aggregated_results)==0)
      next

    if(statistic_parameter!="")
      aggregated_results <- unique(aggregated_results[,c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST","SAMPLES_SQL_CONDITION",statistic_parameter, pvalue_column)])
    else
      aggregated_results <- unique(aggregated_results[,c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST","SAMPLES_SQL_CONDITION", pvalue_column)])

    # reshape to a table with SAMPLES_SQL_CONDITION as columns, area as rows and values as cell without aggreagtion
    aggregated_results_table <- reshape2::dcast(aggregated_results, AREA + SUBAREA + MARKER + FIGURE + AREA_OF_TEST ~ SAMPLES_SQL_CONDITION, value.var = pvalue_column)
    #
    studies_to_comb <- na.omit(unique(aggregated_results$SAMPLES_SQL_CONDITION))
    # calculate the mean of the pvalues
    aggregated_results_table[, pvalue_column] <- mean(aggregated_results_table[,studies_to_comb], na.rm = TRUE)

    # rename columns
    studies_to_comb_cols <- paste0(studies_to_comb,"_",pvalue_column)
    colnames(aggregated_results_table) <- c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST",studies_to_comb_cols, pvalue_column)


    if(statistic_parameter != "")
    {
      aggregated_results_table_statistic_parameter <- reshape2::dcast(aggregated_results, AREA + SUBAREA + MARKER + FIGURE + AREA_OF_TEST ~ SAMPLES_SQL_CONDITION, value.var = statistic_parameter)
      aggregated_results_table_statistic_parameter[, statistic_parameter] <- mean(aggregated_results_table_statistic_parameter[,studies_to_comb], na.rm = TRUE)

      studies_to_comb_cols <- paste0(studies_to_comb,"_",statistic_parameter)
      colnames(aggregated_results_table_statistic_parameter) <- c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST",studies_to_comb_cols, statistic_parameter)

      aggregated_results_table <- merge(aggregated_results_table, aggregated_results_table_statistic_parameter, by = c("AREA", "SUBAREA", "MARKER", "FIGURE", "AREA_OF_TEST"))
    }

    # AI-106 (2026-06-09): removed legacy gsub round-trip (same fix as
    # inter_study_association_overlaps.R). Post-AI-106 all CSVs preserve
    # raw names from upstream annotation; no transformation needed.
    markers <- unique(aggregated_results_table[, c("MARKER")])
    # for (i in seq_along(markers))
    # {
    #   marker <- markers[i]
    #   filename <- io_inference_file_name(inference_detail, marker,ssEnv$result_folderInference,file_extension = "csv",suffix = "", prefix = "")
    #   aggregated_results_table_marker <- subset(aggregated_results_table, MARKER == marker)
    #   aggregated_results_table_marker$DEPTH <- 3
    #   write.csv2(aggregated_results_table_marker, filename, row.names = F)
    # }
    filename <- io_inference_file_name(inference_detail, paste0(markers, collapse = "_") ,dest_folder,file_extension = "csv",
      suffix = "AGGREGATED", prefix = ifelse(signif,"SIGNIFICANT","NOT_SIGNIFICANT"))
    write.csv2(aggregated_results_table, filename, row.names = FALSE)
    core_log_event("INFO: ",format(Sys.time(), "%a %b %d %X %Y"),"  aggregated files saved!")
    # for( j in 2:length(studies_to_comb))
    {
      # studies_comb <- combinat::combn(studies_to_comb, j)
      # if (j==length(studies_to_comb))
      #   studies_comb <- data.frame("st"=studies_comb)
      # message(studies_comb)
      # for(k in 1: ncol(studies_comb))
      {
        # results_inference_comb <- subset(aggregated_results, SAMPLES_SQL_CONDITION %in% studies_comb[,k])
        results_inference_comb <- aggregated_results
        keys <- unique(results_inference_comb[, c("SUBAREA", "AREA", "MARKER", "FIGURE")])
        # sample_condition_count <- length(na.omit(unique(results_inference_comb$SAMPLES_SQL_CONDITION)))

        # Load library
        # library(VennDiagram)
        for (i in seq_len(nrow(keys)))
        {
          # i <- 1
          area_set <- results_inference_comb[
            results_inference_comb$AREA == keys[i, ]$AREA &
              results_inference_comb$SUBAREA == keys[i, ]$SUBAREA &
              results_inference_comb$MARKER == keys[i, ]$MARKER &
              results_inference_comb$FIGURE == keys[i, ]$FIGURE
            , c("AREA_OF_TEST", "SAMPLES_SQL_CONDITION"), ]

          if(nrow(area_set)==0)
            next
          # replace samples_sql_condition with alternative names
          area_set$SAMPLES_SQL_CONDITION <- as.character(area_set$SAMPLES_SQL_CONDITION)
          area_set$SAMPLES_SQL_CONDITION <- core_name_cleaning(area_set$SAMPLES_SQL_CONDITION)
          inference_details$samples_sql_condition <- as.character(inference_details$samples_sql_condition)
          area_set$SAMPLES_SQL_CONDITION <- as.character(area_set$SAMPLES_SQL_CONDITION)
          for (j in seq_along(new_label))
          {
            area_set[core_name_cleaning(area_set$SAMPLES_SQL_CONDITION)==core_name_cleaning(old_label[j]),"SAMPLES_SQL_CONDITION"] <- core_name_cleaning(new_label[j])
            inference_details[core_name_cleaning(inference_details$samples_sql_condition)==core_name_cleaning(old_label[j]),"samples_sql_condition"] <- core_name_cleaning(new_label[j])
          }
          # table(is.na(results_inference_comb))
          area_set <- na.omit(area_set)
          SPLIT <- split(area_set$AREA_OF_TEST, area_set$SAMPLES_SQL_CONDITION)

          categories <- core_name_cleaning(unique(inference_details$samples_sql_condition))
          # put ALL where categories==""
          categories <- ifelse(categories == "", "ALL", categories)
          if(length(categories)<2)
            next
          # AI-044 (2026-06-09): use shared `util_pretty_label()` helper.
          categories <- util_pretty_label(categories)
          dest_folder <- io_dir_check_and_create(ssEnv$result_folderChart,c("ASSOCIATION_CROSS_SAMPLE_ANALYSIS",core_name_cleaning(unique(inference_details$association_results_sql_condition),"ALL"), core_name_cleaning(paste(family_test, pvalue_column, run_prefix))))
          filename <-
            paste(
              dest_folder,  "/",
              keys[i, ]$AREA,
              "_",
              keys[i, ]$SUBAREA,
              "_",
              keys[i, ]$MARKER,
              "_",
              keys[i, ]$FIGURE,
              "_",
              independent_variable,
              "_",
              core_name_cleaning(pvalue_column),
              "_",
              ifelse(signif,"SIGNIFICANT","NOT_SIGNIFICANT"),
              "_",
              alpha,
              "_",
              core_name_cleaning(family_test),
              "_UPSET_PLOT.",
              ssEnv$plot_format,
              sep = ""
            )
          # Set up the Venn diagram parameters
          # color_palette <- color_palette[length(SPLIT)]

          # # Create a Venn diagram using ggvenn
          # venn_plot <- ggvenn::ggvenn(
          #   data = SPLIT,
          #   fill_color = color_palette[seq_along(SPLIT)],
          #   set_name_color = "white"
          # )
          #
          # # Customize the plot with ggplot2
          # venn_plot <- venn_plot +
          #   ggplot2::theme_void() +
          #   ggplot2::theme(
          #     legend.position = "bottom",
          #     legend.title = ggplot2::element_blank(),
          #     legend.text = ggplot2::element_text(size = 10, face = "bold")
          #   ) +
          #   ggplot2::scale_fill_identity() +
          #   ggplot2::scale_fill_manual(
          #     values = color_palette[seq_along(SPLIT)],
          #     labels = categories
          #   ) +
          #   ggplot2::guides(fill = ggplot2::guide_legend(title = NULL, override.aes = list(color = NA)))
          # ggplot2::ggsave(filename, venn_plot, width = 2048, height = 2048, dpi = as.numeric(ssEnv$plot_resolution_ppi), unit="px")



          area_set <- UpSetR::fromList(SPLIT)
          if (is.null(dim(area_set)))
            next
          # remove from column name _PVALUE_ADJ
          to_remove <- paste0("_",pvalue_column)
          # AI-044 (2026-06-09): use shared `util_pretty_label()` helper.
          categories <- util_pretty_label(categories)
          colnames(area_set) <- gsub(to_remove, "", colnames(area_set))
          colnames(area_set) <- gsub("_", " ", colnames(area_set))
          missed_columns <- setdiff(categories, colnames(area_set))
          for (missed_column in missed_columns)
            area_set[[missed_column]] <- 0

          # set all na to 0
          area_set[is.na(area_set)] <- 0

          # set all values not 0 to 1
          area_set[,categories] <- as.numeric(area_set[,categories] > 0)


          # sum by row
          degree <- rowSums(area_set[,categories])
          intersect_length <- nrow(unique(area_set))
          intersect_length <- min(intersect_length, 30)
          # Save the plot to a file (e.g., PNG)

          if (any(degree == length(categories)))
            # upset(area_set, sets = sets, decreasing = c(FALSE,TRUE), nintersects = 30)
            my_plot <- UpSetR::upset(data = area_set, sets = categories, decreasing = c(FALSE,TRUE),
              nintersects = 30,
              order.by = c("freq", "degree"),
              main.bar.color = ssEnv$color_palette[1], matrix.color = ssEnv$color_palette[2], sets.bar.color = ssEnv$color_palette[3],
              queries = list(list(query = UpSetR::intersects, params = list(categories), color = ssEnv$color_palette[4])))
          else
            my_plot <- UpSetR::upset(data = area_set, sets = categories, decreasing = c(FALSE,TRUE), nintersects = 30,
              main.bar.color = ssEnv$color_palette[1], matrix.color = ssEnv$color_palette[2], sets.bar.color = ssEnv$color_palette[3])


          # First, open a file device for plotting (PNG, PDF, etc.)
          grDevices::png(filename, width = 9, height = 9, units="in", res = as.numeric(ssEnv$plot_resolution_ppi))
          print(my_plot)
          # Close the file device to save the plot
          grDevices::dev.off()

          # ggplot2::ggsave(filename, upset_plot, width = 2048, height = 2048, dpi = as.numeric(ssEnv$plot_resolution_ppi), unit="px")

          # Chart
          # Set up the Venn diagram parameters
          # Plot Venn diagram
          # suppressMessages(
          #   VennDiagram::venn.diagram(
          #     x = SPLIT,
          #     # fill = color_palette[seq_along(SPLIT)],
          #     alpha = 0.5,
          #     category.names = categories,
          #     cat.col = rep("black", length(SPLIT)),
          #     cat.cex = 1.2,
          #     cat.fontface = "bold",
          #     cex = 1.5,
          #     output = TRUE,
          #     individuals.in.intersections = TRUE,
          #     disable.logging = TRUE,
          #     filename = filename,
          #     resolution = 600,
          #     units = "px",
          #     height = 2048,
          #     width = 2048
          #   )
          # )

          core_log_event("DEBUG: ",format(Sys.time(), "%a %b %d %X %Y"),"  venn diagram completed !")

          overlaps <- Reduce(intersect, SPLIT)
          if(length(overlaps)>0)
          {
            filename <- io_inference_file_name(inference_detail, paste(keys[i, ]$AREA, keys[i, ]$SUBAREA, keys[i, ]$MARKER, keys[i, ]$FIGURE , sep="_"),dest_folder,file_extension = "csv",suffix = "OVERLAPS", prefix = ifelse(signif,"SIGNIFICANT","NOT_SIGNIFICANT"))
            overlaps <- data.frame("AREA_OF_TEST" = overlaps)
            utils::write.csv2(overlaps, filename, append = TRUE)
          }
        }
      }
    }

    core_log_event("INFO: ",format(Sys.time(), "%a %b %d %X %Y"),"  job completed !")
    # unlink(dest_folder, recursive = TRUE)
    # remove all files ending with .log
    # files <- list.files(ssEnv$result_folderInference, pattern = ".log")
    # for (f in seq_along(files))
    # {
    #   file <- files[f]
    #   file.remove(paste(ssEnv$result_folderInference, "/", file, sep = ""))
    # }
  }
  core_close_env()
}
