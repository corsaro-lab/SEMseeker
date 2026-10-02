enrich_lollipop_plot <- function(data, rules, file_prfx,path, disease,  top=50, column_sorting="SS_RANK")
{
  ssEnv <- core_get_session_info()

  column_of_id <- rules["column_of_id"]
  column_of_adj_pvalue <- rules["column_of_adj_pvalue"]
  column_of_description <- rules["column_of_description"]
  column_of_enrichment <- rules["column_of_enrichment"]

  colnames(data)[which(colnames(data)==column_of_id[1,1])] <- "ID"
  colnames(data)[which(colnames(data)==column_of_adj_pvalue[1,1])] <- "P_Value"
  colnames(data)[which(colnames(data)==column_of_description[1,1])] <- "Description"
  colnames(data)[which(colnames(data)==column_of_enrichment[1,1])] <- "Enrichment"

  if(min(data$P_Value, na.rm=TRUE)==0)
  {
    min_val <- min(data[data$P_Value>0,"P_Value"], na.rm = TRUE)
    data[data$P_Value==0,"P_Value"] <- min_val
  }

  # Transform highest_p to -log10(highest_p)
  data <- data %>%
    dplyr::mutate(log_fdr = -log10(.data$P_Value))

  # sort by P_Value desc
  data <- data[order(data$log_fdr, decreasing = TRUE),]

  figures <- unique(data$FIGURE)
  categories <- unique(data$SS_CATEGORY)
  for (c in categories)
  {
    data_to_plot <- data[data$SS_CATEGORY == c, ]

    # if(any(data_to_plot$by_keyword))
    #   data_to_plot <- subset(data_to_plot, by_keyword == TRUE)

    if(nrow(data_to_plot) == 0)
      next

    data_to_plot$key <- data_to_plot$MARKER

    # divide in 5 bins the enrichment
    # data_to_plot$Fold_Enrichment_Shape <- cut(data_to_plot$Enrichment, breaks = 5, labels = c("21", "22", "23", "24", "25"))

    # divide in 5 quantile
    # more_than_one <- length(unique(data_to_plot$Enrichment))
    # if (more_than_one > 1)
    #   data_to_plot$Fold_Enrichment_Shape <- cut(data_to_plot$Enrichment, quantile(data_to_plot$Enrichment, probs = seq(0, 1, 0.2)), labels = c("21", "22", "23", "24", "25"))
    # else
    #   data_to_plot$Fold_Enrichment_Shape <- "21"


    # # Create the lollipop plot with vertical dodge and shapes based on Enrichment
    # ggplot(data_to_plot, aes(x = log_fdr, y = Description, size = Enrichment, color = key, shape = factor(Fold_Enrichment_Shape))) +
    #   geom_point(aes(size = Enrichment), fill = NA, position = position_dodgev(height = 0.5), stroke = 1.5) +
    #   ggplot2::scale_shape_manual(values = c(21, 22, 23, 24,25)) + # Customize shape values as needed
    #   ggplot2::scale_size_continuous(range = c(3, 10)) + # Adjust the range for point sizes
    #   labs(
    #     title = "",
    #     x = '-log10(P_Value)',
    #     y = 'Description'
    #   ) +
    #   ggplot2::guides(shape = FALSE) + # Hide the shape legend
    #   ggplot2::theme_minimal() +
    #   ggplot2::theme(
    #     axis.title.x = ggplot2::element_text(size = 12),
    #     axis.title.y = ggplot2::element_text(size = 12),
    #     axis.text.x = ggplot2::element_text(size = 10),
    #     axis.text.y = ggplot2::element_text(size = 10),
    #     plot.title = ggplot2::element_text(size = 14, face = 'bold')
    #   )

    # Define your custom colors
    custom_colors <- ssEnv$color_palette


    # get for each Description the mean of log_fdr
    data_to_plot <- data_to_plot %>%
      dplyr::group_by(.data$Description) %>%
      dplyr::mutate(mean_log_fdr = mean(.data$log_fdr)) %>%
      dplyr::ungroup()
    data_to_plot <- as.data.frame(data_to_plot)
    # Ensure the data_to_plot$Description is a factor and reorder it based on log_fdr
    # data_to_plot$Description <- reorder(data_to_plot$Description, data_to_plot$mean_log_fdr, decreasing = TRUE)
    # order data by mean_log_fdr descending
    data_to_plot <- data_to_plot[order(data_to_plot[,column_sorting], decreasing = FALSE),]
    data_to_plot[,performance_category] <- data_to_plot[,performance_category] %>% as.factor()

    unique_description <- unique(data_to_plot$Description)[seq_len(top)]
    data_to_plot <- data_to_plot[data_to_plot$Description %in% unique_description,]

    data_to_plot <- as.data.frame(data_to_plot)
    # count the number of unique values in the performance_category column
    n_unique_values <- length(unique(data_to_plot[,performance_category]))
    shape_values <- 21:(22+n_unique_values)

    data_to_plot$Description <- ifelse(data_to_plot$by_keyword, paste(data_to_plot$Description, " *"), as.character(data_to_plot$Description))

    # Ensure Description is a factor and preserve its order
    data_to_plot$Description <- factor(data_to_plot$Description, levels = unique(data_to_plot$Description))

    # invert the order of the Description
    data_to_plot$Description <- factor(data_to_plot$Description, levels = rev(levels(data_to_plot$Description)))

    # Create the lollipop plot with vertical dodge and shapes based on key
    ggplot2::ggplot(data_to_plot, ggplot2::aes(x = .data$log_fdr, y = .data$Description, size = .data$Enrichment, color = eval(parse(text=performance_category)))) +
      ggplot2::geom_point(ggplot2::aes(shape = eval(parse(text=performance_category))),
        fill = NA,
        # Was ggstance::position_dodgev(height = 0.5), from a package that was
        # never declared as a dependency: the call could only ever fail with
        # "could not find function". ggplot2 dodges along whichever axis is
        # discrete, which here is y, so position_dodge2() is the vertical dodge
        # this plot asked for, without adding an archived dependency.
        position = ggplot2::position_dodge2(width = 0.5),
        stroke = 1.5) +
      ggplot2::scale_shape_manual(values = shape_values) + # Customize shape values as needed
      ggplot2::scale_size_continuous(range = c(3, 10)) + # Adjust the range for point sizes
      ggplot2::scale_color_manual(values = custom_colors) + # Add custom colors for points
      ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = ssEnv$color_palette[1], size = 1) + # Add vertical line at FDR = 0.05
      ggplot2::geom_vline(xintercept = -log10(0.01), linetype = "dashed", color = ssEnv$color_palette[2], size = 1) + # Add vertical line at FDR = 0.05
      labs(
        title = "",
        x = '-log10(P_Value)',
        y = 'Description',
        shape = performance_category,  # Set the shape legend title dynamically
        color = performance_category   # Set the color legend title dynamically
      ) +
      ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 5))) + # Show the shape legend
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 12),
        axis.title.y = ggplot2::element_text(size = 12),
        axis.text.x = ggplot2::element_text(size = 10),
        axis.text.y = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_text(size = 14, face = 'bold')
      ) +
      xlim(c( min(data_to_plot$log_fdr) - 1 , max(data_to_plot$log_fdr) + 0.5)) # Extend the x-axis range

    # ggplot(data_to_plot, aes(x = log_fdr, y = Description, size = Enrichment, color = key)) +
    #   geom_point(aes(shape = factor(cut(Enrichment, breaks = 5, labels = FALSE))),
    #     fill = NA,
    #     position = position_dodgev(height = 0.5),
    #     stroke = 1.5) +
    #   ggplot2::scale_shape_manual(values = c(21, 22, 23, 24, 25)) + # Customize shape values as needed
    #   ggplot2::scale_size_continuous(range = c(3, 10)) + # Adjust the range for point sizes
    #   ggplot2::scale_color_manual(values = custom_colors) + # Add custom colors for points
    #   ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = ssEnv$color_palette[1], size = 1) + # Add vertical line at FDR = 1
    #   labs(
    #     title = "",
    #     x = '-log10(P_Value)',
    #     y = 'Description'
    #   ) +
    #   ggplot2::guides(shape = FALSE) + # Hide the shape legend
    #   ggplot2::theme_minimal() +
    #   ggplot2::theme(
    #     axis.title.x = ggplot2::element_text(size = 12),
    #     axis.title.y = ggplot2::element_text(size = 12),
    #     axis.text.x = ggplot2::element_text(size = 10),
    #     axis.text.y = ggplot2::element_text(size = 10),
    #     plot.title = ggplot2::element_text(size = 14, face = 'bold')
    #   )

    # Create the lollipop plot with vertical dodge and shapes based on Enrichment
    # ggplot(data_to_plot, aes(x = log_fdr, y = Description, size = Enrichment, color = key)) +
    #   geom_point(aes(shape = factor(cut(Enrichment, breaks = 5, labels = FALSE))),
    #     fill = NA,
    #     position = position_dodgev(height = 0.5),
    #     stroke = 1.5) +
    #   ggplot2::scale_shape_manual(values = c(21, 22, 23, 24, 25)) + # Customize shape values as needed
    #   ggplot2::scale_size_continuous(range = c(3, 10)) + # Adjust the range for point sizes
    #   ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) + # Add vertical line at FDR = 1
    #   labs(
    #     title = "",
    #     x = '-log10(P_Value)',
    #     y = 'Description'
    #   ) +
    #   ggplot2::guides(shape = FALSE) + # Hide the shape legend
    #   ggplot2::theme_minimal() +
    #   ggplot2::theme(
    #     axis.title.x = ggplot2::element_text(size = 12),
    #     axis.title.y = ggplot2::element_text(size = 12),
    #     axis.text.x = ggplot2::element_text(size = 10),
    #     axis.text.y = ggplot2::element_text(size = 10),
    #     plot.title = ggplot2::element_text(size = 14, face = 'bold')
    #   )

    # # Create the lollipop plot with vertical dodge and shapes based on Enrichment
    # ggplot(data_to_plot, aes(x = log_fdr, y = Description, size = Enrichment, color = key,)) +
    #   geom_point(aes(size = Enrichment), fill = NA, position = position_dodgev(height = 0.5), stroke = 1.5) +
    #   ggplot2::scale_shape_manual(values = c(21, 22, 23, 24,25)) + # Customize shape values as needed
    #   ggplot2::scale_size_continuous(range = c(3, 10)) + # Adjust the range for point sizes
    #   labs(
    #     title = "",
    #     x = '-log10(P_Value)',
    #     y = 'Description'
    #   ) +
    #   ggplot2::guides(shape = FALSE) + # Hide the shape legend
    #   ggplot2::theme_minimal() +
    #   ggplot2::theme(
    #     axis.title.x = ggplot2::element_text(size = 12),
    #     axis.title.y = ggplot2::element_text(size = 12),
    #     axis.text.x = ggplot2::element_text(size = 10),
    #     axis.text.y = ggplot2::element_text(size = 10),
    #     plot.title = ggplot2::element_text(size = 14, face = 'bold')
    #   )

    fname <- core_name_cleaning(paste(file_prfx,"_lollipop_plot_",c,"_",rules["label"],ifelse(disease=="","", paste("_", disease, sep="")), sep="_"))
    filename <- paste(path, "/",fname, ".",ssEnv$plot_format,sep = "")
    ggplot2::ggsave(filename,  width = 16, height = 9, dpi = as.numeric(ssEnv$plot_resolution_ppi))
  }
}
