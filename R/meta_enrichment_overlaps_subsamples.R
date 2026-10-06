#' @title Enrichment stability across random subsamples of a study
#' @description Tests the stability of enrichment (pathway/term) results under
#'   random subsampling of a single cohort, quantifying how robust the enriched
#'   terms are to sample variability. The intra-study counterpart of
#'   [meta_enrichment_compare_studies()].
#' @param inference_details Data frame describing the enrichment analyses to run
#'   (one row per analysis).
#' @param pathways_sql_selection Optional SQL WHERE fragment to restrict the
#'   enriched terms considered.
#' @param old_label_samples_sql_condition,new_label_samples_sql_condition SQL
#'   conditions selecting the two sample sets to contrast.
#' @param old_label_association_results_sql_condition,new_label_association_results_sql_condition
#'   SQL conditions selecting the association results to contrast.
#' @param run_prefix Optional prefix for the run outputs.
#' @param pathway_package Enrichment backend to use (e.g. `"WebGestalt"`);
#'   `""` compares all available backends.
#' @param association_pvalue_column Name of the p-value column driving
#'   significance.
#' @param significance Logical; keep only significant terms when `TRUE`.
#' @param result_folder Path to the folder holding the study results.
#' @param ... Additional arguments forwarded to [core_init_env()].
#' @return Invisibly `NULL`; stability tables and plots are written under
#'   `result_folder`.
#' @seealso [meta_enrichment_compare_studies()]
#' @examples
#' # See vignette("pathway-analysis", package = "SEMseeker") for a runnable
#' # enrichment-stability (subsampling) workflow on real result folders.
#' invisible(NULL)
#' @export
meta_enrichment_overlaps_subsamples <- function(inference_details,pathways_sql_selection="",
  old_label_samples_sql_condition = "", new_label_samples_sql_condition = "",
  old_label_association_results_sql_condition = "", new_label_association_results_sql_condition = "",
  run_prefix = "",pathway_package="", association_pvalue_column = "", significance=TRUE,
  result_folder, ...)
{

  # join all samples_sql_condition
  samples_sql_folder <- core_name_cleaning(paste0(unique(strsplit(paste0(core_name_cleaning(sort(inference_details$samples_sql_condition)), collapse = "_"),"_")[[1]]),collapse = "_"))

  association_pvalue_column <- core_name_cleaning(association_pvalue_column)
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

  core_log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), " SEMseeker will perform the pathway analysys for sub samples overlaps.")

  color_palette <- ssEnv$color_palette
  localKeys <- unique(ssEnv$keys_for_pathway)

  pathway_results <- data.frame()
  association_results_sql_conditions <- unique(inference_details$association_results_sql_condition)
  key_enrichment_format <- ssEnv$key_enrichment_format
  key_enrichment_format <- subset(key_enrichment_format, key_enrichment_format[,"label"]!="phenolyzer")
  if(pathway_package!="")
    key_enrichment_format <- subset(key_enrichment_format, key_enrichment_format[,"label"]==pathway_package)

  for (pt in seq_len(nrow(key_enrichment_format)))
  {
    column_of_id <- key_enrichment_format[pt,"column_of_id"]
    column_of_adj_pvalue <- key_enrichment_format[pt,"column_of_adj_pvalue"]
    column_of_description <- key_enrichment_format[pt,"column_of_description"]
    column_of_enrichment <- key_enrichment_format[pt,"column_of_enrichment"]
    enrichment_package <- key_enrichment_format[pt,"label"]
    signal_suffixes <-c("","with_signal","without_signal")
    for (s in 1:3)
    {
      for (z in seq_along(association_results_sql_conditions))
      {
        association_results_sql_condition <- core_name_cleaning(association_results_sql_conditions[z])
        for (a in seq_len(nrow(localKeys)))
        {
          # for each SAMPLES_SQL_CONDITION in inference_details
          for (i in seq_len(nrow(inference_details)))
          {
            inference_detail <- inference_details[i,]
            enrich_phenotype_analysis_name <- enrich_phenotype_analysis_name(inference_detail, localKeys[a,],prefix ="", suffix= signal_suffixes[s] , association_pvalue_column, ssEnv$alpha, significance)
            path <- io_dir_check_and_create(ssEnv$result_folderEnrichment,c(enrichment_package, core_name_cleaning(inference_detail$areas_sql_condition), core_name_cleaning(inference_detail$samples_sql_condition), core_name_cleaning(association_results_sql_condition)))
            pathway_report_path <- io_file_path_build(path,enrich_phenotype_analysis_name,"csv")
            if(file.exists(pathway_report_path))
            {
              temp_res <- utils::read.csv2(pathway_report_path)
              if(nrow(temp_res) != 0)
                temp_res <- assoc_filter_sql(pathways_sql_selection,temp_res)

              if(nrow(temp_res) != 0)
              {
                temp_res$MARKER <- localKeys[a, "MARKER"]
                temp_res$FIGURE <- localKeys[a, "FIGURE"]
                temp_res$AREA <- localKeys[a, "AREA"]
                temp_res$SUBAREA <- localKeys[a, "SUBAREA"]
                temp_res$SAMPLES_SQL_CONDITION <- core_name_cleaning(inference_details[i,"samples_sql_condition"])
                temp_res$association_results_sql_condition <- core_name_cleaning(inference_details[i,"association_results_sql_condition"])
                ff <- new_label_samples_sql_condition
                gg <- old_label_samples_sql_condition
                for (j in seq_along(ff))
                {
                  temp_res[core_name_cleaning(temp_res$SAMPLES_SQL_CONDITION)==core_name_cleaning(gg[j]),"SAMPLES_SQL_CONDITION"] <- core_name_cleaning(ff[j])
                }
                ff <- new_label_association_results_sql_condition
                gg <- old_label_association_results_sql_condition
                for (j in seq_along(ff))
                {
                  temp_res[core_name_cleaning(temp_res$association_results_sql_condition)==core_name_cleaning(gg[j]),"association_results_sql_condition"] <- core_name_cleaning(ff[j])
                }
                temp_res$SAMPLES_SQL_CONDITION[temp_res$SAMPLES_SQL_CONDITION == ""] <- "ALL"
                temp_res$association_results_sql_condition[temp_res$association_results_sql_condition == ""] <- "ALL"
                pathway_results <- plyr::rbind.fill(pathway_results, temp_res)
              }
            }
          }
        }

        if(nrow(pathway_results)==0)
          next

        pathway_results <- subset(pathway_results, pathway_results[,column_of_adj_pvalue]<=as.numeric(ssEnv$alpha))
        # Using with() for cleaner code
        pathway_results$KEY_SELECTOR <- with(pathway_results,
          as.factor(paste0(SAMPLES_SQL_CONDITION,"_",association_results_sql_condition))
        )
        core_log_event("INFO: ",format(Sys.time(), "%a %b %d %X %Y")," pathway files aggregated!")

        keys <- unique(pathway_results[, c("SUBAREA", "AREA", "MARKER", "FIGURE")])
        for (i in seq_len(nrow(keys)))
        {
          # i <- 1
          pathway_set <- pathway_results[
            pathway_results$AREA == keys[i, ]$AREA &
              pathway_results$SUBAREA == keys[i, ]$SUBAREA &
              pathway_results$MARKER == keys[i, ]$MARKER &
              pathway_results$FIGURE == keys[i, ]$FIGURE
            , c(column_of_id, "KEY_SELECTOR"), ]

          if(nrow(pathway_set)==0)
            next
          # replace samples_sql_condition with alternative names
          pathway_set$KEY_SELECTOR <- as.character(pathway_set$KEY_SELECTOR)
          pathway_set <- na.omit(pathway_set)
          SPLIT <- split(pathway_set[,column_of_id], pathway_set$KEY_SELECTOR)

          categories <- unique(pathway_set$KEY_SELECTOR)
          # put ALL where categories==""
          categories <- ifelse(categories == "", "ALL", categories)
          if(length(categories)<2)
            next

          dest_folder <- io_dir_check_and_create(ssEnv$result_folderChart,c("PATHWAYS_CROSS_SAMPLE_ANALYSIS",association_results_sql_condition,
            core_name_cleaning(enrichment_package),core_name_cleaning(pathways_sql_selection),core_name_cleaning(samples_sql_folder)))
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
              core_name_cleaning(column_of_adj_pvalue),
              "_",
              core_name_cleaning(ssEnv$alpha),
              "_",
              core_name_cleaning(family_test),
              "_UPSET_PLOT.",
              ssEnv$plot_format,
              sep = ""
            )

          pathway_set <- UpSetR::fromList(SPLIT)
          if (is.null(dim(pathway_set)))
            next

          # sum by row
          degree <- rowSums(pathway_set[,categories])
          intersect_length <- nrow(unique(pathway_set))
          intersect_length <- min(intersect_length, 30)
          # Save the plot to a file (e.g., PNG)

          if (any(degree == length(categories)))
            # upset(pathway_set, sets = sets, decreasing = c(FALSE,TRUE), nintersects = 30)
            my_plot <- UpSetR::upset(data = pathway_set, sets = categories, decreasing = c(FALSE,TRUE),
              nintersects = 30,
              order.by = c("freq", "degree"),
              main.bar.color = ssEnv$color_palette[1], matrix.color = ssEnv$color_palette[2], sets.bar.color = ssEnv$color_palette[3],
              queries = list(list(query = UpSetR::intersects, params = list(categories), color = ssEnv$color_palette[4])))
          else
            my_plot <- UpSetR::upset(data = pathway_set, sets = categories, decreasing = c(FALSE,TRUE), nintersects = 30,
              main.bar.color = ssEnv$color_palette[1], matrix.color = ssEnv$color_palette[2], sets.bar.color = ssEnv$color_palette[3])

          # First, open a file device for plotting (PNG, PDF, etc.)
          grDevices::png(filename, width = 9, height = 9, units="in", res = as.numeric(ssEnv$plot_resolution_ppi))
          print(my_plot)
          # Close the file device to save the plot
          grDevices::dev.off()

          core_log_event("DEBUG: ",format(Sys.time(), "%a %b %d %X %Y"),"  venn diagram completed !")

          dest_folder <- io_dir_check_and_create(ssEnv$result_folderEnrichment,c("PATHWAYS_CROSS_SAMPLE_ANALYSIS",association_results_sql_condition,
            core_name_cleaning(enrichment_package),core_name_cleaning(pathways_sql_selection),core_name_cleaning(samples_sql_folder)))
          inference_detail$samples_sql_condition <- ""
          overlaps <- Reduce(intersect, SPLIT)
          if(length(overlaps)>0)
          {
            filename <- io_inference_file_name(inference_detail, paste(keys[i, ]$AREA, keys[i, ]$SUBAREA, keys[i, ]$MARKER, keys[i, ]$FIGURE , sep="_"),dest_folder,file_extension = "csv",suffix = "OVERLAPS", prefix = "")
            overlaps <- data.frame(overlaps)
            colnames(overlaps) <- column_of_id
            overlaps <- merge(overlaps, pathway_results, by = column_of_id)
            utils::write.csv2(overlaps, filename)
          }

          # create a pivot table
          if (length(categories) > 1)
          {
            filename <- io_inference_file_name(inference_detail, paste(keys[i, ]$AREA, keys[i, ]$SUBAREA, keys[i, ]$MARKER, keys[i, ]$FIGURE , sep="_"),dest_folder,file_extension = "csv",suffix = "PIVOT", prefix = "")
            # Create a simple pivot table
            # Columns named as strings: the count column dplyr::count() adds is
            # "n", and naming it here rather than reading a bare n keeps the
            # pipeline free of names R CMD check cannot bind.
            pivot_table <- pathway_results %>%
              dplyr::count(dplyr::across(dplyr::all_of(c(column_of_id, "KEY_SELECTOR")))) %>%
              tidyr::pivot_wider(
                names_from = "KEY_SELECTOR",        # one column per key
                values_from = "n",                  # filled with the count
                values_fill = 0                     # missing combinations are 0
              )
            colnames(pivot_table)[1] <- column_of_id
            pivot_table <- merge(unique(pathway_results[,c(column_of_description,column_of_id)]),pivot_table, by = column_of_id, all.y = TRUE)
            utils::write.csv2(pivot_table, filename)
          }

          core_log_event("INFO: ",format(Sys.time(), "%a %b %d %X %Y"),"  job completed !")
        }
      }
    }
  }
  core_close_env()
}
