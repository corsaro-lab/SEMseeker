#' @importFrom doRNG %dorng%
#' @importFrom doFuture %dofuture%
# NOTE: currently internal — could be re-exported once @examples are added
marker_quantization_metric <- function(result_folder, maxResources = 90, parallel_strategy  = "multisession", ...)
{
  # result_folder, maxResources = 90, parallel_strategy  = "multisession", ...
  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)
  # ssEnv <- get_session_info()

  keys <- ssEnv$keys_areas_subareas_markers_figures
  keys <- keys[keys$AREA == "PROBE",]
  if(nrow(keys) == 0)
  {
    log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y"), " For this analisys the PROBE area is required")
    close_env()
    return()
  }


  keys <- keys[keys$MARKER != "DELTAS",]
  keys <- keys[keys$MARKER != "DELTAR",]
  keys <- keys[keys$MARKER != "LESIONS",]
  keys <- keys[keys$MARKER != "SIGNAL",]

  keys <- keys[complete.cases(keys),]
  nkeys <- nrow(keys)

  #
  #

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nkeys)
  else
    progress_bar <- ""


  to_export <- c("keys","result_folderPivot","ssEnv","file_path_build","progress_bar", "progression_index","progression","progressor_uuid","owner_session_uuid","trace")
  result_temp <- data.frame()
  scores <- data.frame()


  result_temp <- foreach::foreach(k = 1:nkeys, .combine =  plyr::rbind.fill, .export = to_export) %dorng%
    # for (k in 1:nkeys)
    {
      # workers must NOT saveRDS on every iteration.
      update_session_info(ssEnv, save_to_disk = FALSE)

      # k <- 1
      key <- keys [k,]
      log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y"), "Processing key: ", key$MARKER," ", key$FIGURE," ", key$AREA," ", key$SUBAREA)


      res_temp <- data.frame("MARKER" = key$MARKER, "FIGURE" = key$FIGURE, "AREA" = key$AREA, "SUBAREA" = key$SUBAREA)

      if(is.na(key$MARKER))
      {
        return(NULL)
      }

      if(key$MARKER=="DELTARP" | key$MARKER=="DELTARQ")
        original_marker <- "DELTAR"
      else
        original_marker <- "DELTAS"


      res_temp$ORIGINAL_MARKER <- original_marker

      pivot_subfolder <- dir_check_and_create(result_folderPivot, original_marker)
      fname <- file_path_build( pivot_subfolder ,c(original_marker, key$FIGURE, key$AREA,key$SUBAREA),"csv", add_gz=TRUE)
      if(!file.exists(fname))
      {
        log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y")," File not found: ", fname)
        return(NULL)
      }
      original <- as.matrix(utils::read.csv2(fname, header = TRUE, row.names = 1, skip = 1))


      pivot_subfolder <- dir_check_and_create(result_folderPivot, key$MARKER)
      fname <- file_path_build( pivot_subfolder ,c(key$MARKER, key$FIGURE, key$AREA,key$SUBAREA),"csv", add_gz=TRUE)
      if(!file.exists(fname))
      {
        log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y")," File not found: ", fname)
        return(NULL)
      }
      quantized <- as.matrix(utils::read.csv2(fname, header = TRUE, row.names = 1, skip = 1))

      tryCatch(
        {
          mdl_perf <- model_performance(original, quantized,c(),c())
          res_temp <- cbind(res_temp, mdl_perf)

          # SSIM con pracma
          ssim_value <- ssim(original, quantized)
          res_temp$Structural_Similarity_Index <- ssim_value

          # vi_value <- variation_of_information(original, quantized)
          # res_temp$variation_of_information <- vi_value

          # calculate JSD
          # Combine the unique elements from both samples to create a common event space
          common_events <- unique(c(original, quantized))

          # Create adjusted frequency tables for both samples
          frequency_table1_adjusted <- tabulate(match(original, common_events), nbins = length(common_events))
          frequency_table2_adjusted <- tabulate(match(quantized, common_events), nbins = length(common_events))

          # Convert adjusted frequencies to probabilities
          probability_distribution1_adjusted <- frequency_table1_adjusted / sum(frequency_table1_adjusted)
          probability_distribution2_adjusted <- frequency_table2_adjusted / sum(frequency_table2_adjusted)

          # Calculate the Jensen-Shannon distance
          jsd <- suppressMessages(suppressWarnings(philentropy::JSD(rbind(probability_distribution1_adjusted, probability_distribution2_adjusted))))
          res_temp$JSD <- jsd
        },
        error = function(e)
        {
          log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y"), "Error in calculate_metrics: ", e)
        }
      )

      if(ssEnv$showprogress)
        progress_bar(sprintf("I'm doing metrics calculation for the quantization process."))

      res_temp
      # result_temp <- plyr::rbind.fill(result_temp, res_temp)

      # Errore Assoluto Medio (MAE)
      # •   Willmott, C. J., & Matsuura, K. (2005). Advantages of the mean absolute error (MAE) over the root mean square error (RMSE) in assessing average model performance. Climate Research, 30(1), 79-82.
      # 2.  Indice di Similarità Strutturale (SSIM)
      # •   Wang, Z., Bovik, A. C., Sheikh, H. R., & Simoncelli, E. P. (2004). Image quality assessment: From error visibility to structural similarity. IEEE Transactions on Image Processing, 13(4), 600-612.
      # •   Avanaki, M. R. N. (2009). Exact calculation of the structural similarity index for image quality assessment. Journal of Mathematical Imaging and Vision, 34(2), 121-131.
      # 3.  Variazione dell’Informazione (VI)
      # •   Meilă, M. (2007). Comparing clusterings—an information based distance. Journal of Multivariate Analysis, 98(5), 873-895.
      # •   Vinh, N. X., Epps, J., & Bailey, J. (2010). Information theoretic measures for clusterings comparison: Variants, properties, normalization and correction for chance. Journal of Machine Learning Research, 11(10), 2837-2854.


    }

  # post-foreach end-of-batch snapshot (matches pattern).
  update_session_info(ssEnv, save_to_disk = TRUE)

  colnames(result_temp) <- toupper(colnames(result_temp))
  dataFolder <- dir_check_and_create(ssEnv$result_folderData,c("Distributions"))
  filename  =  file_path_build(dataFolder,c("DISTRIBUTION", "ANALYSIS"),"csv")
  utils::write.csv2(result_temp, file = filename, row.names = FALSE)


  cols_to_cycke <- toupper(c("JSD", "Structural_Similarity_Index", "variation_of_information","MAPE","R-SQUARED"))
  cols_to_cycke <- which(colnames(result_temp) %in% cols_to_cycke)
  scores <- data.frame()
  keys <- unique(result_temp[,c("FIGURE","MARKER")])
  for(metric in cols_to_cycke)
  {
    metric <- colnames(result_temp)[metric]
    scores_temp <- data.frame()
    scores_temp <- metrics_ranking(metric,result_temp, column_to_rank =metric)
    scores_temp$FIGURE <- result_temp$FIGURE
    scores <- plyr::rbind.fill(scores, scores_temp)
  }
  # for (k in 1:nrow(keys))
  # {
  #   key <- keys[k,]
  #   res <- result_temp[result_temp$FIGURE == key$FIGURE & result_temp$MARKER == key$MARKER,]
  #   fig <- key$FIGURE
  #   for(metric in cols_to_cycke)
  #   {
  #     metric <- colnames(res)[metric]
  #     scores <- metrics_ranking(metric = metric, data_frame = res,figure = fig, scores = scores)
  #   }
  # }

  scores$SCORE <- round(scores$SCORE, 2)
  # aggregate scores by MARKER and sum RANK
  scores_agg <- aggregate(scores$SCORE, by = list(scores$MARKER), FUN = sum)
  # calculate TOTAL per marker
  scores_agg <- scores_agg[order(scores_agg$x, decreasing = TRUE),]
  colnames(scores_agg) <- c("MARKER","TOTAL")

  # create a pivot table with MARKER as rows and FIGURE as columns
  scores_agg_fig <- reshape2::dcast(scores, MARKER ~ FIGURE, value.var = "SCORE", fun.aggregate = sum)
  scores_agg_fig <- merge(scores_agg_fig,scores_agg, by="MARKER")
  # sort by SCORE descending
  scores_agg_fig <- scores_agg_fig[order(scores_agg_fig$TOTAL, decreasing = TRUE),]
  # save scores
  filename  =  file_path_build(dataFolder,c("DISTRIBUTION", "ANALYSIS","SCORE"),"csv")
  utils::write.csv2(scores_agg_fig, file = filename, row.names = FALSE)

}



