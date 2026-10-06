#' Diagnostic performance (sensitivity & specificity) of SEMseeker mutations and lesions
#'
#' @param samples_sql_selection character. SQL-style filter expression to
#'   restrict which samples are included (default \code{""}, all samples).
#' @param combinations list of character vectors. Each element is a vector of
#'   sample group labels to compare (e.g. \code{list(c("CASE_A", "CASE_B"))}).
#' @param result_folder character. Path to the SEMseeker result folder.
#' @param independent_variable character. Sample sheet column defining the
#'   grouping variable (default \code{"Sample_Group"}).
#' @param maxResources numeric. Maximum percentage of CPU cores to use
#'   (default 90).
#' @param parallel_strategy character. Parallelisation backend; one of
#'   \code{"multisession"}, \code{"sequential"}, \code{"cluster"}
#'   (default \code{"multisession"}).
#'   Asking for \code{"multicore"} is accepted and converted to
#'   \code{"multisession"}: it means fork(), which is unsafe with this
#'   package's native thread pool on every platform that offers it, and
#'   absent on Windows. The conversion is logged.
#' @param ... Additional arguments passed to \code{core_init_env()}.
#'
#' @return Invisibly \code{NULL}. Sensitivity / specificity tables are written
#'   to the \code{Euristic/} sub-folder of \code{result_folder}.
#' @importFrom doRNG %dorng%
#' @examples
#' result_dir <- tempdir()
#' \dontrun{
#' diagnostic_performance(
#'   combinations         = list(c("CASE_A", "CASE_B")),
#'   result_folder        = "~/semseeker_results/",
#'   independent_variable = "Sample_Group"
#' )
#' }
#' @export
diagnostic_performance <-
  function(samples_sql_selection="",combinations,result_folder,independent_variable = "Sample_Group",maxResources = 90,parallel_strategy  = "multisession",...)
  {
    j <- 0
    k <- 0
    z <- 0
    arguments <- list(...)
    ssEnv <- core_init_env(result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)


    arguments <- list(...)
    areas_selection <- c()
    if (!is.null(arguments[["areas_selection"]]))
    {
      areas_selection <- arguments$areas_selection
    }

    localKeys <- ssEnv$keys_markers_figures
    sample_groups <- c("Reference", "Control", "Case")

    study_summary <-   sem_study_summary_get()
    study_summary <- assoc_filter_sql(samples_sql_selection, study_summary)

    study_summary <- study_summary[, c("Sample_ID",independent_variable)]

    markers <- unique(ssEnv$keys_markers_figures$MARKER)
    # as vector of character
    result_colnames <- c("MARKER", "FIGURE", "AREA", "SUBAREA", "AREA_OF_TEST","Case","Control", "SENSITIVITY", "SPECIFICITY",
      "P_to_be_Case_cond_to_be_Epimutated","P_to_be_Control_cond_to_be_Not_Epimutated","SCORE","BURDEN","JSD")
    for (a in seq_along(markers))
    {
      # a <- 2
      dest_path <- io_dir_check_and_create(ssEnv$result_folderEuristic, samples_sql_selection)
      fileNameResults <- io_file_path_build(baseFolder =  dest_path,detailsFilename =  c( as.character(markers[a]),independent_variable, "sensitivity_specificity_analisys"),extension = "csv")
      fileNameResultsTemp <- io_file_path_build(baseFolder =  dest_path,detailsFilename =  c( as.character(markers[a]),independent_variable,"sensitivity_specificity_analisys","tmp"),extension = "csv")

      localKeys_1 <- ssEnv$keys_areas_subareas_markers_figures
      keys <- localKeys_1[localKeys_1$MARKER == markers[a], ]


      # clean keys from already done analysis
      if (file.exists(fileNameResults))
      {
        results <- utils::read.csv2(fileNameResults, header  =  TRUE, sep=";")
        keys_markers_figures_areas_done <- unlist(apply(unique(results [, c("MARKER", "FIGURE", "AREA", "SUBAREA")]), 1, function(x) paste(x, collapse  =  "_", sep  =  "")))
        keys_to_be_done <-
          unlist(apply(keys[, c("MARKER", "FIGURE", "AREA", "SUBAREA")], 1, function(x)
            paste(
              x, collapse  =  "_", sep  =  ""
            )))
        keys <- keys[!(keys_to_be_done %in% keys_markers_figures_areas_done),]
      }

      nkeys <- nrow(keys)
      if (nrow(keys) > 0)
        for (k in seq_len(nkeys))
        {
          # k <- 3
          key <- keys [k, ]
          # the diagnostic reads the same artefact the inference tested,
          # so it has to name the same aggregation. Defaulting to the produced
          # one keeps the historical behaviour for callers that do not pass an
          # inference detail.
          aggregation <- util_aggregations_allowed(key$MARKER, key$FIGURE,
                                                   discrete = isTRUE(key$DISCRETE),
                                                   default  = TRUE,
                                                   scope    = "INSTANCE",
                                                   area     = key$AREA)[1]
          fname <- io_pivot_file_name_parquet(key$MARKER, key$FIGURE, key$AREA,
                                              key$SUBAREA,
                                              aggregation = aggregation,
                                              scope = "INSTANCE")
          if (file.exists(fname))
          {
            core_log_event("INFO: ",
              format(Sys.time(), "%a %b %d %X %Y"),
              " Starting to read pivot:",
              fname,
              ".")
            # pivot <- readr::read_delim(io_pivot_file_name,
            #   col_types = readr::cols(
            #     .default = readr::col_double(),
            #     AREA = readr::col_character(),
            #   ),
            #   show_col_types=FALSE, progress=FALSE)
            # read via unified dispatcher.
            ss_pivot_lazy <- io_read_pivot(key$MARKER, key$FIGURE, key$AREA, key$SUBAREA)
            tempDataFrame <- as.data.frame(ss_pivot_lazy$collect())
            row.names(tempDataFrame) <- tempDataFrame$AREA
            if (length(areas_selection) > 0)
              tempDataFrame <-tempDataFrame[tempDataFrame[, 1] %in% areas_selection,]

            #removes area_of_test name (eg. gene...)
            tempDataFrame <- tempDataFrame[, -1]

            if (is.null(dim(tempDataFrame)))
              next
            if (plyr::empty(tempDataFrame) | nrow(tempDataFrame) == 0)
              next
            tempDataFrame <- t(tempDataFrame)
            tempDataFrame <- as.data.frame(tempDataFrame)

            core_log_event(
              "INFO: ",
              format(Sys.time(), "%a %b %d %X %Y"),
              " Read pivot:",
              fname,
              " with ",
              ncol(tempDataFrame),
              " rows."
            )
            tempDataFrame$Sample_ID <- rownames(tempDataFrame)
            tempDataFrame <- merge(study_summary,tempDataFrame, by  =  "Sample_ID", all.x=TRUE)
            tempDataFrame <- tempDataFrame[,colnames(tempDataFrame) != "Sample_ID"]
            # check only two values are in the independent variable
            # if (length(unique(tempDataFrame[, independent_variable]))!=2)
            #   stop("ERROR: I'm stopping here, the grouping variable should be only two!")
            tempDataFrame <- as.data.frame(tempDataFrame)
            tempDataFrame[is.na(tempDataFrame)] <- 0
            # 2026-06-09: same sanitize+memo+counter-rename pattern
            # as apply_stat_model.R. Colnames hold AREA_OF_TEST gene/CpG
            # island names that may carry ' ', '-', ':', '/', "'" - all
            # invalid as R identifiers. Sanitise here only for internal
            # formula safety; the result's AREA_OF_TEST is reverse-mapped
            # back to the raw name before writing the CSV (preserves the
            # upstream-annotated identifier for enrichment + resume match).
            real_cols <- colnames(tempDataFrame)
            safe_cols <- gsub("[^A-Za-z0-9_.]", "_", real_cols)
            if (anyDuplicated(safe_cols)) {
              safe_cols <- make.unique(safe_cols, sep = "_")
            }
            safe_to_real <- stats::setNames(real_cols, safe_cols)
            colnames(tempDataFrame) <- safe_cols

            tempDataFrame[,independent_variable] <- as.character(tempDataFrame[,independent_variable])

            labels <- unique(tempDataFrame[,independent_variable])
            # get all the combinations of two over all labels
            # combinations <- combinat::combn(labels, 2)
            for ( comb in seq_len(ncol(combinations)))
            {
              control_label <- combinations[rownames(combinations)=="Control",comb]
              case_label <- combinations[rownames(combinations)=="Case",comb]
              # control_label <- as.character(control_label)
              # case_label <- as.character(case_label)

              # reduce the dataframe to the two labels
              tempDataFrameComb <- tempDataFrame[tempDataFrame[,independent_variable] %in% c(control_label, case_label),]
              actual_labels <- tempDataFrameComb[,independent_variable] == case_label
              phenotype <- tempDataFrameComb[,independent_variable] == case_label

              n_case <- sum(phenotype)
              n_control <- nrow(tempDataFrameComb) - n_case

              if(ssEnv$showprogress)
                progress_bar <- progressr::progressor(along = 2:ncol(tempDataFrameComb))
              else
                progress_bar <- ""

              # 2026-06-09: safe_to_real reaches each worker so the
              # AREA_OF_TEST written in the result_temp data.frame is the
              # raw upstream name (no '-'→'_' rewrite).
              var_to_export <- c("tempDataFrameComb","ssEnv","progress_bar","actual_labels","keys","k","progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace","safe_to_real")
              # for (c in 2:ncol(tempDataFrameComb))
              results_temp <- foreach::foreach(c  =  2:ncol(tempDataFrameComb), .combine  =  rbind, .export  =  var_to_export) %dorng%
              {
                # workers must NOT saveRDS on every iteration.
                core_update_session_info(ssEnv, save_to_disk = FALSE)
                area_of_test <- names(tempDataFrameComb)[c]
                if(ssEnv$showprogress)
                {
                  progress_bar(sprintf("genomic area of test: %s", stringr::str_pad( area_of_test, 20, side=c('left'), pad=' ')))
                }
                tempDataFrameComb[,c] <- as.numeric(tempDataFrameComb[,c])
                tempDataFrameComb[,c] <- ifelse(is.na(tempDataFrameComb[,c]), 0, tempDataFrameComb[,c])
                predicted_labels <- tempDataFrameComb[,c]!=0
                burden <- sum(as.numeric(tempDataFrameComb[,c]))

                # Create a confusion matrix
                conf_matrix <- table(Actual = actual_labels, Predicted = predicted_labels)
                if(sum(predicted_labels) == 0)
                  conf_matrix <- cbind(conf_matrix, "TRUE"=0)
                if(sum(predicted_labels) == length(predicted_labels))
                  conf_matrix <- cbind(conf_matrix, "FALSE"=0)

                # check if conf matrix as 2x2 shape
                if (nrow(conf_matrix) != 2 | ncol(conf_matrix) != 2)
                # Calculate sensitivity and specificity
                sensitivity <- conf_matrix["TRUE","TRUE"] / sum(conf_matrix["TRUE", ])
                specificity <- conf_matrix["FALSE", "FALSE"] / sum(conf_matrix["FALSE", ])

                epimutated <- as.logical(tempDataFrameComb[,c]!=0)
                # apply bayes theorem
                # P(A|B) = P(B|A) * P(A) / P(B)
                # P(A) = probability to be a case
                # P(B) = probability to have a mutation in a specific area_of_test
                # P(B|A) = probability to have a mutation in a specific area_of_test given that the sample is a case
                # P(A|B) = probability to be a case given that the sample has a mutation in a specific area_of_test
                P_to_be_Epimutated <- sum(epimutated) / length(epimutated)

                P_to_be_Case <- n_case / (n_case + n_control)
                P_to_be_Epimutated_cond_to_be_Case <- sum(epimutated & phenotype) / sum(phenotype)
                if (P_to_be_Epimutated == 0)
                  P_to_be_Case_cond_to_be_Epimutated <- 0
                else
                  P_to_be_Case_cond_to_be_Epimutated <- (P_to_be_Epimutated_cond_to_be_Case * P_to_be_Case) / P_to_be_Epimutated

                P_to_be_Control <- n_control / (n_case + n_control)
                P_to_be_Not_Epimutated_cond_to_be_Control <- sum(!epimutated & (!phenotype)) / sum(!phenotype)
                P_to_be_Control_cond_to_be_Not_Epimutated <- (P_to_be_Not_Epimutated_cond_to_be_Control * P_to_be_Control) / (1 -P_to_be_Epimutated)

                sample1 <- as.numeric(tempDataFrameComb[phenotype,c])
                sample2 <- as.numeric(tempDataFrameComb[!phenotype,c])

                # Combine the unique elements from both samples to create a common event space
                common_events <- unique(c(sample1, sample2))

                # Create adjusted frequency tables for both samples
                frequency_table1_adjusted <- tabulate(match(sample1, common_events), nbins = length(common_events))
                frequency_table2_adjusted <- tabulate(match(sample2, common_events), nbins = length(common_events))

                # Convert adjusted frequencies to probabilities
                probability_distribution1_adjusted <- frequency_table1_adjusted / sum(frequency_table1_adjusted)
                probability_distribution2_adjusted <- frequency_table2_adjusted / sum(frequency_table2_adjusted)

                # Calculate the Jensen-Shannon distance
                jsd <- suppressMessages(suppressWarnings(philentropy::JSD(rbind(probability_distribution1_adjusted, probability_distribution2_adjusted))))

                # data.frame("MARKER"= as.character(keys[k, "MARKER"]),"FIGURE"= as.character(keys[k, "FIGURE"]),
                #   "AREA"= as.character(keys[k, "AREA"]), "SUBAREA"= as.character(keys[k, "SUBAREA"]),
                #   "AREA_OF_TEST" = area_of_test,"SENSITIVITY"=sensitivity, "SPECIFICITY"=specificity,
                #   "P_to_be_Case_cond_to_be_Epimutated"= P_to_be_Case_cond_to_be_Epimutated,
                #   "P_to_be_Control_cond_to_be_Not_Epimutated"= P_to_be_Control_cond_to_be_Not_Epimutated,
                #   "SCORE"=0,"BURDEN"=burden,"JSD"=jsd)

                # replace Na or Nan with 0
                # 2026-06-09: reverse-map AREA_OF_TEST back to the
                # raw upstream name so the CSV preserves it for enrichment
                # / resume match. Defensive fallback if the mapping is
                # missing (should not happen).
                area_of_test_raw <- if (area_of_test %in% names(safe_to_real)) {
                  safe_to_real[[area_of_test]]
                } else {
                  area_of_test
                }
                results_temp <- data.frame("MARKER"= as.character(keys[k, "MARKER"]),"FIGURE"= as.character(keys[k, "FIGURE"]),
                  "AREA"= as.character(keys[k, "AREA"]), "SUBAREA"= as.character(keys[k, "SUBAREA"]),
                  "AREA_OF_TEST" = area_of_test_raw, "Case"= case_label, "Control"= control_label,
                  "SENSITIVITY"=sensitivity, "SPECIFICITY"=specificity,
                  "P_to_be_Case_cond_to_be_Epimutated"= P_to_be_Case_cond_to_be_Epimutated,
                  "P_to_be_Control_cond_to_be_Not_Epimutated"= P_to_be_Control_cond_to_be_Not_Epimutated,
                  "SCORE"=0,"BURDEN"=burden,"JSD"=jsd)
                # results_temp[is.na(results_temp) | is.nan(results_temp)] <- 0
                results_temp
              }

              # post-foreach end-of-batch snapshot (matches pattern).
              core_update_session_info(ssEnv, save_to_disk = TRUE)

              results_temp <- as.data.frame(results_temp)
              colnames(results_temp) <- result_colnames
              if (exists("results"))
              {
                if (exists("results_temp"))
                  results_temp <- plyr::rbind.fill(results, results_temp)
                else
                  results_temp <- results
                rm(results)
              }
              # append to csv file
              if(!file.exists(fileNameResultsTemp))
                utils::write.table(x =  results_temp, file = fileNameResultsTemp, append = FALSE, row.names = FALSE, sep="\t")
              else
                utils::write.table(x =  results_temp, file = fileNameResultsTemp, append = TRUE, row.names = FALSE, sep="\t", col.names = FALSE)
            }
          }
        }


      if(file.exists(fileNameResultsTemp))
        results <-  utils::read.csv2(file = fileNameResultsTemp, header = TRUE , stringsAsFactors = FALSE, sep="\t")

      # drop column if all is na
      if (!exists("results"))
        next
      results <- results[, colSums(is.na(results)) != nrow(results)]
      results <- results[,result_colnames]
      results$SPECIFICITY <- round(as.numeric(results$SPECIFICITY),4)
      results$SENSITIVITY <- round(as.numeric(results$SENSITIVITY),4)

      results$P_to_be_Case_cond_to_be_Epimutated <- round(as.numeric(results$P_to_be_Case_cond_to_be_Epimutated),4)
      results$P_to_be_Control_cond_to_be_Not_Epimutated <- round(as.numeric(results$P_to_be_Control_cond_to_be_Not_Epimutated),4)
      results$JSD <- round(as.numeric(results$JSD),4)

      results$SCORE <- round(results$SPECIFICITY,4) +
        round(results$SENSITIVITY,4) +
        round(results$P_to_be_Control_cond_to_be_Not_Epimutated,4) +
        round(results$P_to_be_Case_cond_to_be_Epimutated,4) +
        round(results$JSD,4)

      results <- results[,result_colnames]

      # sort by SCORE descending
      results <- unique(results[order(-results$SCORE),])
      utils::write.csv2(x = results, file = fileNameResults,row.names = FALSE)
      rm(results)
      if(file.exists(fileNameResultsTemp))
        file.remove(fileNameResultsTemp)

    }
    core_close_env()
  }
