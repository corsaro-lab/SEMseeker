#' Summarise an association result over its metric columns
#'
#' Reduces a result table to the metrics it carries, which means deciding what
#' counts as a metric. That decision is \code{metrics_properties}, a dataset
#' shipped with the package: it names the metrics and says, for each, whether
#' scaling affects it and whether higher is better. Columns not in that registry
#' are dropped.
#'
#' @section Two kinds of column are carried past the registry:
#' A p-value is not a metric to rank and there is one per adjustment family and
#' per method, so registering them is neither possible nor useful; they are
#' matched by name instead. The bounds of a confidence interval are the same kind
#' of thing - an attribute of a metric rather than a metric - so they travel the
#' same way, matched as \code{_CI_LOWER} and \code{_CI_UPPER}.
#'
#' Without that second exception the AUC would appear in a summary while its
#' interval vanished from it, silently, because \code{C_STATISTIC_AUC} is
#' registered and its two bounds are not.
#'
#' @section A limitation of the registry, not of this function:
#' \code{Higher_the_Better} is a boolean, which presupposes that goodness is
#' monotone in the value. For an AUC it is not - 0.1 is as informative as 0.9,
#' and what is monotone is the distance from 0.5 - and the same holds for a
#' correlation. No value of the flag is right for those, so a ranking over them
#' is wrong whichever way it is set. Fixing it means changing the schema of a
#' bundled dataset and the behaviour of \code{sem_metrics_ranking()}, so it is
#' tracked rather than done here.
#'
#' @param inference_details The request rows whose results are summarised. The
#'   result table itself is read by \code{assoc_data_extractor()}, not passed in.
#' @param destination_folder Where the summary is written.
#' @param result_folder The results to read. When given, the session is reopened
#'   on it with \code{start_fresh = FALSE}; when empty the current session is
#'   used. A reopen applies the defaults for whatever the call does not name, so
#'   a summary of a folder whose run used non-default coordinates has to name
#'   them.
#' @param ... Passed to the session setup and to the extractor.
#' @return The summary table.
#'
#' @keywords internal
#' @noRd
assoc_analysis_summary <- function(inference_details,destination_folder="", result_folder="", ...)
{
  if(result_folder!="")
    ssEnv <- core_init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  else
    ssEnv <- core_get_session_info()
  association_data <- assoc_data_extractor(inference_details, destination_folder, result_folder, ...)


  sem_available_metrics <- toupper(SEMseeker::metrics_properties[,"Metric"])

  # The registry says which metrics are comparable and in which direction, so a
  # p-value column is carried past it rather than registered: there is one per
  # adjustment family and per method, and they are not metrics to rank. The
  # bounds of a confidence interval are the same kind of thing - an attribute of
  # a metric, not a metric - so they travel the same way. Without this the AUC
  # would appear in a summary and its interval would vanish from it, silently.
  carried <- grepl("PVALUE|_CI_LOWER$|_CI_UPPER$", colnames(association_data))
  if(any(carried))
    sem_available_metrics <- c(sem_available_metrics, colnames(association_data)[carried])

  # remove not numeric columns metrics
  metrics_to_remove <- colnames(association_data[,!vapply(association_data, is.numeric, logical(1))])
  sem_available_metrics <- sem_available_metrics[!(sem_available_metrics %in% metrics_to_remove)]

  sem_available_metrics <- sem_available_metrics[sem_available_metrics %in% colnames(association_data)]

  #  create a summary table for the association analysis grouping by AREA,SUBAREA,MARKER,FIGURE and SAMPLES_SQL_CONDITION if exists
  if(any("SAMPLES_SQL_CONDITION" %in% colnames(association_data))) {
    summary_table <- association_data %>%
      dplyr::group_by(.data$AREA, .data$SUBAREA, .data$MARKER, .data$FIGURE, .data$SAMPLES_SQL_CONDITION) %>%
      dplyr::summarise(dplyr::across(sem_available_metrics, list(
        max= ~max(., na.rm=TRUE),
        min= ~min(., na.rm=TRUE),
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE),
        count_below_0.05 = ~sum(. < 0.05, na.rm = TRUE))))
  } else {
    summary_table <- association_data %>%
      dplyr::group_by(.data$AREA, .data$SUBAREA, .data$MARKER, .data$FIGURE) %>%
      dplyr::summarise(dplyr::across(sem_available_metrics, list(
        max= ~max(., na.rm=TRUE),
        min= ~min(., na.rm=TRUE),
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE),
        count_below_0.05 = ~sum(. < 0.05, na.rm = TRUE))))
  }

  summary_table <- as.data.frame(summary_table)
  return(summary_table)
}
