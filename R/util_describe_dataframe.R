#' Describe a data frame with summary statistics per column
#'
#' Pure helper: no I/O, no side effects. Returns a data frame with one row per
#' column of \code{df} and the following fields: Variable, Class, Missing_Values,
#' Missing_Values_Percent, Unique_Values, Mean, Median, Min, Max.
#'
#' @param df A data frame to describe.
#' @return A data frame with summary statistics.
#'
util_describe_dataframe <- function(df) {
  data.frame(
    Variable              = names(df),
    Class                 = vapply(df, function(x) class(x)[1], character(1)),
    Missing_Values        = vapply(df, function(x) sum(is.na(x)), numeric(1)),
    Missing_Values_Percent = round(vapply(df, function(x) sum(is.na(x)) / length(x) * 100, numeric(1)), 2),
    Unique_Values         = vapply(df, function(x) length(unique(x)), numeric(1)),
    Mean   = round(vapply(df, function(x) if (is.numeric(x)) as.numeric(mean(x,   na.rm = TRUE)) else NA_real_, numeric(1)), 2),
    Median = round(vapply(df, function(x) if (is.numeric(x)) as.numeric(stats::median(x, na.rm = TRUE)) else NA_real_, numeric(1)), 2),
    Min    = round(vapply(df, function(x) if (is.numeric(x)) as.numeric(min(x,    na.rm = TRUE)) else NA_real_, numeric(1)), 2),
    Max    = round(vapply(df, function(x) if (is.numeric(x)) as.numeric(max(x,    na.rm = TRUE)) else NA_real_, numeric(1)), 2)
  )
}
