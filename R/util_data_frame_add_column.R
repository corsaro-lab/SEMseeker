#' Add or update a column in a data frame
#'
#' If \code{col_name} already exists in \code{df}, its values are overwritten.
#' If it does not exist, a new column is appended. Handles empty data frames
#' (zero rows) by returning the new column as a single-column data frame.
#'
#' @param df A \code{data.frame}.
#' @param col_name Character scalar: name of the column to add or update.
#' @param value Vector of values to assign to the column.  Must be compatible
#'   with the number of rows in \code{df} (or any length when \code{df} is
#'   empty).
#'
#' @return The modified \code{data.frame} with the column added or updated.
#'
#' @examples
#' df <- data.frame(a = 1:3)
#' SEMseeker:::util_data_frame_add_column(df, "b", c(4, 5, 6))
#'
util_data_frame_add_column <- function(df,col_name, value)
{
  if(any(col_name %in% colnames(df)))
    df[,col_name] <- value
  else
  {
    tmp <- stats::setNames(data.frame(value, stringsAsFactors = FALSE), col_name)
    if (nrow(df)==0)
      df <- tmp
    else
      df <- cbind(df,tmp)
  }
  return(df)
}
