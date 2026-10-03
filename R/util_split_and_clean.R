#' Split a \code{"+"}-separated request field into its parts
#'
#' The request fields that carry several values - \code{covariates},
#' \code{covariates_dummy}, \code{covariates_transformation},
#' \code{independent_variable_order} - are single strings joined with
#' \code{"+"}, and this is the one place that takes them apart. Whitespace
#' around each part is dropped, so \code{"AGE + BMI"} and \code{"AGE+BMI"} are
#' the same request.
#'
#' @section Why \code{unique_values} exists:
#' The default drops duplicates, which is right for a SET of names: naming the
#' same covariate twice means nothing.
#'
#' It is wrong for a vector that is paired with another by POSITION. With
#' \code{covariates_transformation = "exp + exp + none"} over three covariates,
#' de-duplicating returns two entries, every entry after the first duplicate
#' shifts by one, and a legitimate request - two covariates transformed the same
#' way - is refused for a length mismatch or, worse, silently misaligned. Pass
#' \code{unique_values = FALSE} for those.
#'
#' This was found by a test, not by reading: the length check on
#' \code{covariates_transformation} rejected \code{"exp + exp"} as one entry
#' against two covariates.
#'
#' @param x The field, a string or something coercible to one.
#' @param split Regular expression of the separator. Default \code{"\\+"}.
#' @param unique_values Drop duplicates. \code{TRUE} by default, which is the
#'   behaviour every caller had before this argument existed. \code{FALSE} for a
#'   positional vector.
#' @return Character vector of the trimmed, non-empty parts.
#'
#' @keywords internal
#' @noRd
util_split_and_clean <- function(x, split = "\\+", unique_values = TRUE) {
  # unique_values = FALSE for a vector that is paired with another by POSITION,
  # where duplicates are normal and dropping one shifts every entry after it:
  # "exp + exp + none" for three covariates is three entries, not two.

  x <- as.character(x)
  # Check if the input is NULL, NA, empty string, or character NA
  if(rlang::is_empty(x))
    return("")

  x <- as.character(x)
  if (grepl(split, x)) {
    # Split the string by the specified delimiter
    parts <- unlist(strsplit(x, split))
  }
  else
    parts <- c(x)

  # Remove leading and trailing whitespace from each part
  cleaned_parts <- trimws(parts)

  # remove any empty strings
  cleaned_parts <- cleaned_parts[cleaned_parts != ""]

  cleaned_parts <- cleaned_parts[!is.na(cleaned_parts)]

  # Return the cleaned parts as a vector
  if (!unique_values)
    return(stats::na.omit(cleaned_parts))
  return(na.omit(unique(cleaned_parts)))
}
