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
