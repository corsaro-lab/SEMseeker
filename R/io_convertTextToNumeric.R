io_convertTextToNumeric <- function(text) {

  # Get the system's locale settings
  locale_settings <- Sys.localeconv()

  # Extract the decimal separator from the locale settings
  decimal_separator <- paste("[",locale_settings["decimal_point"],"]", sep="")

  # Example usage:
  # core_log_event(paste("Decimal separator:", decimal_separator))

  if(length(gregexpr("\\.", text)[[1]])>1 | length(gregexpr("\\,", text)[[1]])>1)
  {
    # Check if the text contains a decimal separator
    if (grepl("[.,]", text)) {

      # Multiple dots → European format: dots are thousands sep, comma is decimal
      # e.g. "1.200.000,34" → remove dots → "1200000,34" → replace comma → 1200000.34
      converted_value <- suppressWarnings(as.numeric(gsub(",", ".", gsub("\\.", "", text))))
      if (!is.na(converted_value))
        return(converted_value)

      # Multiple commas → US format: commas are thousands sep, dot is decimal
      # e.g. "1,200,000.34" → remove commas → "1200000.34" → 1200000.34
      converted_value <- suppressWarnings(as.numeric(gsub(",", "", text)))
      if (!is.na(converted_value))
        return(converted_value)
    }
  }
  else
    {
      # Try to convert with "." as the decimal separator
      # and  "," as thousands separator
      converted_value <- suppressWarnings(as.numeric(gsub(",", ".", text)))

      # If conversion successful, return the numeric value
      if (!is.na(converted_value))
        return(converted_value)

      # Try to convert with "," as the decimal separator
      converted_value <- suppressWarnings(as.numeric(gsub("\\.", ",", text)))

      # If conversion successful, return the numeric value
      if (!is.na(converted_value))
        return(converted_value)
    }

  # If no decimal separator found or conversion unsuccessful, return NA
  return(NA)
}



# Examples (see unit tests test-sem_misc.R):
# io_convertTextToNumeric("1,200,000.34") == 1200000.34
# io_convertTextToNumeric("1.200.000,34") == 1200000.34
