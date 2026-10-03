#' The column-name suffix of a transformation
#'
#' A transformed covariate is a new column, and the name says which
#' transformation made it. This is not decoration. The coefficient of a model is
#' named after the column it was fitted on, so a covariate transformed in place
#' would be reported as \code{AGE_ESTIMATE} while the number is the coefficient
#' of \code{exp(AGE)} - the same defect that had a polynomial estimate carrying a
#' name that was not its own.
#'
#' \code{"scale"} answers \code{"SCALED"} and not \code{"SCALE"}, because that
#' suffix was already in use before this function existed and
#' \code{io_inference_file_name()} and \code{enrich_phenotype_analysis_name()}
#' both strip and re-add it by name.
#'
#' @param transformation One name from the vocabulary of
#'   \code{io_transform_apply()}.
#' @return The suffix, with no leading underscore, or \code{NA_character_} for
#'   \code{"none"} - which makes no column and leaves the name alone.
#'
#' @keywords internal
#' @noRd
io_transform_suffix <- function(transformation) {

  transformation <- trimws(as.character(transformation))
  if (!length(transformation) || all(is.na(transformation)) ||
      !any(nzchar(transformation)) || identical(transformation[1L], "none"))
    return(NA_character_)

  transformation <- transformation[1L]
  if (identical(transformation, "scale")) return("SCALED")

  # pow_2 -> POW2, quantile_3 -> QUANTILE3, log10 -> LOG10. The underscore goes
  # because a covariate name can carry one of its own and the suffix has to be
  # recognisable at the end of the result.
  toupper(gsub("_", "", transformation))
}
