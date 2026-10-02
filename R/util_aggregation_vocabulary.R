#' Every aggregation name the taxonomy knows (internal)
#'
#' AI-248, extended by AI-255. The full vocabulary of the AGGREGATION axis,
#' independent of any artefact: every name that is a name at all.
#'
#' This answers a different question from [util_aggregations_allowed()], and the
#' two must not be confused. This one is a **spell-check**: is the string the
#' caller wrote the name of an aggregation, or a typo? The other one is the
#' **semantic** check: is this operator meaningful for *this* artefact. A
#' validator needs both, in that order - a helpful error for `"MEDAIN"` ("not a
#' name we know") is not the same as the error for `MODELOW` per instance ("a
#' real name, wrong scope").
#'
#' It exists as its own function because the callers used to obtain the union by
#' asking `util_aggregations_allowed("SIGNAL", "BETA", ...)` - the one pair that
#' happens to admit everything. That worked, but it read as an assumption about
#' the scale of the session, which it never was.
#'
#' @return character vector of aggregation names.
#' @keywords internal
#' @noRd
util_aggregation_vocabulary <- function() {
  # The union over both classes, plus the identity. SIGNAL/BETA at SAMPLE scope
  # is what widens the set - it is the only combination carrying the two modes -
  # but here it names the axis, not the run: the vocabulary of a language does
  # not depend on the sentence.
  c(unique(unlist(lapply(c(TRUE, FALSE), function(discrete)
      util_aggregations_allowed("SIGNAL", "BETA", discrete = discrete,
                                default = FALSE, scope = "SAMPLE")))),
    "VALUE")
}
