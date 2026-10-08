#--- util_absent_as_zero ---------------------------------------------------------
# The one rule for an empty cell in a pivot that reaches a model.
#
# For every marker SEMseeker counts or derives from a deviation (MUTATIONS,
# LESIONS, DELTAS, DELTAR and the quantised DELTAP, DELTAQ, DELTARP, DELTARQ) an
# empty cell means that nothing was found there: it is a zero. The derived
# pivots are empty exactly where DELTAS is 0, and the per-sample path writes a
# position only where some sample has a value.
#
# For SIGNAL an empty cell is a missing measurement, and a beta of 0 is an
# extreme value, not an absence: the gap is refused, never filled.
#
# Only pivot values pass through here. Sample-sheet columns (the independent
# variable, the covariates) keep their NA: a missing age is not an age of 0.
#
# @param values matrix or data.frame of pivot values only, no key columns.
# @param marker the marker the values belong to.
# @return values with NA replaced by 0.
# @keywords internal
util_absent_as_zero <- function(values, marker) {
  if (!anyNA(values))
    return(values)
  if (identical(toupper(as.character(marker)), "SIGNAL"))
    stop("ERROR: I'm stopping here, the SIGNAL pivot has missing values: a ",
         "missing measurement cannot be read as 0. Impute the signal (inpute) ",
         "or remove the missing positions before the association.",
         call. = FALSE)
  values[is.na(values)] <- 0
  values
}
