#' Reduce a set of values to one number, by aggregation name (internal)
#'
#' AI-248, revised by AI-255. The other half of the AGGREGATION axis:
#' [util_aggregations_allowed()] declares which names are legal for an artefact,
#' this one is the only place that says what each name *computes*. Keeping the
#' two together in one vocabulary is what stops a column called `MEDIAN` from
#' holding a mean.
#'
#' `VALUE` is the identity - the block is already a single position, there is
#' nothing to reduce. It is not a degenerate `SUM`: naming it `SUM` would invite
#' the reader to believe a reduction happened.
#'
#' `MODELOW` / `MODEHIGH` are the only distribution-shape operators: they
#' estimate the highest density peak on each side of 0.5, so they need the whole
#' distribution rather than a streaming reduce, and enough of it - they delegate
#' to [util_signal_descriptors()], which refuses below a minimum numerosity.
#'
#' AI-255 removed the `N_PROBES` branch that used to live here. The number of
#' usable positions is a property of the **imputation** - how many probes of that
#' sample survived the treatment of missing values - not an aggregation of a
#' marker, and it belongs in `SAMPLE_SHEET_RESULT` with the other descriptive
#' properties of the sample. Note that nothing is lost for the density: the
#' `MEAN` of a binary marker *is* the density, denominator included.
#'
#' @param values numeric vector. Non-finite values are dropped.
#' @param aggregation one name from [util_aggregation_vocabulary()].
#' @return a single numeric, `NA_real_` when the values do not support it.
#' @keywords internal
#' @noRd
util_aggregate_values <- function(values, aggregation) {

  aggregation <- toupper(as.character(aggregation))
  values <- as.numeric(values)
  values <- values[is.finite(values)]

  if (length(values) == 0L)
    return(NA_real_)

  if (identical(aggregation, "VALUE")) {
    if (length(values) > 1L)
      stop("util_aggregate_values(): aggregation 'VALUE' means the block holds ",
           "a single position, but ", length(values), " values were passed. ",
           "Name the reduction that produced the one number.", call. = FALSE)
    return(values)
  }

  switch(
    aggregation,
    SUM       = sum(values),
    MEAN      = mean(values),
    MEDIAN    = stats::median(values),
    VARIANCE  = stats::var(values),
    IQR       = stats::IQR(values),
    MODELOW   = util_signal_descriptors(values, beta = TRUE)$MODELOW,
    MODEHIGH  = util_signal_descriptors(values, beta = TRUE)$MODEHIGH,
    stop("util_aggregate_values(): unknown aggregation '", aggregation,
         "'. Legal names: ",
         paste(util_aggregation_vocabulary(), collapse = ", "), ".")
  )
}
