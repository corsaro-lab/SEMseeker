#' Reduce one sample's signal vector to its descriptors
#'
#' Generic descriptive statistics on a numeric vector: no SEM domain knowledge,
#' no I/O. AI-223 calls it once per sample on the whole probe set (scope
#' `SAMPLE`); the region scopes reuse it unchanged on the subset of positions of
#' an area/subarea.
#'
#' The returned keys are built from [io_signal_stats()] rather than written out
#' here, so the computed statistics and the vocabulary declared to the consumer
#' cannot diverge. A statistic this function does not fill stays `NA`.
#'
#' Split at 0.5 for the two beta modes: methylation beta values are bimodal with
#' peaks near 0 (unmethylated) and near 1 (methylated), so the highest density
#' peak on each side of 0.5 is the natural estimate of each mode. On M-values the
#' distribution is unimodal and the split is meaningless, which is why
#' [io_signal_stats()] does not declare the mode columns there.
#'
#' **The numerosity guard (AI-255).** The old guard asked only for two distinct
#' points, which is what `density()` needs to run - not what the estimate needs
#' to mean something. On a handful of values the kernel bandwidth dominates and
#' "the highest peak below 0.5" is noise, but the function would still return a
#' number, and a plausible-looking number is worse than a refusal because nobody
#' goes looking for it. The taxonomy already keeps the modes at `SCOPE = SAMPLE`
#' for this reason (see [util_aggregations_allowed()]); this guard is the second
#' line, for any caller that reaches here with a small vector.
#'
#' @param values numeric vector, already filtered to finite values.
#' @param beta logical. TRUE when the signal is on the [0,1] beta scale.
#' @return named list of statistics, one element per [io_signal_stats()] entry.
#' @keywords internal
#' @noRd
util_signal_descriptors <- function(values, beta = TRUE) {

  wanted <- io_signal_stats(beta)
  out    <- stats::setNames(vector("list", length(wanted)), wanted)
  out[]  <- NA_real_

  # AI-255: N_PROBES describes the INPUT - how many finite values were reduced -
  # and is reported here for callers that want it. It is deliberately not in
  # io_signal_stats(), which is the vocabulary of columns an artefact carries,
  # because no artefact carries it any more: it is a property of the imputation
  # and it travels on the sample sheet.
  out$N_PROBES <- length(values)

  if (length(values) == 0)
    return(out)

  out$MEDIAN   <- stats::median(values)
  out$MEAN     <- mean(values)
  out$VARIANCE <- stats::var(values)
  out$IQR      <- stats::IQR(values)

  # Off the beta scale the modes are not part of the vocabulary at all, so
  # there is nothing left to estimate.
  if (!all(c("MODELOW", "MODEHIGH") %in% wanted))
    return(out)

  # Not enough of a distribution to have two identifiable peaks: leave NA rather
  # than return the bandwidth's opinion.
  if (length(values) < .MODE_MIN_N || length(unique(values)) < 2L)
    return(out)

  dens <- tryCatch(stats::density(values, from = 0, to = 1),
                   error = function(e) NULL)
  if (is.null(dens))
    return(out)

  low <- dens$x < 0.5
  if (any(low))
    out$MODELOW <- dens$x[low][which.max(dens$y[low])]
  if (any(!low))
    out$MODEHIGH <- dens$x[!low][which.max(dens$y[!low])]

  out
}

#' Minimum numerosity for a two-peak density estimate (internal)
#'
#' AI-255. A deliberate, conservative default rather than a derived quantity:
#' below roughly a hundred values a gaussian-kernel density on [0,1] is shaped
#' more by its bandwidth than by the data, and the two-peak structure the modes
#' claim to find is not identifiable. Callers that legitimately want the modes
#' work at `SCOPE = SAMPLE`, where the vector is the whole probe set of a sample
#' (tens to hundreds of thousands of values) or a region class of it (thousands),
#' so the guard never bites there.
#'
#' @keywords internal
#' @noRd
.MODE_MIN_N <- 100L
