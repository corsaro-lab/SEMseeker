#' Post-hoc power of a two-sample comparison from a standardised effect size
#'
#' \code{pwr::pwr.t2n.test()} takes \code{d}, Cohen's standardised difference,
#' and this is the one place that converts to it, because the two callers in
#' \code{assoc_test_model()} were each passing something else.
#'
#' The wilcoxon branch passed the common-language effect size A, which is bounded
#' in \[0, 1\]: at A = 0.5, the exact null, that is read as \code{d = 0.5}, a
#' medium effect, and the reported power was 0.53 where the truth is the
#' significance level. The t.test branch passed the raw difference of the means,
#' which carries the units of the data: the same difference on the beta scale and
#' on the M-value scale gave 0.053 and 0.171 for an effect whose true power was
#' 1.000. Standardising is precisely what removes that.
#'
#' @param d Cohen's d. Non-finite values are clamped: \code{pwr} saturates at
#'   \code{1} well before \code{|d| = 10}, so an infinite effect - two groups
#'   with no overlap and no within-group variance - answers 1 rather than
#'   failing. \code{NA} answers \code{NA}, because a group of one has no
#'   standardised difference and inventing one would be worse than saying so.
#' @param n1,n2 The two group sizes.
#' @param alpha Significance level.
#' @return The power, or \code{NA_real_}.
#'
#' @keywords internal
#' @noRd
assoc_power_from_d <- function(d, n1, n2, alpha = 0.05) {

  d <- suppressWarnings(as.numeric(d)[1L])
  if (length(d) == 0L || is.na(d)) return(NA_real_)
  if (n1 < 2 || n2 < 2) return(NA_real_)

  if (!is.finite(d)) d <- sign(d) * 10
  d <- max(min(abs(d), 10), 0)

  # pwr.t2n.test() calls a quantile function with `lower=` where the formal is
  # `lower.tail`, three times per call. SEMseeker turns R's partial-match
  # warnings on at load - deliberately, to catch that class of bug in its own
  # dispatchers - so a dependency's internal partial match becomes three
  # warnings in the user's session for every position tested. Muffled by
  # message, not with suppressWarnings(), so a real warning from pwr still
  # reaches the caller.
  out <- tryCatch(
    withCallingHandlers(
      pwr::pwr.t2n.test(d = d, n1 = n1, n2 = n2,
                        sig.level = alpha, power = NULL)$power,
      warning = function(w) {
        if (grepl("partial argument match", conditionMessage(w), fixed = TRUE))
          invokeRestart("muffleWarning")
      }),
    error = function(e) NA_real_)
  as.numeric(out)
}
