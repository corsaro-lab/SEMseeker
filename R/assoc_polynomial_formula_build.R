#' Build a polynomial regression formula, adjusting for covariates
#'
#' Pure helper: no I/O, no side effects. Constructs a formula of the form
#' \code{y ~ I(x^1) + I(x^2) + ... + I(x^degree) + cov1 + cov2 + ...}
#'
#' @section Why the covariates are main effects:
#' Adjusting for a covariate means holding it constant, and the term that does
#' that is the covariate on its own. The polynomial coefficients are then the
#' shape of the curve at a fixed covariate value, which is the quantity the
#' request is asking for.
#'
#' Through 0.99.5 the covariates entered only as interactions with each
#' polynomial term, \code{I(x^k):cov}, with no \code{cov} term of its own. That is not a
#' different amount of adjustment, it is none: writing the model out,
#'
#' \preformatted{
#'   was      y = b0 + b1 x + b2 x^2 + c1 (x cov) + c2 (x^2 cov)
#'   is       y = b0 + b1 x + b2 x^2 + a cov
#' }
#'
#' the covariate is at the first power in both, so neither model is more linear
#' in it than the other. What the old one lacked was a term able to shift the
#' level: its derivative with respect to the covariate was \code{c1 x + c2 x^2},
#' which is \strong{zero at x = 0} and grows with x. Nobody chose that
#' constraint; it is what is left when the main effect is missing. And with no
#' other place for the covariate to enter, it entered the slopes, so the
#' polynomial coefficients came out as the curve at \code{cov = 0} - age zero,
#' reference sex - an extrapolation outside the data presented as an adjusted
#' estimate.
#'
#' @section What is assumed, and what is not:
#' For a two-valued covariate, a dummy among them, nothing is assumed: one
#' coefficient on a 0/1 column is the whole difference between the groups.
#'
#' For a continuous covariate it is an assumption, and it is worth naming: the
#' covariate's effect on the response is taken to be a straight line. A caller
#' who doubts it does not need new code, because \code{covariates} is a list of
#' column names: supplying both \code{AGE} and a derived \code{AGE_SQ} column
#' gives a quadratic in age. That was unreachable while the terms were
#' interactions only.
#'
#' The shape of the curve in x is assumed to be the same at every covariate
#' value. Letting it vary is a different question - whether the covariate
#' modifies the effect rather than confounds it - and it is deliberately not
#' reachable from here.
#'
#' The asymmetry is deliberate: flexible in the variable under study, degree 2
#' or 3, and rigid in the covariates, degree 1. The covariates are nuisance to
#' be held constant, not the object of the request.
#'
#' @param dependent_variable  Name of the response column (string).
#' @param independent_variable Name of the predictor column (string).
#' @param degree              Polynomial degree (positive integer).
#' @param covariates          Character vector of covariate column names. Empty
#'   strings are dropped: \code{assoc_sig_formula_vars()} can return one, and it
#'   used to reach the formula as a bare \code{I(x^k):} and stop the run on a
#'   parse error.
#' @return A \code{formula} object.
#'
assoc_polynomial_formula_build <- function(dependent_variable, independent_variable, degree, covariates) {
  polynomial_terms <- paste0("I(", independent_variable, "^", seq_len(degree), ")")
  x_part <- paste(polynomial_terms, collapse = " + ")

  covariates <- covariates[nzchar(as.character(covariates))]

  if (length(covariates) > 0)
    formula_string <- paste(dependent_variable, "~", x_part, "+",
                            paste(covariates, collapse = " + "))
  else
    formula_string <- paste(dependent_variable, "~", x_part)

  stats::as.formula(formula_string)
}
