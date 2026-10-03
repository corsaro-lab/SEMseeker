#' The plain names of the transformation vocabulary
#'
#' The parameterised forms, \code{pow_<n>} and \code{quantile_<n>}, are not in
#' this list because they are families and not names: carrying their argument is
#' the point. \code{io_transform_known()} answers for all of them, and both the
#' door check and \code{io_transform_apply()} go through it so the two cannot
#' come to disagree about what is legal.
#'
#' @return Character vector of names.
#' @keywords internal
#' @noRd
io_transform_vocabulary <- function()
  c("none", "scale", "log", "log2", "log10", "exp", "factor")

#' Is this a transformation the package knows?
#'
#' @param transformation One name, possibly parameterised.
#' @return \code{TRUE} or \code{FALSE}. An empty string and \code{NA} are
#'   \code{TRUE}: they mean \code{"none"}.
#' @keywords internal
#' @noRd
io_transform_known <- function(transformation) {

  transformation <- trimws(as.character(transformation))
  if (!length(transformation) || all(is.na(transformation)) ||
      !any(nzchar(transformation)))
    return(TRUE)

  transformation <- transformation[1L]
  if (transformation %in% io_transform_vocabulary()) return(TRUE)

  # pow_<n>: any number, fractions included, so pow_0.5 is the square root.
  if (grepl("^pow_-?[0-9]+(\\.[0-9]+)?$", transformation)) return(TRUE)
  # quantile_<n>: an integer of two or more, because one tile is not a tiling.
  if (grepl("^quantile_[0-9]+$", transformation))
    return(as.numeric(sub("^quantile_", "", transformation)) >= 2)

  FALSE
}

#' Apply one declared transformation to values
#'
#' The vocabulary of transformations, in one place. It was written twice, once
#' for the dependent variable and once for the independent one, with
#' \code{quantile_<n>} handled beside the second copy rather than inside it, so
#' the two lists had already drifted and a third copy for the covariates would
#' have guaranteed they stayed apart.
#'
#' @section The vocabulary:
#' \describe{
#'   \item{\code{"none"}}{The values, unchanged. An empty string and \code{NA}
#'     mean this.}
#'   \item{\code{"scale"}}{Centred and scaled to unit variance.}
#'   \item{\code{"log"}, \code{"log2"}, \code{"log10"}}{The logarithm in that
#'     base.}
#'   \item{\code{"exp"}}{The exponential.}
#'   \item{\code{"pow_<n>"}}{Raised to the power \code{n}, e.g. \code{"pow_2"}
#'     for the square. \code{n} may be fractional: \code{"pow_0.5"} is the
#'     square root.}
#'   \item{\code{"quantile_<n>"}}{Replaced by the index of its \code{n}-tile. A
#'     column with fewer than \code{n} distinct values becomes zeros, which is
#'     the behaviour this had where it was written before and is kept so that a
#'     request that used it reads the same.}
#'   \item{\code{"factor"}}{Coerced to a factor. Admissible for the dependent
#'     and independent variables; for a covariate the request says
#'     \code{covariates_dummy} instead, which encodes it properly rather than
#'     relabelling it.}
#' }
#'
#' @section An unknown name is refused:
#' It used to be a silent no-op. \code{switch()} fell through to its default and
#' returned the values untouched, while the request's own string was recorded in
#' the result as \code{TRANSFORMATION_Y} or \code{TRANSFORMATION_X}: a result
#' that says \code{"lgo10"} over data that was never transformed. Nothing
#' validated the vocabulary anywhere, because until now there was no one place
#' that held it.
#'
#' @param values A vector or a data.frame. A data.frame is transformed
#'   column-wise, which is what the dependent variable needs: it arrives as one
#'   column per sample.
#' @param transformation One name from the vocabulary above.
#' @return The transformed values, in the shape they arrived in.
#'
#' @keywords internal
#' @noRd
io_transform_apply <- function(values, transformation) {

  transformation <- trimws(as.character(transformation))
  if (!length(transformation) || all(is.na(transformation)) ||
      !any(nzchar(transformation)))
    transformation <- "none"
  transformation <- transformation[1L]

  if (identical(transformation, "none"))
    return(values)

  # The two parameterised names carry their argument, the way the family strings
  # do (polynomial_2_1, limma_2). Parsed before the switch, so the switch stays a
  # list of plain names.
  if (grepl("^pow_", transformation)) {
    n <- suppressWarnings(as.numeric(sub("^pow_", "", transformation)))
    if (is.na(n))
      stop("transformation \"", transformation, "\": pow_<n> needs a number, ",
           "e.g. pow_2 for the square or pow_0.5 for the square root.",
           call. = FALSE)
    return(values ^ n)
  }

  if (grepl("^quantile_", transformation)) {
    qq <- suppressWarnings(as.numeric(sub("^quantile_", "", transformation)))
    if (is.na(qq) || qq < 2)
      stop("transformation \"", transformation, "\": quantile_<n> needs an ",
           "integer of 2 or more, e.g. quantile_3 for tertiles.", call. = FALSE)
    ntile_one <- function(x) {
      if (length(unique(x)) >= qq) as.numeric(dplyr::ntile(x, n = qq))
      else rep(0, length(x))
    }
    if (is.data.frame(values) || is.matrix(values))
      return(as.data.frame(apply(values, 2, ntile_one)))
    return(ntile_one(values))
  }

  switch(
    transformation,
    "scale"  = scale(values),
    "log"    = log(values),
    "log2"   = log2(values),
    "log10"  = log10(values),
    "exp"    = exp(values),
    "factor" = if (is.data.frame(values)) as.data.frame(lapply(values, as.factor))
               else as.factor(values),
    stop("transformation \"", transformation, "\" is not one this package knows.\n",
         "  Known: ", paste(io_transform_vocabulary(), collapse = ", "),
         ", pow_<n>, quantile_<n>.\n",
         "  It used to be accepted and ignored, which left the result naming a ",
         "transformation the data had never had.", call. = FALSE)
  )
}
