#' Validate the transformations of each request against the vocabulary
#'
#' The door check for `transformation_y`, `transformation_x` and
#' `covariates_transformation`, run beside [assoc_validate_scope()] and
#' [assoc_validate_aggregation()] in `association_analysis()` and for the same
#' reason: a request that cannot be honoured should stop before anything is
#' computed, rather than be half-answered.
#'
#' @section Why the door and not the point of use:
#' Nothing validated the vocabulary anywhere, because until the vocabulary had
#' one place there was nothing to validate against. `switch()` fell through to
#' its default and returned the values untouched, so a request naming `"lgo10"`
#' was computed on untransformed data; what saved it from being invisible was a
#' guard elsewhere that renames the recorded transformation to `NA_<name>` when
#' the values come back unchanged - a weak signal, in a column of a result file,
#' and only after the whole analysis has run.
#'
#' Refusing inside [io_transform_apply()] is not enough on its own either: the
#' calls that transform the dependent and independent variables sit inside
#' `try()`, which swallows the error and leaves exactly the old behaviour. The
#' refusal has to happen before that, which means here.
#'
#' @section What is checked:
#' Every value is a name [io_transform_known()] recognises, which includes the
#' parameterised families `pow_<n>` and `quantile_<n>`. `"factor"` is refused for
#' a covariate and the message says what to use instead: `covariates_dummy`
#' encodes a categorical covariate, which relabelling it as a factor does not.
#' The length of `covariates_transformation` against `covariates` is checked
#' where the pairing happens, in [assoc_covariates_model()], because that is
#' where both have been split.
#'
#' @param inference_details The request table.
#' @return `inference_details`, unchanged, or the call stops.
#'
#' @keywords internal
#' @noRd
assoc_validate_transformation <- function(inference_details) {

  if (is.null(inference_details) || nrow(inference_details) == 0)
    return(inference_details)

  legal <- paste0(paste(io_transform_vocabulary(), collapse = ", "),
                  ", pow_<n>, quantile_<n>")

  refuse <- function(z, field, value, extra = "") {
    stop("inference_details row ", z, ": '", field, "' is \"", value,
         "\", which is not a transformation this package knows.\n",
         "  Known: ", legal, "\n", extra,
         "  It used to be accepted and ignored, so the analysis ran on ",
         "untransformed data under the name of a transformation.",
         call. = FALSE)
  }

  for (z in seq_len(nrow(inference_details))) {

    for (field in c("transformation_y", "transformation_x")) {
      value <- inference_details[[field]]
      if (is.null(value)) next
      value <- trimws(as.character(value[z]))
      if (!io_transform_known(value)) refuse(z, field, value)
    }

    # transformation_x = "factor" asks for a categorical independent variable.
    # The two-group and k-group families already take it as one; the regression
    # families fit one slope and report one p-value per coefficient, with no
    # overall test for a k-level factor, so they would answer a question the
    # request did not ask. It used to be accepted there and ignored.
    tx <- inference_details[["transformation_x"]]
    family <- inference_details[["family_test"]]
    if (!is.null(tx) && !is.null(family) &&
        identical(trimws(as.character(tx[z])), "factor") &&
        !isTRUE(tryCatch(assoc_is_family_dicotomic(as.character(family[z])),
                         error = function(e) FALSE)))
      stop("inference_details row ", z, ": transformation_x is \"factor\" and ",
           "family_test is \"", family[z], "\", which fits the independent ",
           "variable as a number.\n",
           "  A categorical independent variable is taken as such by the ",
           "two-group and k-group families (t.test, wilcoxon, kruskal.test, ",
           "binomial, ...); for a regression family declare its order with ",
           "independent_variable_order, or use one of those families.\n",
           "  It used to be accepted and ignored, so the model ran on the ",
           "variable as a number under the name of a factor.", call. = FALSE)

    value <- inference_details[["covariates_transformation"]]
    if (is.null(value)) next
    # Positional against covariates, so duplicates are kept: "exp + exp" is two
    # entries and dropping one would shift every entry after it.
    each <- util_split_and_clean(value[z], unique_values = FALSE)
    for (one in each) {
      if (identical(trimws(one), "factor"))
        refuse(z, "covariates_transformation", one,
               paste0("  \"factor\" is not one of them here: a categorical ",
                      "covariate goes in covariates_dummy, which encodes it, ",
                      "rather than being relabelled.\n"))
      if (!io_transform_known(one)) refuse(z, "covariates_transformation", one)
    }
  }

  inference_details
}
