#' Validate the scope of each request against the run
#'
#' AI-308. `SCOPE` is the coordinate that says over what extent one number is
#' valid: `SAMPLE` reduces the positions of a region class to one number per
#' sample, `INSTANCE` to one number per instance of that class: one per gene,
#' island, cytoband or probe. They are **two branches of aggregation**, and
#' [io_pivot_build()] has always kept them apart; what was missing upstream was
#' anyone choosing between them, so every request produced both and the two were
#' silently summed into one result file.
#'
#' This is the door check, run beside [assoc_validate_aggregation()] in
#' `association_analysis()` and for the same reason: a request that cannot be
#' honoured must be refused before any result exists, not discovered after part
#' of the output has been written.
#'
#' **Why `scope` is required and has no default.** It is the same argument that
#' made `aggregation` required in AI-248. Any default would silently give a
#' request half of what the previous release gave it: an inference CSV that
#' looks complete and simply never tested one of the two branches. A missing
#' coordinate is a request that does not identify what it wants.
#'
#' **Why `SAMPLE` is refused for the batch families.** `limma_`/`voom_` estimate
#' a prior variance across the instances they are handed; at `SCOPE = SAMPLE`
#' the artefact holds one row, the fit degenerates to OLS and contaminates that
#' pool. While the two branches were additive this could be an `INFO` that
#' skipped the collapsed keys and carried on with the per-instance ones; with
#' one branch per request there is nothing left to carry on with, so it is an
#' error rather than an empty file.
#'
#' @param inference_details validated data.frame of requests.
#' @return `inference_details` with `scope` cleaned to its canonical spelling.
#'   Callers must use the returned value.
#' @keywords internal
#' @noRd
assoc_validate_scope <- function(inference_details) {

  if (is.null(inference_details) || nrow(inference_details) == 0)
    return(inference_details)

  legal <- io_scope_vocabulary()

  for (z in seq_len(nrow(inference_details))) {
    requested <- inference_details$scope[z]

    if (is.null(requested) || length(requested) == 0 || all(is.na(requested)) ||
        !any(nzchar(as.character(requested))))
      stop("inference_details row ", z, ": 'scope' is required. Name which of ",
           "the two aggregations to run: \"SAMPLE\" (one number per sample, ",
           "over the positions of each region class of the run) or ",
           "\"INSTANCE\" (one number per instance of the class: per gene, ",
           "per island, per probe). They are alternatives, not addends: a ",
           "request that wants both writes two rows.", call. = FALSE)

    # io_scope_validate() owns the vocabulary and the message for a value that
    # is not in it; prefixing the row number is all this adds.
    scope <- tryCatch(io_scope_validate(core_name_cleaning(as.character(requested)[1])),
                      error = function(e)
                        stop("inference_details row ", z, ": ",
                             conditionMessage(e), call. = FALSE))

    family_test <- util_split_and_clean(inference_details$family_test[z])
    if (identical(scope, "SAMPLE") && any(grepl("^(limma|voom)_", family_test)))
      stop("inference_details row ", z, ": family_test = '",
           paste(family_test, collapse = ", "), "' cannot be fitted at ",
           "scope = \"SAMPLE\". Those families estimate a prior variance ",
           "across the instances they are handed, and a collapsed artefact ",
           "holds one row: the fit degenerates to OLS and contaminates the ",
           "eBayes prior. Use scope = \"INSTANCE\" for a batch family, or a ",
           "per-sample family test for the collapsed burden.", call. = FALSE)

    inference_details$scope[z] <- scope
  }

  inference_details
}
