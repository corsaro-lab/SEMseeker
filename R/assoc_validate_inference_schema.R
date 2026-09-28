#' Validate the schema of inference_details (strict, with helpful errors)
#'
#' Two responsibilities:
#'
#' 1. **Fill missing optional columns** with NA, so downstream code can
#'    always reference them without `is.null()` checks.
#' 2. **Reject unknown columns** with a clear diagnostic — including a
#'    fuzzy-match suggestion when the unknown name is close to an
#'    expected one (typical case: typo in setup like `phenolyser` vs
#'    `phenolyzer`, or `areas_sql_condtion` vs `areas_sql_condition`).
#'
#' The vocabulary of legal columns is the source of truth for what the
#' package accepts from the user. Any new user-input column must be
#' registered here.
#'
#' @param inference_details data.frame from the user.
#' @param strict Whether to error on unknown columns (`TRUE`, default,
#'   recommended) or just warn (`FALSE`, lenient mode for back-compat).
#' @return data.frame with all expected columns present (NA where missing).
#' @keywords internal
assoc_validate_inference_schema <- function(inference_details, strict = TRUE) {
  inference_details <- as.data.frame(inference_details)

  # Canonical user-input vocabulary for inference_details.
  # Source of truth: this list. New user-input fields go HERE.
  # (Internal runtime-stamped fields like node_name, session_id, start_time,
  # processed_items are NOT in this list — they are added by the engine
  # after the row has been validated.)
  expected_columns <- c(
    "independent_variable",
    "family_test",
    "covariates",
    "covariates_dummy",
    "covariates_pca",
    "collinearity_check",
    "transformation_y",
    "transformation_x",
    # AI-308: which of the two aggregation branches the request wants. SAMPLE
    # reduces the positions to one number per sample, INSTANCE to one number per
    # instance of the region class. They are alternatives, not addends, so the
    # request has to name one. The region classes are NOT named here: they are
    # the (AREA, SUBAREA) pairs of the run, declared with
    # association_analysis(areas =, subareas =) and built at runtime.
    "scope",
    # AI-248: which aggregation of the feature is tested (SUM, MEAN, MEDIAN,
    # VARIANCE, IQR, MODELOW, MODEHIGH). Mandatory: with several aggregations
    # over the same scope, a request that does not name one does not identify
    # what it wants tested.
    "aggregation",
    "filter_p_value",
    "samples_sql_condition",
    "areas_sql_condition",
    "association_results_sql_condition"
  )

  # Engine-stamped columns that may appear when validate is called on a
  # already-run inference (e.g. for re-validation in a resume). Allowed
  # without flagging as unknown.
  runtime_stamped <- c("node_name", "session_id", "start_time", "end_time",
                       "processed_items", "processed_time")

  allowed_columns <- c(expected_columns, runtime_stamped)

  # ---- (a) Reject unknown columns with helpful diagnostic ----
  unknown_cols <- setdiff(colnames(inference_details), allowed_columns)

  # AI-255: a column that was REMOVED is not a typo, and telling the user to
  # register it as legal would send them the wrong way. Name it, say it is gone,
  # and say what replaced it.
  retired <- c(
    depth_analysis = paste0(
      "removed. The granularity of an artefact is said by SCOPE, AREA and ",
      "SUBAREA, which say more: those pairs are only partially ordered, so an ",
      "integer scale projected a lattice onto a line. A request that used ",
      "depth_analysis = 1 now writes scope = \"SAMPLE\"; anything above 1 ",
      "wrote scope = \"INSTANCE\", and the region classes it ranges over are ",
      "the (AREA, SUBAREA) pairs of the run"),
    scopes = paste0(
      "removed. It named the region classes a second time, and they are ",
      "already the (AREA, SUBAREA) pairs of the run: declare them with ",
      "association_analysis(areas =, subareas =) and they are built at ",
      "runtime. What the request names now is `scope`: SAMPLE or INSTANCE, ",
      "which of the two aggregation branches to run. The two are mutually ",
      "exclusive: a request that wants both writes two rows"))
  hit <- intersect(names(retired), unknown_cols)
  if (length(hit) > 0)
    stop("inference_details carries column(s) that no longer exist:\n  ",
         paste(sprintf("'%s' — %s", hit, retired[hit]), collapse = "\n  "),
         call. = FALSE)

  if (length(unknown_cols) > 0) {
    suggestions <- vapply(unknown_cols, function(uk) {
      d <- utils::adist(uk, expected_columns, ignore.case = TRUE)[1, ]
      close_idx <- which(d <= max(2, nchar(uk) %/% 4))
      if (length(close_idx) > 0) {
        best <- expected_columns[close_idx[which.min(d[close_idx])]]
        sprintf("'%s' (did you mean '%s'?)", uk, best)
      } else {
        sprintf("'%s'", uk)
      }
    }, character(1))

    msg <- sprintf(
      paste0(
        "inference_details has unknown column(s):\n  %s\n",
        "Expected columns: %s\n",
        "Add new fields to assoc_validate_inference_schema()$expected_columns ",
        "if they should be legal."
      ),
      paste(suggestions, collapse = "\n  "),
      paste(expected_columns, collapse = ", ")
    )
    if (strict) {
      stop(msg, call. = FALSE)
    } else {
      warning(msg, call. = FALSE)
    }
  }

  # ---- (b) Fill missing optional columns with NA ----
  missing_cols <- setdiff(expected_columns, colnames(inference_details))
  for (ev in missing_cols) {
    inference_details[, ev] <- NA
  }

  inference_details
}
