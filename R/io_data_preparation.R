#' Title
#'
#' @param family_test test or regression to apply
#' @param transformation_y transformation_y to apply to data
#' @param tempDataFrame data frame to use for test/regression
#' @param independent_variable regressor
#' @param g_start index of the first burden column in tempDataFrame. Checked,
#'   not trusted: the burden columns are found by name (every column but the
#'   independent variable and the covariates) and the run stops if
#'   g_start:g_end does not point at exactly those.
#' @param g_end index of the last burden column in tempDataFrame
#' @param covariates vector of covariates to be found in the sample sheet
#' @param key named list with AREA, SUBAREA, MARKER and FIGURE identifiers, used
#'   to name the artefact in the log when degenerate columns are dropped
#' @param transformation_x not applied here: the independent variable is
#'   transformed upstream, in assoc_covariates_model(), into a column of its own.
#'   Kept so the callers' signature does not change.
#' @param independent_variable_order optional order of the levels of an ordinal
#'   independent variable, "+"-separated. See util_level_order().
#'   before the fit; "none" leaves it untouched
#'
#' @return A named list: \code{tempDataFrame} (the prepared and optionally
#'   transformed data.frame), \code{independent_variableLevels} (the factor
#'   levels of the independent variable, or \code{NULL} for continuous outcomes)
#'   and \code{burden_columns} (the names of the burden columns that survived,
#'   which is what the callers iterate over).
#'
io_data_preparation <- function(family_test,transformation_y,tempDataFrame, independent_variable, g_start, g_end, covariates, key, transformation_x = "none",
                                independent_variable_order = NULL)
{

  #
  ssEnv <- core_get_session_info()

  transformation_y <- as.character(transformation_y)
  transformation_x <- as.character(transformation_x)
  independentVariableIsFactor <- FALSE
  independent_variableLevels <- NULL
  if (assoc_is_family_dicotomic(family_test))
  {
    # The order comes from util_level_order(), which is also what the other
    # factor conversion in sem_prepare_study_for_analysis() uses: the two have
    # different family lists, so an order imposed here only would disagree with
    # the one imposed there. Declared when the request declares it, mixedsort
    # otherwise - never the plain sort that used to be attempted on the line
    # below, which read "10" as smaller than "2".
    level_order <- util_level_order(tempDataFrame[, independent_variable],
                                    independent_variable_order)

    independent_variableLevels <- NA
    independentVariableIsFactor <- FALSE
    if(length(level_order) > 1)
    {
      tempDataFrame[, independent_variable] <- factor(
        as.character(tempDataFrame[, independent_variable]), levels = level_order)
      tempDataFrame[, independent_variable] <- droplevels(tempDataFrame[, independent_variable])
    }
    if(is.factor(tempDataFrame[, independent_variable]))
    {
      independentVariableIsFactor <- TRUE
      independentVariableData <- tempDataFrame[, independent_variable]
      independent_variableLevels <- levels(tempDataFrame[, independent_variable])
      # independent_variable1stLevel <- levels(tempDataFrame[, independent_variable])[1]
      # independent_variable2ndLevel <- levels(tempDataFrame[, independent_variable])[2]
    }
  }
  else
  {
    # The regression families coerce the whole frame to numeric, so a labelled
    # ordinal variable became NA here and the caller had to supply a second
    # column holding the same thing as a number. With a declared order the
    # labels have positions, so the position IS the number and the hand-built
    # column is not needed. Without one, nothing changes: a numeric column
    # coerces to itself.
    if (length(util_split_and_clean(independent_variable_order)) > 1L &&
        independent_variable %in% colnames(tempDataFrame) &&
        !is.numeric(tempDataFrame[, independent_variable]))
    {
      level_order <- util_level_order(tempDataFrame[, independent_variable],
                                      independent_variable_order)
      tempDataFrame[, independent_variable] <- match(
        as.character(tempDataFrame[, independent_variable]), level_order)
    }
    tempDataFrame <- as.data.frame(vapply(tempDataFrame, as.numeric, numeric(nrow(tempDataFrame))))
  }

  originalDataFrame <- tempDataFrame
  if (independentVariableIsFactor)
    tempDataFrame[, independent_variable] <- independentVariableData

  # The columns are split by NAME: the independent variable and the covariates
  # head the table, every other column is a burden. g_start and g_end are the
  # callers' positional description of the same split; they are checked against
  # it, so a layout that drifted stops here instead of testing a covariate as an
  # area or an area as a covariate.
  lead_cols <- c(independent_variable, covariates)
  absent <- setdiff(lead_cols, colnames(tempDataFrame))
  if (length(absent) > 0)
    stop("ERROR: I'm stopping here, the data to associate lack the column(s) ",
         paste(absent, collapse = ", "), ".", call. = FALSE)
  # Service columns (core_service_columns()) identify the row and are neither a
  # predictor nor a burden; they travel in the head, untouched.
  service_cols <- setdiff(intersect(core_service_columns(), colnames(tempDataFrame)),
                          lead_cols)
  burden_cols <- core_data_columns(tempDataFrame, also = lead_cols)
  by_position <- if (g_end >= g_start) colnames(tempDataFrame)[g_start:g_end] else character(0)
  if (!identical(by_position, burden_cols))
    stop("ERROR: I'm stopping here, g_start:g_end points at ",
         length(by_position), " column(s) but ", length(burden_cols),
         " burden column(s) are found by name; the table is not laid out as ",
         "independent variable, covariates, burdens.", call. = FALSE)

  df_head <- tempDataFrame[, c(service_cols, lead_cols), drop = FALSE]

  # drop = FALSE and lapply: with a single burden column (always the case at
  # scope SAMPLE) sapply on a 1-column slice returned a bare vector and the
  # area's name became "burden_values".
  burden_values <- as.data.frame(lapply(tempDataFrame[, burden_cols, drop = FALSE], as.numeric),
                                 check.names = FALSE)

  df_colnames <- colnames(tempDataFrame)

  # the TOTAL column is gone, and with it the last use of depth.
  #
  # It was `apply(burden_values, 1, sum)` - the sum of the per-area aggregates -
  # and `depth_analysis == 2` meant "test only that". But the partition into
  # areas is NOT disjoint: the annotation maps one probe onto several genes and
  # the explode replicates it, so a probe on three genes entered that total three
  # times. It was the composition of aggregates the taxonomy forbids, hiding in
  # the one place where nobody looked for it - and the old depth scale gave it a
  # rung of its own.
  #
  # The quantity it meant to compute - "the whole region class, one number per
  # sample" - is an artefact in its own right now: SCOPE = SAMPLE on that
  # (AREA, SUBAREA), derived from the masked position pivot where every position
  # counts once. It is produced by the same consumer as everything else, so
  # there is nothing left to synthesise here.

  if(grepl("log",transformation_y))
  {
    burden_values <- burden_values + min(burden_values[burden_values>0])
  }
  transformation_y <- as.character(transformation_y)
  if(is.null(transformation_y) | length(transformation_y)==0 | is.na(transformation_y))
    transformation_y <- "none"


  burden_values <- as.data.frame(burden_values)
  df_values_orig <- burden_values
  try(
    {
      # The vocabulary lives in io_transform_apply(). It was written out here and
      # again for the independent variable below, with quantile_<n> handled after
      # this block rather than inside it, so the two copies had already drifted.
      burden_values <- io_transform_apply(burden_values, transformation_y)
    }
  )
  burden_values <- as.data.frame(burden_values)

  if(setequal(burden_values,df_values_orig) & transformation_y !="none")
    transformation_y <- paste0("NA_", transformation_y, sep="")

  # 2026-06-08: universal degenerate-burden filter.
  # Burden columns where the response Y is constant (variance == 0) across
  # samples carry no signal - they produce NaN/garbage stats in every model:
  #   - binomial GLM: MLE diverges (intercept-only fit, NaN coeffs/p-values)
  #   - limma/voom lmFit: 0-effect, NaN t-statistic (misleading "no signal")
  #   - polynomial / gaussian glm: rank-deficient, NaN coefficients
  # Critical for LESIONS @ PROBE where ~92% of probes are all-zero across
  # samples (manifest-aligned pivot, retained for positional join with
  # annotations). Filtering here covers ALL downstream callers in one place:
  # apply_stat_model.R (per-probe foreach) and apply_stat_model_batch.R
  # (limma/voom batch) - their existing variance checks become safety nets.
  if (ncol(burden_values) > 0L) {
    is_degenerate <- vapply(burden_values, function(x) {
      u <- unique(stats::na.omit(as.numeric(x)))
      length(u) < 2L
    }, logical(1))
    if (any(is_degenerate)) {
      core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " io_data_preparation [", family_test, " / ", key$MARKER, " ",
                key$FIGURE, " ", key$AREA, " ", key$SUBAREA,
                "]: dropping ", sum(is_degenerate), "/", length(is_degenerate),
                " degenerate burden columns (var==0).")
      burden_values <- burden_values[, !is_degenerate, drop = FALSE]
    }
  }


  # The independent variable is not transformed here. It used to be, and the
  # result was written into tempDataFrame and then lost when the table was
  # rebuilt from df_head below, so transformation_x never reached a model.
  # assoc_covariates_model() now transforms it upstream into a column of its
  # own, the way it does for each covariate.

  # 2026-06-08: rebuild df_colnames after the degenerate-burden
  # filter above. df_head columns are unchanged (sample-level: IV +
  # covariates); burden_values may now have fewer columns. This replaces
  # the prior strict length check, which fired any time we dropped probes.
  df_colnames <- c(colnames(df_head), colnames(burden_values))
  tempDataFrame <- data.frame(df_head, burden_values)
  if(ncol(tempDataFrame)!=length(df_colnames))
    stop("ERROR: I'm stopping here data are not the same size, file a bug!")

  colnames(tempDataFrame) <- df_colnames
  burden_cols <- colnames(burden_values)

  # After transformation_y some burden values can be NaN (log of 0, ...). They
  # are looked for in the burden columns only and column by column: apply() over
  # the whole table turned it into a character matrix as soon as the independent
  # variable was a factor, and is.nan() then found nothing. The independent
  # variable and the covariates keep their NA, for the models to drop.
  lost_cols <- burden_cols[vapply(tempDataFrame[burden_cols],
                                  function(x) any(is.nan(x)), logical(1))]
  if (length(lost_cols) != 0)
    utils::write.csv2(lost_cols, file.path(ssEnv$session_folder,paste("lost_data_",transformation_y,"_",stringi::stri_rand_strings(1, 12, pattern = "[A-Za-z0-9]"),".log", sep="")))
  for (cl in lost_cols)
    tempDataFrame[is.nan(tempDataFrame[[cl]]), cl] <- 0
  if(family_test=="binomial" | family_test=="binomial_bulk")
    tempDataFrame[, independent_variable] <- as.factor(tempDataFrame[, independent_variable])

  # # remove rows with all NA
  # tempDataFrame <- tempDataFrame[,colSums(is.na(tempDataFrame)) != nrow(tempDataFrame)]

  # 2026-06-09: no more colname sanitisation here. Names stay
  # pass-through from the upstream annotation (HLA-A, chr10:...-..., etc).
  # The per-gene foreach in assoc_apply_stat_model() applies its own LOCAL
  # safe<->real memoised mapping ONLY for the duration of the formula
  # machinery (R formula identifiers cannot contain '-' or ':'), then
  # reverses the mapping before assigning AREA_OF_TEST in the result.
  # CSV ends up with raw names → enrichment downstream resolves HGNC
  # correctly, resume match is exact.

  result <- list(tempDataFrame, independent_variableLevels, burden_cols)
  names(result) <- c("tempDataFrame", "independent_variableLevels", "burden_columns")

  return (result)
}
