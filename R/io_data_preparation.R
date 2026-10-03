#' Title
#'
#' @param family_test test or regression to apply
#' @param transformation_y transformation_y to apply to data
#' @param tempDataFrame data frame to use for test/regression
#' @param independent_variable regressor
#' @param g_start index of the first burden column in tempDataFrame
#' @param g_end index of the last burden column in tempDataFrame
#' @param covariates vector of covariates to be found in the sample sheet
#' @param key named list with AREA, SUBAREA, MARKER and FIGURE identifiers, used
#'   to name the artefact in the log when degenerate columns are dropped
#' @param transformation_x transformation to apply to the independent variable
#' @param independent_variable_order optional order of the levels of an ordinal
#'   independent variable, "+"-separated. See util_level_order().
#'   before the fit; "none" leaves it untouched
#'
#' @return A named list with two elements: \code{tempDataFrame} (the prepared
#'   and optionally transformed data.frame) and \code{independent_variableLevels}
#'   (the factor levels of the independent variable, or \code{NULL} for continuous
#'   outcomes).
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

  # drop = FALSE: when g_start == 2 (only one head column = IV) the default
  # 1D slice returns a vector, colnames(vec) = NULL, and the length check at
  # the bottom (`ncol(tempDataFrame) != length(df_colnames)`) fires the
  # "data are not the same size" stop. Forcing data.frame keeps the rebuild
  # symmetric regardless of how many sample-level columns there are.
  df_head <- tempDataFrame[, seq_len(g_start - 1), drop = FALSE]

  burden_values <- sapply(tempDataFrame[,g_start:g_end], as.numeric)
  burden_values <- as.data.frame(burden_values)

  df_colnames <- colnames(tempDataFrame)

  # AI-255: the TOTAL column is gone, and with it the last use of depth.
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
      burden_values <- switch(
        as.character(transformation_y),
        "scale"  = scale(burden_values),
        "log"    = log(burden_values),
        "log2"   = log2(burden_values),
        "log10"  = log10(burden_values),
        "exp"    = exp(burden_values),
        "factor" = as.data.frame(lapply(burden_values, as.factor)),   # AI-044 binomial Y as factor (0/1)
        "none"   = burden_values,
        burden_values
      )
    }
  )
  if(grepl("quantile", transformation_y))
  {
    qq <- as.numeric(unlist(strsplit(transformation_y,"\\_"))[2])
    burden_values <- as.data.frame(apply(burden_values,2,function(x){
      if(length(unique(x))>=qq)
        as.numeric(dplyr::ntile(x, n=qq))
      else
        rep(0,length(x))
    }))
  }
  burden_values <- as.data.frame(burden_values)

  if(setequal(burden_values,df_values_orig) & transformation_y !="none")
    transformation_y <- paste0("NA_", transformation_y, sep="")

  # AI-044 (2026-06-08): universal degenerate-burden filter.
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


  if(family_test!="binomial" & family_test!="binomial_bulk" & family_test!="wilcoxon" & family_test!="jsd" & family_test!="t.test" & family_test!="poisson" &
      family_test!="chisq.test" & family_test!="fisher.test" & family_test!="kruskal.test")
  {
    variable_to_transform <- independent_variable
    if(length(covariates)>0)
    {
      variable_to_transform <- c(independent_variable,covariates)
      independent_variableValues <-as.data.frame(apply(tempDataFrame[,variable_to_transform] ,2, as.numeric))
    }
    else
    {
      independent_variableValues <-as.data.frame(as.numeric(tempDataFrame[,variable_to_transform]))
    }
    independent_variableValuesOrig <- independent_variableValues
    suppressWarnings(
      try(
        {
          # AI-044 (2026-06-08): questo switch ora usa `transformation_x`
          # (era `transformation_y` per legacy reuse) -> separazione semantica
          # Y vs X. Aggiunto case "factor" per encode l'IV come categorical
          # (utile per glm binomial: OR per livello vs reference).
          independent_variableValues <- switch(
            as.character(transformation_x),
            "scale"  = scale(independent_variableValues),
            "log"    = log(independent_variableValues),
            "log2"   = log2(independent_variableValues),
            "log10"  = log10(independent_variableValues),
            "exp"    = exp(independent_variableValues),
            "factor" = as.data.frame(lapply(independent_variableValues, as.factor)),
            "none"   = independent_variableValues,
            independent_variableValues
          )
        }
      )
    )



    # if(grepl("quantile", transformation_y))
    # {
    #   qq <- unlist(strsplit(transformation_y,"_")[2])
    #   df_values_temp <- as.data.frame(apply( burden_values,2,function(x) dplyr::ntile(x, n=qq)))
    #   colnames(df_values_temp) <- colnames(burden_values)
    # }

    if(setequal(burden_values,df_values_orig) & transformation_y !="none")
      transformation_y <- paste0("NA_", transformation_y, sep="")
    else
      tempDataFrame[, variable_to_transform] <- independent_variableValues
  }


  # AI-044 (2026-06-08): rebuild df_colnames after the degenerate-burden
  # filter above. df_head columns are unchanged (sample-level: IV +
  # covariates); burden_values may now have fewer columns. This replaces
  # the prior strict length check, which fired any time we dropped probes.
  df_colnames <- c(colnames(df_head), colnames(burden_values))
  tempDataFrame <- data.frame(df_head, burden_values)
  if(ncol(tempDataFrame)!=length(df_colnames))
    stop("ERROR: I'm stopping here data are not the same size, file a bug!")

  colnames(tempDataFrame) <- df_colnames
  # after the transformation_y some data could be missed
  lost_cols <- colSums(apply(tempDataFrame,2,is.nan))!=0
  lostDataFrame <-  colnames(tempDataFrame)[lost_cols]
  if(sum(lost_cols)!=0)
    utils::write.csv2(lostDataFrame, file.path(ssEnv$session_folder,paste("lost_data_",transformation_y,"_",stringi::stri_rand_strings(1, 12, pattern = "[A-Za-z0-9]"),".log", sep="")))

  #  we want to preserve the NA in the independent variables to be removed by the models
  tempDataFrame[apply(tempDataFrame,2,is.nan)] <- 0
  if(family_test=="binomial" | family_test=="binomial_bulk")
    tempDataFrame[, independent_variable] <- as.factor(tempDataFrame[, independent_variable])

  # # remove rows with all NA
  # tempDataFrame <- tempDataFrame[,colSums(is.na(tempDataFrame)) != nrow(tempDataFrame)]

  # AI-106 (2026-06-09): no more colname sanitisation here. Names stay
  # pass-through from the upstream annotation (HLA-A, chr10:...-..., etc).
  # The per-gene foreach in assoc_apply_stat_model() applies its own LOCAL
  # safe<->real memoised mapping ONLY for the duration of the formula
  # machinery (R formula identifiers cannot contain '-' or ':'), then
  # reverses the mapping before assigning AREA_OF_TEST in the result.
  # CSV ends up with raw names → enrichment downstream resolves HGNC
  # correctly, resume match is exact.

  result <- list(tempDataFrame, independent_variableLevels)
  names(result) <- c("tempDataFrame", "independent_variableLevels")

  return (result)
}
