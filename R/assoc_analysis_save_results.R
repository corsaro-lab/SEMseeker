assoc_analysis_save_results <- function(results=NULL,fileNameResults, family_test, filter_p_value, append=FALSE ){

  if(nrow(results)==0)
    return()

  ssEnv <- core_get_session_info()

  # C-06: stamp provenance columns so any downstream CSV stack retains build+tech
  genome_build_val <- if (!is.null(ssEnv$genome_build) && nzchar(ssEnv$genome_build))
    as.character(ssEnv$genome_build) else "hg19"
  tech_val <- if (!is.null(ssEnv$tech) && nzchar(ssEnv$tech))
    as.character(ssEnv$tech) else ""
  results$GENOME_BUILD <- genome_build_val
  results$TECH         <- tech_val

  utils::write.csv2(.assoc_key_first(results), fileNameResults, row.names = FALSE)
  multiple_test_adj <- core_name_cleaning(ssEnv$multiple_test_adj)
  # there is a bug which mantain more family test in the same results file
  # so we need to filter the results
  #
  colnames(results) <- core_name_cleaning(colnames(results))
  results <- subset(results, FAMILY_TEST==as.character(family_test))

  # check if results is empty
  if(is.null(results))
    return()

  if(nrow(results)==0)
    return()

  if (!append)
  {
    results <- results[,!grepl("SAMPLES_SQL_CONDITION", colnames(results))]
    results <- unique(results)

    pvalue_columns <- colnames(results)[grepl("PVALUE", colnames(results)) & !grepl("_ADJ", colnames(results))]

    # AI-257: drop every adjusted column of a previous pass before recomputing.
    # `PVALUE_ADJ` is in the list because it is the name the three levels
    # replaced: a file written before this release carries it, and leaving it
    # there would put a column whose family nobody can name next to three that
    # say theirs.
    stale_adj <- grepl("_ADJ_KEY_|_ADJ_SCOPE_|_ADJ_ALL_", colnames(results)) |
      colnames(results) == "PVALUE_ADJ"
    results <- unique(results[, !stale_adj, drop = FALSE])

    # AI-257: the FDR is controlled over three nested families, and each column
    # says which one it belongs to. See .assoc_adjust_levels() for why there are
    # three and not one.
    results <- .assoc_adjust_levels(results, method = ssEnv$multiple_test_adj,
                                    method_label = multiple_test_adj,
                                    alpha = ssEnv$alpha)

    if (exists("results") & length(pvalue_columns)>0)
    {
      for (p in seq_along(pvalue_columns))
      {
        col_p <- core_name_cleaning(paste0(pvalue_columns[p], "_ADJ_ALL_", multiple_test_adj))
        results[, col_p] <- .assoc_adjust_vector(results[, pvalue_columns[p]],
                                                 method = ssEnv$multiple_test_adj,
                                                 alpha = ssEnv$alpha,
                                                 what = col_p)
        colnames(results) <- core_name_cleaning(colnames(results))
      }

      # AI-257: order by the widest level, named outright. It used to be
      # `grepl(multiple_test_adj, colnames)[1]` - the first column whose name
      # merely CONTAINS the method - which with three adjusted levels picks
      # whichever one happens to come first in the frame.
      pvalue_adj_colname <- core_name_cleaning(paste0("PVALUE_ADJ_ALL_", multiple_test_adj))

      if (nrow(results)>0 && pvalue_adj_colname %in% colnames(results))
        results <- results[order(results[,pvalue_adj_colname]),]

    }

    if(nrow(results)==0)
      return()

    # AI-255: DEPTH is not stamped any more. It was a number standing in for
    # what SCOPE, AREA and SUBAREA now say outright, and it stood in badly: the
    # 1/2/3 ladder projected a partial order onto a line, and its rung 2 marked
    # rows produced by composing aggregates - a quantity that no longer exists.
    # Nothing reads it to decide anything, so writing it would only invite
    # someone to start.
    #
    # Gone with it: `results[is.na(results$SUBAREA), "SUBAREA"] <- "TOTAL"`.
    # TOTAL was the label of that synthesis, and it is not a value of the
    # SUBAREA vocabulary - filling a missing coordinate with an invented one
    # hides the defect inside the key instead of showing it.
    # replace empty with NA
    results[results == ""] <- NA
    results[results == " "] <- NA
    # remove columns where all rows are NA
    results <- results[, colSums(is.na(results)) < nrow(results)]

    # check if exists at least a column with PVALUE
    if(!any(grepl("PVALUE", colnames(results))))
      return()
  }

  # AI-257: select the columns by name, not by `grepl(method, colnames)`.
  # That predicate matched every column whose name merely CONTAINS the method
  # string, so adding PVALUE_ADJ_KEY_BH and PVALUE_ADJ_SCOPE_BH next to
  # PVALUE_ADJ_ALL_BH would silently turn this flag into "significant at all
  # three levels at once" - a change of meaning invisible in the diff, and one
  # nobody asked for. The flag answers for the widest family, and now says so.
  results$SIGNIFICATIVE_ADJ_ALL <- .assoc_all_below(
    results, .assoc_level_columns(results, "ALL", multiple_test_adj), ssEnv$alpha)
  results$SIGNIFICATIVE <- .assoc_all_below(
    results, colnames(results)[grepl("PVALUE", colnames(results)) &
                                 !grepl("_ADJ", colnames(results))], ssEnv$alpha)

  # AI-257: this used to read SIGNIFICATIVE_ADJ, a column this function never
  # creates - it is built by assoc_results_get() on the way out, not on the way
  # in. `filter_p_value` defaults to TRUE (assoc_analysis.R:124-125), so the
  # filter either stopped the run or was never reached. It filters on the flag
  # that exists here.
  if(filter_p_value && "SIGNIFICATIVE_ADJ_ALL" %in% colnames(results))
    results <- subset(results, SIGNIFICATIVE_ADJ_ALL)

  # remove duplicates based on MARKER   FIGURE  AREA    SUBAREA AREA_OF_TEST    FAMILY_TEST TRANSFORMATION_Y    PVALUE  R_MODEL
  # calculating the max of all others columns
  # C-06: include provenance columns in the grouping key so summarise() preserves them
  # AI-248: AGGREGATION is part of the identity. Without it the summarise(max)
  # below would fuse the median and the mean of the same scope into one row.
  # AI-255: SCOPE belongs here too. Without it the collapsed row and the
  # per-instance row of the same region class would be fused by the summarise
  # below - the very thing the aggregation axis was added to prevent, one
  # coordinate further along.
  group_column <- c("MARKER", "FIGURE", "SCOPE", "AGGREGATION", "AREA", "SUBAREA", "AREA_OF_TEST", "FAMILY_TEST",
                    "TRANSFORMATION_Y", "R_MODEL", "TRANSFORMATION_X",
                    "INDEPENDENT_VARIABLE", "COVARIATES",
                    "GENOME_BUILD", "TECH")
  group_column <- group_column[group_column %in% colnames(results)]

  if(ncol(results[,!colnames(results) %in% group_column])>2)
  {
    # AI-061+ (2026-06-09): use base |> pipe (R 4.1+) instead of %>%.
    # The %>% reference was unresolved at runtime (no @importFrom magrittr)
    # and caused "could not find function %>%" mid-association on ewas v31.
    results <- results |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_column))) |>
      dplyr::summarise(dplyr::across(dplyr::everything(),
                                     ~ max(.x, na.rm = TRUE)),
                       .groups = 'drop')

    utils::write.csv2(.assoc_key_first(results), fileNameResults, row.names = FALSE)
  }
}

#' The three nested families a p-value is adjusted over (internal)
#'
#' AI-257. The family over which an FDR is controlled is a statistical choice,
#' and it has to be readable from the column that carries the result. Until this
#' release two things made it something else:
#'
#' \itemize{
#'   \item `assoc_apply_stat_model()` adjusted whatever `result_temp` it was
#'     handed, and at `SCOPE = INSTANCE` that is one chunk of the pivot - the
#'     family was the memory split;
#'   \item the single `PVALUE_ADJ` said nothing about its family, and at
#'     `SCOPE = SAMPLE` that family holds **one row**, so the column was equal to
#'     `PVALUE` by arithmetic while its name promised an adjustment.
#' }
#'
#' Three columns now, three families, nested, each named after its own:
#'
#' \preformatted{
#' PVALUE_ADJ_KEY_<m>    the identity key - MARKER, FIGURE, SCOPE, AREA,
#'                       SUBAREA, AGGREGATION. Members: the instances.
#' PVALUE_ADJ_SCOPE_<m>  every row of the same SCOPE. Members: the region
#'                       classes, figures and aggregations of that scope.
#' PVALUE_ADJ_ALL_<m>    the whole file. Members: everything tested on this
#'                       marker under this model.
#' }
#'
#' **Why the middle one exists.** At `SCOPE = SAMPLE` the key holds a single row,
#' so `KEY` is the raw p-value, and `ALL` puts that row in a family dominated by
#' the tens of thousands of per-instance rows it shares the file with. Neither
#' answers "how many things did I test on this sample"; `SCOPE` does.
#'
#' **Why each family is defined the same way on every row.** The tempting design
#' is a single column whose family adapts - the instances where there are
#' instances, something else where there are none. That is one name meaning two
#' things depending on the row, which is what `depth` did and what AI-255 removed
#' it for. Each of the three is one rule, applied uniformly; at
#' `SCOPE = INSTANCE` the `SCOPE` level pools genes with cytobands and islands,
#' which is a number to read with that in mind, not a number to hide.
#'
#' **They are not a severity ladder.** Benjamini-Hochberg is adaptive: enlarging
#' a family with strongly significant tests lowers the ratio `n/j` for the
#' others, since `(n+k)/(i+k) <= n/i` for every `i <= n`. A wider family is
#' usually but not necessarily more conservative. The three answer three
#' questions, and the analysis has to say which one it used.
#'
#' @param results the accumulated results of the run.
#' @param method the adjustment method as the run declared it (`ssEnv$multiple_test_adj`).
#' @param method_label the same, cleaned, as it appears in column names.
#' @param alpha the significance level, for the `q` estimator.
#' @return `results` with the `KEY` and `SCOPE` levels recomputed.
#' @keywords internal
#' @noRd
.assoc_adjust_levels <- function(results, method, method_label, alpha) {

  if (is.null(results) || nrow(results) == 0 || !("PVALUE" %in% colnames(results)))
    return(results)

  # AI-257: the two narrow levels are computed on PVALUE only, not on every
  # column matching "PVALUE". A file already carries one adjusted column per
  # p-value column at the ALL level, and these files reach 600-880 MB on a
  # 485k-probe run (AI-078); tripling that width to adjust the intercept's
  # p-value buys nothing. PVALUE is the model's own p-value - the one the
  # taxonomy is about and the one the enrichment reads.
  levels <- list(
    list(name = "KEY",
         cols = c("MARKER", "FIGURE", "SCOPE", "AREA", "SUBAREA", "AGGREGATION")),
    list(name = "SCOPE", cols = "SCOPE"))

  pvalues <- suppressWarnings(as.numeric(results$PVALUE))

  for (level in levels) {
    family_cols <- intersect(level$cols, colnames(results))
    if (length(family_cols) == 0)
      next

    col_name <- core_name_cleaning(paste0("PVALUE_ADJ_", level$name, "_",
                                          method_label))
    families <- do.call(paste, c(lapply(family_cols, function(cl)
      as.character(results[[cl]])), list(sep = "\r")))

    adjusted <- rep(NA_real_, length(pvalues))
    for (fam in unique(families)) {
      idx <- which(families == fam)
      adjusted[idx] <- .assoc_adjust_vector(pvalues[idx], method = method,
                                            alpha = alpha, what = col_name)
    }
    results[[col_name]] <- adjusted
  }

  results
}

#' Adjust a vector of p-values by the estimator the run declared (internal)
#'
#' AI-257. One estimator for all three levels, so that the numbers in the three
#' columns differ by their family and by nothing else. Before this, the narrow
#' level had `method = "BH"` written into it while the global level honoured
#' `ssEnv$multiple_test_adj`: a run asking for `q` got a column named for BH's
#' sibling and computed as BH, with nothing to show for it.
#'
#' `qvalue()` needs enough p-values to estimate pi0, and the narrow families are
#' small by construction - a `SCOPE = SAMPLE` key holds one row. When it cannot
#' estimate, the answer is `NA` and a line in the log: a family too small for the
#' estimator is a fact about the request, and substituting a different estimator
#' would make the column name a lie.
#'
#' @param pvalues numeric vector.
#' @param method the method as declared (`BH`, `fdr`, `q`, ...).
#' @param alpha significance level, used by the `q` estimator.
#' @param what the column being computed, for the log line.
#' @return numeric vector of the same length, `NA` where it cannot be estimated.
#' @keywords internal
#' @noRd
.assoc_adjust_vector <- function(pvalues, method, alpha, what = "") {

  pvalues <- suppressWarnings(as.numeric(pvalues))
  if (length(pvalues) == 0)
    return(numeric(0))

  if (identical(core_name_cleaning(method), "Q")) {
    estimated <- try(qvalue::qvalue(pvalues, fdr.level = as.numeric(alpha),
                                    pi0.method = "bootstrap",
                                    na.rm = TRUE)$qvalues, silent = TRUE)
    if (inherits(estimated, "try-error")) {
      core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " ", what, ": qvalue could not estimate pi0 on a family of ",
                length(pvalues), " - left NA rather than substituting another ",
                "estimator, which would make the column name wrong.")
      return(rep(NA_real_, length(pvalues)))
    }
    return(as.numeric(estimated))
  }

  stats::p.adjust(pvalues, method = method)
}

#' The adjusted columns of one level (internal)
#'
#' AI-257. `PVALUE_ADJ_ALL_BH` and `I_AGE_PVALUE_ADJ_ALL_BH` both belong to the
#' ALL level; `PVALUE_ADJ_KEY_BH` does not. Selecting them by the level's own
#' infix is what keeps a new level from being swept into an existing flag.
#'
#' @param results the results data.frame.
#' @param level `"KEY"`, `"SCOPE"` or `"ALL"`.
#' @param method_label the method as it appears in column names.
#' @return the names of the matching columns, possibly none.
#' @keywords internal
#' @noRd
.assoc_level_columns <- function(results, level, method_label) {
  infix <- core_name_cleaning(paste0("_ADJ_", level, "_", method_label))
  colnames(results)[endsWith(colnames(results), infix)]
}

#' Is every one of these columns below alpha, row by row (internal)
#'
#' AI-257. Returns `NA` when there is no column to answer with, instead of the
#' `TRUE` that `all(logical(0))` yields - "every one of no columns is
#' significant" is how an empty selection used to pass for a positive result.
#'
#' @param results the results data.frame.
#' @param columns the columns to test.
#' @param alpha significance level.
#' @return logical vector, one entry per row.
#' @keywords internal
#' @noRd
.assoc_all_below <- function(results, columns, alpha) {

  if (is.null(results) || nrow(results) == 0)
    return(logical(0))

  columns <- intersect(columns, colnames(results))
  if (length(columns) == 0)
    return(rep(NA, nrow(results)))

  values <- as.data.frame(results[, columns, drop = FALSE])
  apply(values, 1, function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (all(is.na(x))) NA else all(x < as.numeric(alpha), na.rm = TRUE)
  })
}

#' Put the taxonomy key first, and refuse a key with a hole in it (internal)
#'
#' AI-255. The six coordinates plus `AREA_OF_TEST` - the instance within the
#' artefact - are the identity of a result row, so they lead the file. A reader
#' opening the CSV sees what the row *is* before seeing what was measured on it,
#' and the column order stops depending on the order in which the models happened
#' to add their fields.
#'
#' The NA check is the other half. Filling a missing coordinate used to be normal
#' here - `results[is.na(results$SUBAREA), "SUBAREA"] <- "TOTAL"` invented a
#' value that is not in the SUBAREA vocabulary - and that is exactly how a defect
#' hides inside a key: two rows that cannot be told apart, and nothing to show
#' for it. With every coordinate composed by one function there is no legitimate
#' NA left, so an NA means something upstream did not set what it was supposed
#' to, and it stops here rather than travelling into a dedup or a join.
#'
#' @param results the results data.frame.
#' @return `results` with the key columns first.
#' @keywords internal
#' @noRd
.assoc_key_first <- function(results) {

  if (is.null(results) || nrow(results) == 0)
    return(results)

  key_cols <- c("MARKER", "FIGURE", "SCOPE", "AREA", "SUBAREA", "AGGREGATION",
                "AREA_OF_TEST")
  present <- intersect(key_cols, colnames(results))
  if (length(present) == 0)
    return(results)

  # An empty SUBAREA has always meant "the whole area" - the same normalisation
  # io_pivot_file_name() applies. Settle it before the check, so the convention
  # is honoured and what remains missing is genuinely missing.
  if ("SUBAREA" %in% present) {
    blank <- is.na(results$SUBAREA) | !nzchar(trimws(as.character(results$SUBAREA)))
    if (any(blank)) results$SUBAREA[blank] <- "WHOLE"
  }

  holed <- present[vapply(present, function(cl)
    any(is.na(results[[cl]]) | !nzchar(trimws(as.character(results[[cl]])))),
    logical(1))]
  if (length(holed) > 0)
    stop("the taxonomy key of a result row is incomplete: ",
         paste(holed, collapse = ", "),
         " carries missing values. The key is the identity of the row - two ",
         "rows with a hole in the same place cannot be told apart - so this is ",
         "a coordinate that was never set upstream, not a value to fill in.",
         call. = FALSE)

  results[, c(present, setdiff(colnames(results), present)), drop = FALSE]
}
