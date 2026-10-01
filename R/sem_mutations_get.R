#' sem_mutations_get
#'
#' @param values values of methylation - data.frame with columns CHR, START, END
#'   and a fourth numeric VALUE column.
#' @param figure figure to get Mutations of HYPO or HYPER methylation
#' @param thresholds threshold to use for comparison - data.frame with columns
#'   CHR, START, END, signal_inferior_thresholds, signal_superior_thresholds.
#' @param sampleName name of the sample
#'
#' @return mutations data.frame with columns CHR, START, END, MUTATIONS (0/1).
#'   Only positions present in BOTH values and thresholds are returned
#'   (inner join on CHR, START, END via util_join_values_to_thresholds).
#'
sem_mutations_get <- function(values, figure, thresholds, sampleName) {

  # Polars inner join on (CHR, START, END).
  # Coverage banner is emitted once per batch by sem_analyze_population() before
  # the per-sample loop - not repeated here to avoid log noise.
  joined <- util_join_values_to_thresholds(values, thresholds)

  core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
    " [sem_mutations_get] sample=", sampleName, " figure=", figure,
    " covered=", nrow(joined), "/", nrow(values))

  # ── Empty-result guard ──────────────────────────────────────────────────────
  if (nrow(joined) == 0L) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
      " [sem_mutations_get] No overlapping positions - returning empty result",
      " for sample=", sampleName)
    return(data.frame(CHR = character(), START = integer(),
                      END = integer(),   MUTATIONS = integer(),
                      stringsAsFactors = FALSE))
  }

  # ── Mutation call ────────────────────────────────────────────────────────────
  if (figure == "HYPO") {
    mutation <- as.numeric(joined$VALUE < joined$signal_inferior_thresholds)
  } else {
    mutation <- as.numeric(joined$VALUE > joined$signal_superior_thresholds)
  }

  mutation_annotated_sorted <- anno_sort_by_chr_and_start(data.frame(
    CHR       = joined$CHR,
    START     = joined$START,
    END       = joined$END,
    MUTATIONS = mutation,
    stringsAsFactors = FALSE
  ))

  return(mutation_annotated_sorted)
}
