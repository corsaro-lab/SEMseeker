#' The pathway keys of a run, built from the hypothesis it was given (internal)
#'
#' AI-261. `util_keys_create()` used to build `keys_for_pathway` at session init
#' by rewriting the registry: every `HYPER` or `HYPO` figure became `HYPER_HYPO`
#' and every window other than `WHOLE` became `ALL_SUBAREAS`. Two consequences,
#' neither of them written anywhere. The enrichment restricted to
#' hypermethylation, and the enrichment restricted to promoters, were not
#' reachable from any argument of any exported function - and the line that
#' looked like it dropped single-figure keys on purpose could never fire, because
#' the rewrite above it had already removed every one of them.
#'
#' The keys are built here instead, from the two answers, at the point where they
#' are asked. One key per marker: the hypothesis is a property of the run, not of
#' each key, so it is stored once and stamped once.
#'
#' `AREA = "GENE"` stays an invariant rather than an argument. A pathway is a set
#' of genes; a cytoband or an island is not a member of one.
#'
#' @param epimutation_direction `"HYPER"`, `"HYPO"` or `"ANY"`.
#' @param gene_region an alias or an explicit vector of windows.
#' @param ssEnv the session environment.
#' @return `ssEnv`, with `keys_for_pathway` and `enrichment_hypothesis` set.
#' @keywords internal
#' @noRd
enrich_pathway_keys_build <- function(epimutation_direction, gene_region, ssEnv) {

  figures <- util_epimutation_direction_expand(epimutation_direction)
  windows <- util_gene_region_expand(gene_region)

  registry <- ssEnv$keys_areas_subareas_markers_figures
  if (is.null(registry) || nrow(registry) == 0)
    stop("enrich_pathway_keys_build(): the run has no annotation registry yet.",
         call. = FALSE)

  available <- registry[core_name_cleaning(as.character(registry$AREA)) == "GENE", ,
                        drop = FALSE]
  if (nrow(available) == 0)
    stop("this run has no GENE areas, so there is nothing a pathway analysis ",
         "could be about. A pathway is a set of genes.", call. = FALSE)

  # A window the run never computed is a request that can never be answered, and
  # it stops here rather than producing an empty gene list that reads like a
  # negative result.
  present <- unique(core_name_cleaning(as.character(available$SUBAREA)))
  missing_windows <- setdiff(as.character(windows), present)
  if (length(missing_windows) > 0)
    stop("gene_region asks for ", paste(missing_windows, collapse = ", "),
         ", which this run did not compute. Available for GENE: ",
         paste(sort(present), collapse = ", "),
         ". An enrichment over a window that was never tested would return an ",
         "empty gene list, which reads like a negative result.", call. = FALSE)

  present_figures <- unique(core_name_cleaning(as.character(available$FIGURE)))
  missing_figures <- setdiff(as.character(figures), present_figures)
  if (length(missing_figures) > 0)
    stop("epimutation_direction asks for ",
         paste(missing_figures, collapse = ", "),
         ", which this run did not compute. Available: ",
         paste(sort(present_figures), collapse = ", "), ".", call. = FALSE)

  markers <- unique(as.character(available$MARKER))

  keys <- data.frame(AREA    = "GENE",
                     SUBAREA = attr(windows, "alias"),
                     MARKER  = markers,
                     FIGURE  = attr(figures, "alias"),
                     stringsAsFactors = FALSE)
  keys$COMBINED <- paste(keys$AREA, keys$SUBAREA, keys$MARKER, keys$FIGURE,
                         sep = "_")

  # The expansions are what the gene-set reader filters on, and what the result
  # is stamped with. They belong to the run because the run answers the two
  # questions once.
  ssEnv$enrichment_hypothesis <- list(
    direction         = attr(figures, "alias"),
    figures           = as.character(figures),
    direction_is_union = isTRUE(attr(figures, "union")),
    region            = attr(windows, "alias"),
    windows           = as.character(windows),
    region_is_union   = length(windows) > 1L)

  ssEnv$keys_for_pathway <- keys
  ssEnv
}

#' The hypothesis this run was given, or a refusal that says how to give one
#'
#' @keywords internal
#' @noRd
enrich_hypothesis_get <- function(ssEnv = NULL) {

  if (is.null(ssEnv))
    ssEnv <- core_get_session_info()

  hypothesis <- ssEnv$enrichment_hypothesis
  if (is.null(hypothesis))
    stop("this session carries no enrichment hypothesis. An enrichment is not ",
         "one analysis but a family of them, and which one you ran is decided ",
         "by `epimutation_direction` and `gene_region`. Enter through a ",
         "function that asks for both.", call. = FALSE)

  hypothesis
}
