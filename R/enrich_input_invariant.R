#' What an enrichment can be about (internal)
#'
#' AI-311. An enrichment over pathways takes a **set of genes**. That is not a
#' filter this package applies for convenience, it is what the analysis is: a
#' pathway is a set of genes, so a CpG island is not a member of one and a
#' collapsed per-sample burden is not a gene. The two coordinates are therefore
#' invariants of the enrichment layer: `SCOPE = INSTANCE`, because a gene list
#' needs one row per gene, and `AREA = GENE`.
#'
#' Regions that are not genes are not excluded from enrichment in general: they
#' are mapped onto genes **first**, by a rule that is itself a biological claim
#' ("nearest TSS", "overlaps the promoter", "within N kb" are three different
#' statements). Once mapped, what reaches this layer is `AREA = GENE` again.
#'
#' It exists as a function, and not as two literals, because it used to be two
#' literals, written out by hand in all five backends. Nothing declared the
#' invariant and nothing enforced it: a sixth backend asking for `ISLAND` would
#' have been served, because [assoc_results_get()] accepts any region class (the
#' cross-study overlaps iterate over all of them legitimately).
#'
#' @return list with `scope` and `area`.
#' @keywords internal
#' @noRd
enrich_input_invariant <- function() list(scope = "INSTANCE", area = "GENE")

#' Refuse an enrichment whose input was never computed
#'
#' AI-311. The invariant above used to be applied as a `subset()` on results
#' already read. A result folder holding no `SCOPE = INSTANCE`, `AREA = GENE`
#' row therefore produced *nothing*, quietly: every backend received zero rows
#' and wrote no result. A reader of that folder sees "no significant
#' enrichment", which is a statement about the data, when the truth is "the
#' input was never computed", which is a statement about the pipeline.
#'
#' That gap widened when the two aggregation branches stopped being additive
#' (AI-308): `association_analysis(scope = "SAMPLE")` is now a legitimate and
#' complete request, and it produces a folder on which every enrichment is
#' silently empty.
#'
#' So the check moves to the door and names the remedy. Two failure modes:
#' \itemize{
#'   \item the run has no `GENE` region class at all: nothing a pathway
#'     analysis could be about;
#'   \item it has one, but no inference result carries a per-instance gene row:
#'     the association was run, on something else.
#' }
#' A result file written before `SCOPE` existed as a column cannot be judged and
#' is not blocked: absence of evidence is not the same as evidence of absence,
#' and blocking on it would refuse folders that are merely old.
#'
#' @param inference_details the requests the enrichment is about to run.
#' @param ssEnv session environment; read from the session when omitted.
#' @return invisibly `TRUE`; called for the refusal.
#' @keywords internal
#' @noRd
enrich_input_assert <- function(inference_details, ssEnv = NULL) {

  if (is.null(ssEnv)) ssEnv <- core_get_session_info()
  want <- enrich_input_invariant()

  registry <- ssEnv$keys_areas_subareas
  areas <- if (is.null(registry)) character(0) else
    core_name_cleaning(as.character(registry$AREA))
  if (!(want$area %in% areas))
    stop("this run declares no ", want$area, " region class, so there is ",
         "nothing a pathway analysis could be about: a pathway is a set of ",
         "genes. Declare it on the analysis: ",
         "association_analysis(areas = \"GENE\", scope = \"INSTANCE\", ...), ",
         "and run the enrichment on the result.", call. = FALSE)

  markers <- unique(as.character(ssEnv$keys_markers_figures$MARKER))
  judged <- FALSE

  for (z in seq_len(nrow(inference_details))) {
    for (marker in markers) {
      f <- io_inference_file_name(inference_details[z, ], marker,
                                  ssEnv$result_folderInference)
      if (!file.exists(f) || file.info(f)$size <= 10) next

      header <- tryCatch(utils::read.csv2(f, nrows = 1, stringsAsFactors = FALSE),
                         error = function(e) NULL)
      if (is.null(header) || !all(c("SCOPE", "AREA") %in% colnames(header))) next

      # A file written before the coordinates were columns: it cannot answer,
      # so it does not get to condemn either.
      judged <- TRUE
      hit <- tryCatch({
        lazy <- polars::pl$scan_csv(f, separator = ";", decimal_comma = TRUE,
                                    null_values = "NA",
                                    infer_schema_length = 100000L)
        nrow(as.data.frame(
          lazy$select(c("SCOPE", "AREA"))$filter(
            polars::pl$col("SCOPE")$eq(want$scope)$and(
              polars::pl$col("AREA")$eq(want$area)))$head(1L)$collect()))
      }, error = function(e) NA_integer_)

      if (is.na(hit)) { judged <- FALSE; next }
      if (hit > 0) return(invisible(TRUE))
    }
  }

  if (!judged) return(invisible(TRUE))

  stop("no inference result of this folder carries a gene row: the enrichment ",
       "reads SCOPE = ", want$scope, " and AREA = ", want$area,
       ", and every result here was produced on something else. This is not a ",
       "negative result: the input was never computed. Produce it with ",
       "association_analysis(scope = \"INSTANCE\", areas = \"GENE\", ...); a ",
       "run at scope = \"SAMPLE\" answers a different question and cannot feed ",
       "a pathway analysis.", call. = FALSE)
}
