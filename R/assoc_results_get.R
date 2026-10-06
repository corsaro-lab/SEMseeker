# NOTE: naming convention - `adjustment_method` (singular) is used in internal
# helpers and result-retrieval functions where a single method string is expected.
# High-level public functions (e.g. enrichment_analysis()) use `adjustment_methods`
# (plural) because they accept a vector to iterate over multiple corrections.
#
# `area` and `scope` have no defaults. `area = "GENE"` used to be one,
# and a default is the wrong shape for this: which region class a caller wants is
# never obvious from the outside, and a caller that forgot to say got the genes
# and no indication that it had chosen anything. Every call site now declares
# what it reads. The enrichment does not call this directly at all - it goes
# through enrich_gene_set_get(), where GENE and INSTANCE are invariants rather
# than arguments.
#
# `adjust_per_area` and `adjust_globally` are gone. They re-ran
# p.adjust() at READ time, on top of a column that is already adjusted -
# `pvalue_column` defaults to PVALUE_ADJ_ALL_BH, so switching one on meant BH
# over BH. No call site in the package ever passed TRUE. `adjust_per_area`
# additionally reused the name `area` as its loop counter, shadowing the
# parameter, so the final `subset(AREA == area)` filtered on the last area of
# the loop instead of the requested one - the function would have returned a
# different region class than the one asked for.
assoc_results_get <- function (inference_detail, marker,
  pvalue_column="PVALUE_ADJ_ALL_BH", adjustment_method = "BH", area, scope,
  omit_na = TRUE, significance = NULL)
{

  if (missing(area) || is.null(area) || !nzchar(as.character(area)))
    stop("assoc_results_get(): 'area' is required - name the region class to ",
         "read (e.g. \"GENE\"). It used to default to GENE, which is how a ",
         "caller could read one class while meaning another.", call. = FALSE)
  if (missing(scope) || is.null(scope) || !nzchar(as.character(scope)))
    stop("assoc_results_get(): 'scope' is required - \"INSTANCE\" for one row ",
         "per gene, island or probe, \"SAMPLE\" for the collapsed artefact. ",
         "The two are different quantities and no default can pick between ",
         "them.", call. = FALSE)

  ssEnv <- core_get_session_info()
  resultFolder <- ssEnv$result_folderInference

  inferenceFile <- io_inference_file_name(inference_detail, marker, ssEnv$result_folderInference)

  if(!file.exists(inferenceFile))
  {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),  " Inference file does not exist: ", inferenceFile)
    return(data.frame())
  }


  results_inference <- utils::read.csv2(inferenceFile, row.names = NULL, header = TRUE, stringsAsFactors = FALSE)
  colnames(results_inference) <- core_name_cleaning(colnames(results_inference))
  pvalue_column <- core_name_cleaning(pvalue_column)

  # pvalue_columns in the user setup is a single vector applied to
  # every (inference_detail × marker) combination, but inference CSVs are
  # per-IV: a setup with IVs {TUMOUR_STAGE_N, BIOLOGICAL_RANK} writes
  # 'I_TUMOUR_STAGE_N_..._PVALUE_ADJ_ALL_FDR' to one CSV and
  # 'I_BIOLOGICAL_RANK_..._PVALUE_ADJ_ALL_FDR' to the other. Stopping the
  # whole run when one column happens to belong to the wrong IV kills the
  # entire enrichment pipeline on the first irrelevant lookup. Return an
  # empty result instead so the for-pvalue_columns loop in
  # enrichment_analysis() just skips this combination and keeps going.
  if((!pvalue_column %in% colnames(results_inference)))
  {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " ", pvalue_column,
              " column does not exist in inference file (skipping): ",
              inferenceFile)
    return(data.frame())
  }

  # remove rows ehere pvalue_column is inf or -inf
  results_inference <- results_inference[!is.infinite(results_inference[,pvalue_column]),]

  # sem_metrics_name_collect() removed from here. Its body is commented
  # out in full, so the call did nothing - but what it used to do is the reason
  # it does not belong on a read path: it wrote a metrics registry to disk while
  # a consumer was reading. A reader reads.
  multiple_test_adj <- core_name_cleaning(ssEnv$multiple_test_adj)
  # this reader is free - it returns the artefacts the caller asks for.
  # It replaces `subset(DEPTH == 3)`, which said "per instance" through a number
  # whose meaning had to be remembered, with the coordinates that say it.
  #
  # `area` was declared as a parameter but never filtered a single row: it was
  # used only to rewrite AREA_OF_TEST a few lines below, so every caller got the
  # AREA_OF_TEST of every other class as well. Writing the predicate out is what
  # made that visible.
  #
  # The invariant SCOPE = INSTANCE & AREA = GENE belongs to enrichment, and is
  # enforced where enrichment enters - not here, where the cross-study overlaps
  # legitimately iterate over every area of the registry.
  results_inference <- results_inference[which(results_inference$AREA == area), , drop = FALSE]
  if (!is.null(scope) && "SCOPE" %in% colnames(results_inference))
    results_inference <- results_inference[which(results_inference$SCOPE == scope), , drop = FALSE]
  # name the level instead of matching the method string. `grepl("BH",
  # colnames)` caught every adjusted column at once, so the flag silently became
  # an AND across levels the moment a second one existed. SIGNIFICATIVE_ADJ
  # answers for the widest family, the same one `pvalue_column` defaults to.
  results_inference$SIGNIFICATIVE_ADJ <- .assoc_all_below(
    results_inference,
    .assoc_level_columns(results_inference, "ALL", multiple_test_adj),
    ssEnv$alpha)
  results_inference$SIGNIFICATIVE <- .assoc_all_below(
    results_inference,
    colnames(results_inference)[grepl("PVALUE", colnames(results_inference)) &
                                  !grepl("_ADJ", colnames(results_inference))],
    ssEnv$alpha)
  # results_inference <- results_inference[,c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column,"PVALUE","DEPTH")]

  results_inference[results_inference$AREA==area,"AREA_OF_TEST"] <- gsub(results_inference[results_inference$AREA==area,"AREA_OF_TEST"] , pattern="_", replacement="-")

  # the re-adjustment at read time is gone, and with it the second
  # `subset(AREA == area)` that used to follow it - the first one, above, is the
  # filter, and repeating it here only mattered because the loop in between had
  # overwritten `area`.

  # preserve only subareas selected
  results_inference <- results_inference[results_inference$SUBAREA %in% unique(ssEnv$keys_areas_subareas$SUBAREA),]
  if(!is.null(significance))
  {
    if (significance)
      results_inference <- subset(results_inference, results_inference[,pvalue_column] < as.numeric(ssEnv$alpha))
    else
      results_inference <- subset(results_inference, results_inference[,pvalue_column] >= as.numeric(ssEnv$alpha))
  }

  # remove where pvalue_column is NA
  results_inference <- results_inference[!is.na(results_inference[,pvalue_column]),]
  # remove where pvalue_column is -Inf or +Inf
  results_inference <- results_inference[results_inference[,pvalue_column] != -Inf,]
  results_inference <- results_inference[results_inference[,pvalue_column] != Inf,]

  # if(omit_na)
  #   results_inference <- na.omit(results_inference)

  results_inference <- assoc_filter_sql(inference_detail$areas_sql_condition, results_inference)

  if(nrow(results_inference) == 0)
    return(data.frame())

  entrez_ids <- rep(NA, length(results_inference$AREA_OF_TEST))
  tryCatch({
    entrez_ids <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = results_inference$AREA_OF_TEST, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  }, error = function(e) {

  })


  results_inference$ENTREZID <- entrez_ids
  #
  entrez_ids_not_na <- entrez_ids[!is.na(entrez_ids)]
  if ((length(entrez_ids_not_na) != length(results_inference$AREA_OF_TEST)))
  {
    lost_gene <- unique(results_inference$AREA_OF_TEST[is.na(entrez_ids)])
  }

  results_inference <- assoc_filter_sql(inference_detail$association_results_sql_condition, results_inference)

  core_log_event("DEBUG: ",format(Sys.time(), "%a %b %d %X %Y")," inference file loaded:", inferenceFile)
  return(results_inference)
}
