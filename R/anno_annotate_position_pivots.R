#' @title Materialise the area-level artefacts of the run
#'
#' Rewritten as a thin batch pass over [io_pivot_build()], which is now
#' the single derivation path of the taxonomy. What used to live here - the join
#' with the annotation, the multi-gene explode, the group-by - moved into
#' `io_pivot_build()` and [.anno_area_explode()], so the batch pre-pass and the
#' on-demand construction cannot diverge.
#'
#' What disappeared with it is the rule that chose the operator:
#'
#' \preformatted{
#' if (localKeys[i, "DISCRETE"]) pivot$group_by("AREA")$sum()
#' else                          pivot$group_by("AREA")$mean()
#' }
#'
#' One file per key, the operator picked by a flag and **absent from the name**.
#' Downstream, `inference_details$aggregation` was validated and then ignored:
#' asking for `MEDIAN` on `GENE_TSS1500` returned the mean, silently. The keys
#' now carry the aggregation, so what is written says which operator produced it.
#'
#' This pass is a convenience, not a requirement: anything it does not
#' materialise is built on first read by [io_read_pivot()].
#'
#' @return nothing
#' @importFrom doRNG %dorng%
anno_annotate_position_pivots <- function ()
{
  start_time <- Sys.time()
  ssEnv <- core_get_session_info()
  localKeys <- ssEnv$keys_areas_subareas_markers_figures

  # POSITION is the source, not a destination.
  localKeys <- localKeys[localKeys$AREA != "POSITION", ]

  if (nrow(localKeys) == 0)
    return()

  # Short-circuit: if every destination artefact already exists there is nothing
  # to annotate. Avoids the unconditional anno_probe_features_get() load of the
  # Illumina manifest (~10-30s) and the spurious log line on resume.
  wanted <- lapply(seq_len(nrow(localKeys)), function(i)
    util_aggregations_allowed(localKeys[i, "MARKER"], localKeys[i, "FIGURE"],
                              discrete = isTRUE(localKeys[i, "DISCRETE"]),
                              default  = TRUE,
                              scope    = "INSTANCE",
                              area     = localKeys[i, "AREA"]))

  all_dest_exist <- all(vapply(seq_len(nrow(localKeys)), function(i)
    all(vapply(wanted[[i]], function(agg)
      file.exists(io_pivot_file_name_parquet(
        localKeys[i, "MARKER"], localKeys[i, "FIGURE"],
        localKeys[i, "AREA"],   localKeys[i, "SUBAREA"],
        aggregation = agg, scope = "INSTANCE")), logical(1))),
    logical(1)))

  if (all_dest_exist) {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
      " Annotation skipped: all destination artefacts already exist.")
    return()
  }

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Annotating genomic area.")

  progress_bar <- ""
  if (ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = seq_len(nrow(localKeys)))

  for (i in seq_len(nrow(localKeys)))
  {
    marker  <- as.character(localKeys[i, "MARKER"])
    figure  <- as.character(localKeys[i, "FIGURE"])
    subarea <- as.character(localKeys[i, "SUBAREA"])
    area    <- as.character(localKeys[i, "AREA"])

    missing_aggs <- Filter(function(agg)
      !file.exists(io_pivot_file_name_parquet(marker, figure, area, subarea,
                                              aggregation = agg,
                                              scope = "INSTANCE")),
      wanted[[i]])

    if (length(missing_aggs) > 0) {
      written <- tryCatch(
        io_pivot_build(marker, figure, scope = "INSTANCE", area = area,
                       subarea = subarea, aggregations = missing_aggs),
        error = function(e) {
          core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
            " Annotation of ", marker, "_", figure, " on ", area, "_", subarea,
            " failed: ", conditionMessage(e))
          character(0)
        })

      if (length(written) == 0)
        ssEnv$key_missed_areas_subareas <- unique(rbind(
          ssEnv$key_missed_areas_subareas, localKeys[i, c("AREA", "SUBAREA")]))
    }

    if (ssEnv$showprogress)
      progress_bar(sprintf("Annotating position pivots."))
  }

  # remove missed keys
  selector <- !((ssEnv$keys_areas_subareas_markers_figures$AREA %in% ssEnv$key_missed_areas_subareas$AREA) & (ssEnv$keys_areas_subareas_markers_figures$SUBAREA %in% ssEnv$key_missed_areas_subareas$SUBAREA))
  ssEnv$keys_areas_subareas_markers_figures  <- ssEnv$keys_areas_subareas_markers_figures[selector,]

  core_update_session_info(ssEnv)
  end_time <- Sys.time()
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Annotation genomic areas file finished in ", difftime(end_time,start_time,units = "mins")," minutes.")
  gc()
}
