io_get_pivot_both <- function(marker)
{
  ssEnv <- core_get_session_info()
  is_discrete <- unique(ssEnv$keys_markers_figures_default[ssEnv$keys_markers_figures_default$MARKER==marker,"DISCRETE"])

  # Migrated 2026-06-01 / 2026-06-03 to use the new io_read_pivot() dispatcher (AI-027).
  # Branch (a): BOTH already materialised (cache).
  pivot_both_lazy <- io_read_pivot(marker, "BOTH", area = "PROBE", subarea = "WHOLE")
  if (!is.null(pivot_both_lazy))
    return(as.data.frame(pivot_both_lazy$collect()))

  # Branch (b): build BOTH from HYPER+HYPO. io_read_pivot() transparently handles
  # both the cached-parquet and the bed-streaming-merge backends; if neither
  # the parquet nor the per-sample bed files exist for a figure, it returns
  # NULL and we treat that figure as empty.
  pivot_hyper <- io_read_pivot(marker, "HYPER", area = "PROBE", subarea = "WHOLE")
  pivot_hypo  <- io_read_pivot(marker, "HYPO",  area = "PROBE", subarea = "WHOLE")

  if (is.null(pivot_hyper) && is.null(pivot_hypo))
    return(data.frame())
  # Drop the NULL side(s); polars concat requires non-NULL operands
  parts <- Filter(Negate(is.null), list(pivot_hyper, pivot_hypo))

  # Union (concatenate) the available figures.
  # NOTA polars 1.11 (e già da 1.x): pl$concat() vuole i frame come varargs (...),
  # non come list. Usiamo do.call per espandere la lista in posizionali.
  # how="diagonal" allinea schemi diversi (HYPER e HYPO hanno set di sample
  # diversi: solo i sample con outlier in quella figura compaiono come colonna).
  # I sample mancanti nell'altra figura diventano null e il sum() successivo
  # li ignora correttamente.
  pivot_both <- do.call(polars::pl$concat, c(parts, list(how = "diagonal_relaxed")))$collect()
  pivot_both <- pivot_both$group_by("AREA", .maintain_order=FALSE)$sum()

  # The cache this branch exists to fill. The target path used to be a bare
  # `pivot_file_name_both`, a variable this function never assigned, so the
  # build-from-HYPER-and-HYPO branch raised `object not found` after doing all
  # of its work - and that branch is the normal one, taken whenever no BOTH
  # pivot is already cached. The name comes from io_pivot_file_name_parquet(),
  # the same function io_read_pivot() looks the cache up with, so what is
  # written here is what the fast path above will find next time.
  pivot_file_name_both <- io_pivot_file_name_parquet(marker, "BOTH",
                                                     "PROBE", "WHOLE")
  dir.create(dirname(pivot_file_name_both), recursive = TRUE, showWarnings = FALSE)
  pivot_both$write_parquet(pivot_file_name_both)
  # Sidecar JSON is materialised by core_ensure_sidecars() at pipeline end.

  # pivot_both <- pivot_both$group_by("AREA", maintain_order=FALSE)$agg(
  #   pl$all()$exclude("AREA")$sum()
  # )
  #
  # pivot_both$sink_parquet(pivot_file_name_both)


  # pivot_hyper <- data.frame()
  # if(file.exists(pivot_file_name_hyper))
  #   pivot_hyper <- polars::pl$read_parquet(pivot_file_name_hyper)$to_data_frame()
  # pivot_hypo <- data.frame()
  # if(file.exists(pivot_file_name_hypo))
  #   pivot_hypo <- polars::pl$read_parquet(pivot_file_name_hypo)$to_data_frame()
  #
  # # get the row with the max count of values
  # count_m <- plyr::rbind.fill(as.data.frame(pivot_hyper), as.data.frame(pivot_hypo))
  # rm(pivot_hypo)
  # rm(pivot_hyper)
  #
  #
  # # sort count_m by AREA
  # count_m <- count_m[order(count_m$AREA),]
  #
  # if (nrow(count_m)!=0)
  #   # Group by AREA and sum all other columns
  #   df_grouped <- df$group_by("AREA")$agg(
  #     pl$all()$exclude("AREA")$sum()
  #   )
  #
  # count_m <- as.data.frame(count_m)
  # arrow::write_parquet(count_m, pivot_file_name_both)
  return(as.data.frame(pivot_both))
}
