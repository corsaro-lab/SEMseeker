#' Build an artefact from the position pivot, on demand (internal)
#'
#' AI-255. The single derivation path of the taxonomy:
#'
#' \preformatted{
#' POSITION/PROBE  ->  mask (AREA, SUBAREA)  ->  group  ->  aggregate
#'                                                  |
#'                                  SCOPE = INSTANCE -> one group per instance
#'                                  SCOPE = SAMPLE   -> one group
#' }
#'
#' **No aggregate is ever derived from another aggregate.** Not for elegance:
#' the mean of means is weighted wrong when regions hold different numbers of
#' positions, the median of medians does not exist, and — the case that makes the
#' rule general rather than a precaution — not even the sum of sums, because the
#' partition into genes is **not disjoint**. The annotation maps one probe onto
#' several genes, so adding the per-gene rows counts that probe once per gene.
#'
#' That is also why the two scopes part company here, and only here:
#'
#' * `INSTANCE` **explodes** the multi-gene areas, because a probe belonging to
#'   three genes must contribute to all three — they are three different
#'   questions (AI-050);
#' * `SAMPLE` **de-duplicates positions** instead: the mask says which positions
#'   belong to the region class, each position enters the single group once.
#'
#' Same grammar, same function, one branch — which is exactly the difference
#' between *selecting* and *selecting and grouping*.
#'
#' **The rule forbids composition, not co-computation.** One `group_by` emits
#' `SUM`, `MEAN`, `MEDIAN`, `VARIANCE` and `IQR` in a single pass, so a set of
#' aggregations costs one scan of the position pivot, not one scan each.
#'
#' @param marker,figure the quantity.
#' @param scope `"SAMPLE"` or `"INSTANCE"`.
#' @param area,subarea the region class.
#' @param aggregations character vector; defaults to everything admissible for
#'   the artefact. The algebraic ones are computed lazily in polars; the two
#'   modes are handled column by column in R (see [util_signal_descriptors()]),
#'   and are admissible only at `SCOPE = SAMPLE`.
#' @return named character vector of written paths, one per aggregation.
#' @keywords internal
#' @noRd
io_pivot_build <- function(marker, figure, scope, area, subarea,
                           aggregations = NULL) {

  scope <- io_scope_validate(scope)
  ssEnv <- core_get_session_info()

  if (is.null(aggregations))
    aggregations <- util_aggregations_allowed(marker, figure, default = FALSE,
                                              scope = scope, area = area)
  aggregations <- unique(core_name_cleaning(as.character(aggregations)))

  illegal <- setdiff(aggregations,
                     util_aggregations_allowed(marker, figure, default = FALSE,
                                               scope = scope, area = area))
  if (length(illegal) > 0)
    stop("io_pivot_build(): aggregation(s) ", paste(illegal, collapse = ", "),
         " are not admissible for ", marker, "/", figure, " at scope ", scope,
         " on ", area, "_", subarea, ". Admissible: ",
         paste(util_aggregations_allowed(marker, figure, default = FALSE,
                                         scope = scope, area = area),
               collapse = ", "), ".", call. = FALSE)

  grouped <- .io_pivot_masked_lazy(marker, figure, scope, area, subarea)
  if (is.null(grouped))
    return(character(0))

  key_col     <- grouped$key_col
  sample_cols <- setdiff(names(grouped$lazy), key_col)
  if (length(sample_cols) == 0)
    return(character(0))

  written <- character(0)

  algebraic <- setdiff(aggregations, c("MODELOW", "MODEHIGH"))
  if (length(algebraic) > 0) {
    # One pass, every algebraic aggregation at once. Composition is forbidden;
    # co-computation is what makes the rule affordable.
    for (agg in algebraic) {
      exprs <- lapply(sample_cols, function(cl)
        .io_agg_expr(cl, agg)$alias(cl))
      # polars' R bindings take expressions through `...`, not as a list, so the
      # list has to be splatted with do.call(). Same idiom as sem_deltaX_get().
      out <- if (identical(scope, "SAMPLE"))
        do.call(grouped$lazy$select, exprs)$with_columns(
          polars::pl$lit(key_col_value(agg))$alias("AREA"))
      else
        do.call(grouped$lazy$group_by(key_col, .maintain_order = FALSE)$agg, exprs)

      path <- io_pivot_file_name_parquet(marker, figure, area, subarea,
                                         aggregation = agg, scope = scope)
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      out$collect()$write_parquet(path)
      written[agg] <- path
      core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
                     " io_pivot_build wrote ", basename(path))
    }
  }

  modes <- intersect(aggregations, c("MODELOW", "MODEHIGH"))
  if (length(modes) > 0)
    written <- c(written,
                 .io_pivot_build_modes(grouped, marker, figure, scope, area,
                                       subarea, modes, sample_cols))

  written
}

#' The key value of a collapsed artefact (internal)
#'
#' AI-255. A `SCOPE = SAMPLE` artefact has exactly one row, and the key of that
#' row becomes `AREA_OF_TEST` downstream — the answer to "what was tested".
#'
#' For a collapsed artefact that answer is **which aggregate**, not which region
#' class: the class is already said by the `AREA` and `SUBAREA` columns, and
#' repeating it would leave the mean and the median of the same class carrying
#' the identical `AREA_OF_TEST`, distinguishable only by a column further right.
#'
#' @param aggregation the operator that produced the single row.
#' @keywords internal
#' @noRd
key_col_value <- function(aggregation) {
  core_name_cleaning(as.character(aggregation))
}

#' Polars expression for one aggregation on one column (internal)
#'
#' AI-255. `IQR` has no direct polars reduction and is the difference of two
#' quantiles — which is its definition, not a workaround.
#'
#' @keywords internal
#' @noRd
.io_agg_expr <- function(column, aggregation) {
  cl <- polars::pl$col(column)
  switch(
    aggregation,
    SUM      = cl$sum(),
    MEAN     = cl$mean(),
    MEDIAN   = cl$median(),
    VARIANCE = cl$var(),
    IQR      = (cl$quantile(0.75) - cl$quantile(0.25)),
    VALUE    = cl$first(),
    stop(".io_agg_expr(): '", aggregation, "' has no lazy form; it must be ",
         "handled column by column in R.", call. = FALSE)
  )
}

#' Mask the position pivot onto a region class (internal)
#'
#' Returns `list(lazy = <LazyFrame>, key_col = <character>)`, or `NULL` when the
#' position pivot is not available.
#'
#' @keywords internal
#' @noRd
.io_pivot_masked_lazy <- function(marker, figure, scope, area, subarea) {

  base <- io_read_pivot(marker, figure, area = "POSITION", subarea = "WHOLE",
                        aggregation = "VALUE", scope = "INSTANCE", build = FALSE)
  if (is.null(base))
    return(NULL)

  base <- base$with_columns(
    polars::pl$col("START")$cast(polars::pl$Int32),
    polars::pl$col("END")$cast(polars::pl$Int32),
    polars::pl$col("CHR")$cast(polars::pl$String)$str$replace("^(?i)chr", "")
  )

  # The whole genome, unmasked: no annotation needed, every position takes part.
  if (io_area_is_single_position(area) && identical(toupper(subarea), "WHOLE")) {
    if (identical(scope, "SAMPLE"))
      return(list(lazy = base$drop(intersect(c("CHR", "START", "END", "PROBE"),
                                             names(base))),
                  key_col = "AREA"))

    # A position pivot is keyed by CHR/START/END — three columns — while every
    # consumer downstream wants ONE identifier per row, because that identifier
    # becomes AREA_OF_TEST. Compose it the way the package already does:
    # `<CHR>_<START>`, the synthetic probe id io_probe_id_to_coord() splits back
    # on the last underscore (START is always a pure integer suffix).
    #
    # Without this the three coordinates would fall into `sample_cols` and be
    # aggregated as if they were samples.
    if (!("AREA" %in% names(base)) && all(c("CHR", "START") %in% names(base)))
      base <- base$with_columns(
        (polars::pl$col("CHR")$cast(polars::pl$String) + polars::pl$lit("_") +
           polars::pl$col("START")$cast(polars::pl$String))$alias("AREA")
      )
    base <- base$drop(intersect(c("CHR", "START", "END", "PROBE"), names(base)))
    return(list(lazy = base, key_col = "AREA"))
  }

  area_subarea <- paste0(area, "_", ifelse(subarea == "", "WHOLE", subarea))
  probe_features <- anno_probe_features_get(area_subarea)
  probe_features$CHR <- as.character(probe_features$CHR)
  pf <- polars::as_polars_df(probe_features)$lazy()
  pf <- pf$with_columns(polars::pl$col(area_subarea)$alias("AREA"))$drop(area_subarea)
  pf <- pf$with_columns(
    polars::pl$col("START")$cast(polars::pl$Int32),
    polars::pl$col("END")$cast(polars::pl$Int32),
    polars::pl$col("CHR")$cast(polars::pl$String)$str$replace("^(?i)chr", "")
  )

  if (identical(scope, "SAMPLE")) {
    # SELECT ONLY. The positions of the class, each one once: a probe annotated
    # to three genes is one position of the sample, not three.
    #
    # AI-308: `drop_nulls("AREA")` first, and it is not decoration. On the
    # Illumina path anno_probe_features_get() returns the WHOLE annotation table
    # (every probe of the array) with NA in the column of the class asked for:
    # for K450, GENE_TSS1500 carries 84,808 annotated probes against 401,394 NA.
    # Selecting the coordinates without dropping those rows built the mask on
    # all 485,512 positions whatever the class, so the inner join below excluded
    # nothing and every "restricted" burden came out equal to the burden of the
    # whole sample: identical for GENE_TSS1500, ISLAND_N_SHORE and GENE_WHOLE
    # alike. Silent, because a number was produced for each class.
    #
    # The INSTANCE branch below never had the defect: it drops the same nulls
    # after the join. The two branches part company here, so they have to select
    # the same positions here.
    keys <- pf$drop_nulls("AREA")$select(c("CHR", "START", "END"))$unique()
    masked <- base$join(keys, on = c("CHR", "START", "END"), how = "inner")
    drop_cols <- intersect(c("CHR", "START", "END", "PROBE", "K27", "K450", "K850"),
                           names(masked))
    return(list(lazy = masked$drop(drop_cols), key_col = "AREA"))
  }

  # SELECT AND GROUP. Here the explode is required: the probe must contribute to
  # each of its genes, because those are different questions.
  joined <- pf$join(base, on = c("CHR", "START", "END"), how = "inner")
  drop_cols <- intersect(c("PROBE", "CHR", "START", "END", "K27", "K450", "K850"),
                         names(joined))
  joined <- joined$drop(drop_cols)$drop_nulls("AREA")
  joined <- .anno_area_explode(joined)
  list(lazy = joined, key_col = "AREA")
}

#' The two modes, column by column (internal)
#'
#' AI-255. They have no lazy form: the estimate needs the whole distribution of
#' the group. Admissible only at `SCOPE = SAMPLE` — see
#' [util_aggregations_allowed()] — where a group is one sample, so this walks one
#' column at a time and never materialises the matrix.
#'
#' @keywords internal
#' @noRd
.io_pivot_build_modes <- function(grouped, marker, figure, scope, area, subarea,
                                  modes, sample_cols) {

  if (!identical(scope, "SAMPLE"))
    stop("MODELOW/MODEHIGH are admissible only at scope SAMPLE: per instance a ",
         "region holds a handful of positions (about nineteen probes in a gene, ",
         "two to five in a TSS200 window) and a two-peak density estimate on ",
         "those is the bandwidth's opinion, not the data's.", call. = FALSE)

  values <- lapply(stats::setNames(modes, modes), function(x) numeric(0))
  for (cl in sample_cols) {
    v <- as.data.frame(grouped$lazy$select(polars::pl$col(cl))$collect())[[1]]
    d <- util_signal_descriptors(v[is.finite(v)], beta = TRUE)
    for (m in modes)
      values[[m]] <- c(values[[m]], if (is.null(d[[m]])) NA_real_ else d[[m]])
  }

  written <- character(0)
  for (m in modes) {
    df <- as.data.frame(as.list(stats::setNames(values[[m]], sample_cols)),
                        check.names = FALSE)
    df <- cbind(AREA = key_col_value(m), df, stringsAsFactors = FALSE)
    path <- io_pivot_file_name_parquet(marker, figure, area, subarea,
                                       aggregation = m, scope = scope)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    polars::as_polars_df(df)$write_parquet(path)
    written[m] <- path
  }
  written
}
