#' Build an artefact from the position pivot, on demand (internal)
#'
#' The single derivation path of the taxonomy:
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
#' positions, the median of medians does not exist, and - the case that makes the
#' rule general rather than a precaution - not even the sum of sums, because the
#' partition into genes is **not disjoint**. The annotation maps one probe onto
#' several genes, so adding the per-gene rows counts that probe once per gene.
#'
#' That is also why the two scopes part company here, and only here:
#'
#' * `INSTANCE` **explodes** the multi-gene areas, because a probe belonging to
#'   three genes must contribute to all three - they are three different
#'   questions;
#' * `SAMPLE` **de-duplicates positions** instead: the mask says which positions
#'   belong to the region class, each position enters the single group once.
#'
#' Same grammar, same function, one branch - which is exactly the difference
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
#'   the artefact. The algebraic ones are computed by the data engine, which
#'   copies the result to parquet without bringing it into R; the two modes are
#'   handled column by column in R (see [util_signal_descriptors()]), and are
#'   admissible only at `SCOPE = SAMPLE`.
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
  sample_cols <- setdiff(colnames(grouped$lazy), key_col)
  if (length(sample_cols) == 0)
    return(character(0))

  connection   <- dbplyr::remote_con(grouped$lazy)
  masked_query <- dbplyr::sql_render(grouped$lazy)
  column_types <- .io_query_column_types(connection, masked_query)
  written <- character(0)

  algebraic <- setdiff(aggregations, c("MODELOW", "MODEHIGH"))
  for (agg in algebraic) {
    aggregate_list <- paste(vapply(sample_cols, function(cl)
      paste(.io_agg_sql(connection, cl, agg, column_types[[cl]]), "AS", DBI::dbQuoteIdentifier(connection, cl)),
      character(1)), collapse = ",\n  ")
    key_identifier <- DBI::dbQuoteIdentifier(connection, key_col)
    # SAMPLE: one row, the key last (the column the polars build appended).
    # INSTANCE: one row per instance, the key first; sorted, so the file is the
    # same at every run.
    query_text <- if (identical(scope, "SAMPLE"))
      sprintf("SELECT %s,\n  %s AS %s\nFROM (%s) AS masked",
              aggregate_list, DBI::dbQuoteString(connection, key_col_value(agg)), key_identifier, masked_query)
    else
      sprintf("SELECT %s,\n  %s\nFROM (%s) AS masked\nGROUP BY %s\nORDER BY %s",
              key_identifier, aggregate_list, masked_query, key_identifier, key_identifier)

    path <- io_pivot_file_name_parquet(marker, figure, area, subarea,
                                       aggregation = agg, scope = scope)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    DBI::dbExecute(connection, sprintf("COPY (%s) TO %s (FORMAT parquet, COMPRESSION zstd)",
                                       query_text, DBI::dbQuoteString(connection, path)))
    written[agg] <- path
    core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
                   " io_pivot_build wrote ", basename(path))
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
#' A `SCOPE = SAMPLE` artefact has exactly one row, and the key of that
#' row becomes `AREA_OF_TEST` downstream - the answer to "what was tested".
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

#' Column types of a query, as the data engine reports them (internal)
#'
#' The SQL of an aggregation depends on the type of the column it reduces: a
#' sum of integers is cast to BIGINT, which is what the polars build wrote,
#' where DuckDB would widen it to a 128-bit integer.
#'
#' @param connection the data engine connection the query runs on
#' @param query_text the SQL of the query
#' @return named character vector, column name to DuckDB type
#' @keywords internal
#' @noRd
.io_query_column_types <- function(connection, query_text) {
  described <- DBI::dbGetQuery(connection, paste("DESCRIBE", query_text))
  stats::setNames(described$column_type, described$column_name)
}

#' SQL for one aggregation on one column (internal)
#'
#' Each operator reproduces what the polars build computed, so that moving the
#' engine moves no number:
#' * `SUM` of a column with no value is 0, not NULL;
#' * `MEAN` is `favg`, a compensated sum: polars sums with compensation too,
#'   and DuckDB's plain `avg` drifted by up to 6.6e-15 relative on 18 000
#'   values, where `favg` stays at the last bit;
#' * `VARIANCE` is the sample variance (n - 1). It differs from polars by up
#'   to 3.5e-15 relative, a known limit of the engine change: no DuckDB
#'   formula reproduced polars bit for bit;
#' * `MEDIAN` interpolates between the two central values;
#' * `IQR` is the difference of the 0.75 and 0.25 quantiles taken as the
#'   element of rank `floor((n - 1) q + 0.5)` in the sorted non-missing values,
#'   which is polars' "nearest" rule. DuckDB's `quantile_disc` uses another rank
#'   and differed in half of the cases measured;
#' * `VALUE` is the first value of the group: it is admissible only where a
#'   group holds one position.
#'
#' @param connection the data engine connection, for quoting
#' @param column the column to reduce
#' @param aggregation one of SUM, MEAN, MEDIAN, VARIANCE, IQR, VALUE
#' @param column_type the DuckDB type of the column, from
#'   .io_query_column_types(); decides the cast of a SUM
#' @return the SQL expression, as a character string
#' @keywords internal
#' @noRd
.io_agg_sql <- function(connection, column, aggregation, column_type) {
  quoted <- DBI::dbQuoteIdentifier(connection, column)
  is_integer <- toupper(column_type) %in% c("TINYINT", "SMALLINT", "INTEGER", "BIGINT", "UTINYINT", "USMALLINT", "UINTEGER")
  # SQL for polars' "nearest" quantile, which no DuckDB function computes:
  # the element of rank floor((n - 1) q + 0.5) among the sorted non-missing
  # values. quantile_level: q, a number in [0, 1].
  nearest <- function(quantile_level)
    sprintf("list_extract(list_sort(list(%1$s) FILTER (WHERE %1$s IS NOT NULL)), CAST(floor((count(%1$s) - 1) * %2$s + 0.5) AS BIGINT) + 1)",
            quoted, quantile_level)
  switch(
    aggregation,
    SUM      = if (is_integer) sprintf("CAST(COALESCE(sum(%s), 0) AS BIGINT)", quoted)
               else sprintf("COALESCE(sum(%s), 0)", quoted),
    MEAN     = sprintf("favg(%s)", quoted),
    MEDIAN   = sprintf("median(%s)", quoted),
    VARIANCE = sprintf("var_samp(%s)", quoted),
    IQR      = sprintf("CASE WHEN count(%s) = 0 THEN NULL ELSE %s - %s END", quoted, nearest(0.75), nearest(0.25)),
    VALUE    = sprintf("first(%s)", quoted),
    stop(".io_agg_sql(): '", aggregation, "' has no SQL form; it must be ",
         "handled column by column in R.", call. = FALSE)
  )
}

#' Mask the position pivot onto a region class (internal)
#'
#' Returns `list(lazy = <tbl_lazy>, key_col = <character>)`, or `NULL` when the
#' position pivot is not available. The lazy table holds the key column and
#' one column per sample; nothing is read until it is collected or aggregated.
#'
#' @param marker,figure the quantity.
#' @param scope `"SAMPLE"` or `"INSTANCE"`.
#' @param area,subarea the region class.
#' @keywords internal
#' @noRd
.io_pivot_masked_lazy <- function(marker, figure, scope, area, subarea) {

  base_path <- io_pivot_ensure(marker, figure, area = "POSITION", subarea = "WHOLE",
                               aggregation = "VALUE", scope = "INSTANCE", build = FALSE)
  if (is.null(base_path))
    return(NULL)

  connection <- io_engine_connection()
  # Sample ids are free text: every column name enters the SQL quoted.
  # name: a column name; returns it as a quoted SQL identifier.
  quote_name <- function(name) as.character(DBI::dbQuoteIdentifier(connection, name))
  base_columns <- colnames(io_scan_parquet(base_path))
  # Coordinates normalised as the polars build did: 32-bit START and END, CHR
  # as text without a leading "chr". The join with the annotation needs both
  # sides in the same form.
  base_select <- vapply(base_columns, function(cl) switch(cl,
    START = 'CAST(START AS INTEGER) AS START',
    END   = 'CAST("END" AS INTEGER) AS "END"',
    CHR   = "regexp_replace(CAST(CHR AS VARCHAR), '^(?i)chr', '') AS CHR",
    quote_name(cl)), character(1))
  base_query <- sprintf("SELECT %s FROM read_parquet(%s)", paste(base_select, collapse = ", "),
                        DBI::dbQuoteString(connection, normalizePath(base_path, winslash = "/")))
  # Every column but these is kept, in its order: what is left is the key and
  # the samples.
  # columns: the column names, in order; dropped: the names to leave out.
  without <- function(columns, dropped) setdiff(columns, dropped)
  # The callers aggregate with dplyr verbs, so the masked query is handed
  # back as a lazy table. query_text: the SELECT that masks the pivot.
  as_lazy <- function(query_text) dplyr::tbl(connection, dbplyr::sql(query_text))

  # The whole genome, unmasked: no annotation needed, every position takes part.
  if (io_area_is_single_position(area) && identical(toupper(subarea), "WHOLE")) {
    coordinates <- c("CHR", "START", "END", "PROBE")
    if (identical(scope, "SAMPLE"))
      return(list(lazy = as_lazy(sprintf("SELECT %s FROM (%s) AS base",
                                         paste(quote_name(without(base_columns, coordinates)), collapse = ", "), base_query)),
                  key_col = "AREA"))

    # A position pivot is keyed by CHR/START/END - three columns - while every
    # consumer downstream wants ONE identifier per row, because that identifier
    # becomes AREA_OF_TEST. Compose it the way the package already does:
    # `<CHR>_<START>`, the synthetic probe id io_probe_id_to_coord() splits back
    # on the last underscore (START is always a pure integer suffix).
    #
    # Without this the three coordinates would fall into `sample_cols` and be
    # aggregated as if they were samples.
    kept <- quote_name(without(base_columns, coordinates))
    if (!("AREA" %in% base_columns) && all(c("CHR", "START") %in% base_columns))
      kept <- c(kept, "CHR || '_' || CAST(START AS VARCHAR) AS AREA")
    return(list(lazy = as_lazy(sprintf("SELECT %s FROM (%s) AS base", paste(kept, collapse = ", "), base_query)),
                key_col = "AREA"))
  }

  area_subarea <- paste0(area, "_", ifelse(subarea == "", "WHOLE", subarea))
  probe_features <- anno_probe_features_get(area_subarea)
  probe_features$CHR <- as.character(probe_features$CHR)
  # The class column becomes AREA, last, as the polars build renamed it.
  probe_features$AREA <- probe_features[[area_subarea]]
  probe_features[[area_subarea]] <- NULL
  features_query <- dbplyr::sql_render(io_register_frame(probe_features))
  feature_columns <- colnames(probe_features)
  feature_select <- vapply(feature_columns, function(cl) switch(cl,
    START = 'CAST(START AS INTEGER) AS START',
    END   = 'CAST("END" AS INTEGER) AS "END"',
    CHR   = "regexp_replace(CHR, '^(?i)chr', '') AS CHR",
    quote_name(cl)), character(1))
  features_query <- sprintf("SELECT %s FROM (%s) AS features", paste(feature_select, collapse = ", "), features_query)
  annotation_columns <- c("PROBE", "CHR", "START", "END", "K27", "K450", "K850")

  if (identical(scope, "SAMPLE")) {
    # SELECT ONLY. The positions of the class, each one once: a probe annotated
    # to three genes is one position of the sample, not three.
    #
    # The NULL areas are dropped first, and it is not decoration. On the
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
    kept <- quote_name(without(base_columns, annotation_columns))
    masked_query <- sprintf(
      "SELECT %s FROM (%s) AS base\nSEMI JOIN (SELECT DISTINCT CHR, START, \"END\" FROM (%s) AS features WHERE AREA IS NOT NULL) AS keys\nUSING (CHR, START, \"END\")",
      paste(paste0("base.", kept), collapse = ", "), base_query, features_query)
    return(list(lazy = as_lazy(masked_query), key_col = "AREA"))
  }

  # SELECT AND GROUP. Here the explode is required: the probe must contribute to
  # each of its genes, because those are different questions. The columns keep
  # the order of the polars join: the annotation's first, then the pivot's.
  kept <- c(paste0("features.", quote_name(without(feature_columns, annotation_columns))),
            paste0("base.", quote_name(without(base_columns, annotation_columns))))
  joined_query <- sprintf(
    "SELECT %s FROM (%s) AS features\nJOIN (%s) AS base USING (CHR, START, \"END\")\nWHERE features.AREA IS NOT NULL",
    paste(kept, collapse = ", "), features_query, base_query)
  list(lazy = as_lazy(.anno_area_explode(connection, joined_query)), key_col = "AREA")
}

#' The two modes, column by column (internal)
#'
#' They have no lazy form: the estimate needs the whole distribution of
#' the group. Admissible only at `SCOPE = SAMPLE` - see
#' [util_aggregations_allowed()] - where a group is one sample, so this walks one
#' column at a time and never materialises the matrix.
#'
#' @param grouped the masked pivot, as returned by .io_pivot_masked_lazy()
#' @param marker,figure,scope,area,subarea the artefact, for its file name
#' @param modes the modes to compute, a subset of MODELOW and MODEHIGH
#' @param sample_cols the sample columns of the masked pivot
#' @return named character vector of written paths, one per mode
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
    v <- as.data.frame(dplyr::collect(dplyr::select(grouped$lazy, dplyr::all_of(cl))))[[1]]
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
    io_sink_parquet(io_register_frame(df), path)
    written[m] <- path
  }
  written
}
