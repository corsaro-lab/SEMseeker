#' Compose the on-disk name of a pivot (internal)
#'
#' AI-248, rewritten by AI-255. The file name of a pivot **is** its identity key
#' ([io_artefact_key()]) plus an extension: a pivot copied out of its folder and
#' attached to a mail still says what it is.
#'
#' \preformatted{
#' <MARKER>_<FIGURE>_<SCOPE>_<AREA>_<SUBAREA>_<AGGREGATION>_<build>.parquet
#'
#' MUTATIONS_HYPER_SAMPLE_PROBE_WHOLE_SUM_HG19.parquet
#' DELTAS_HYPO_INSTANCE_GENE_TSS200_MEDIAN_HG19.parquet
#' SIGNAL_BETA_INSTANCE_PROBE_WHOLE_VALUE_HG19.parquet
#' }
#'
#' The aggregation is in the name because before AI-248 the name did **not** say
#' which operator had produced the file: an existing pivot was reused on trust.
#' With one operator per marker that trust was justified; with several it is the
#' worst kind of bug - silent, because the file is there and the run does not
#' complain while a mean is consumed as if it were a sum.
#'
#' The scope is in the name because AI-255 made it a coordinate: a burden over
#' the whole sample and a burden per gene are the same marker reduced over
#' different extents, and used to live in files of different *shape* - a CSV
#' sibling and a pivot - which is what hid the difference.
#'
#' @param marker,figure,area,subarea the region class and the quantity.
#' @param aggregation the AGGREGATION axis. `NULL` is resolved to `"VALUE"` when
#'   the block is already a single position (`PROBE`/`POSITION` at
#'   `SCOPE = INSTANCE`), where `VALUE` is the only meaningful operator, and is
#'   an error anywhere else: a block holding many positions was reduced somehow,
#'   and the name has to say how.
#' @param scope `"INSTANCE"` (default, one row per instance) or `"SAMPLE"` (one
#'   row, the whole extent collapsed).
#' @param add_gz passed through for the CSV flavour.
#' @keywords internal
#' @noRd
io_pivot_file_name <- function(marker, figure, area, subarea, add_gz = TRUE,
                               aggregation = NULL, scope = "INSTANCE")
{
  .io_pivot_path(marker, figure, area, subarea, aggregation, scope, ".csv", add_gz)
}

#' @rdname io_pivot_file_name
#' @keywords internal
#' @noRd
io_pivot_file_name_parquet <- function(marker, figure, area, subarea,
                                       aggregation = NULL, scope = "INSTANCE")
{
  .io_pivot_path(marker, figure, area, subarea, aggregation, scope, ".parquet", FALSE)
}

#' Resolve an omitted aggregation, or refuse to guess (internal)
#'
#' AI-255. The absence of the segment used to mean "POSITION-level pivot, nothing
#' aggregated yet". That reading only ever worked because there was one such
#' case; with `VALUE` in the vocabulary the case has a name, and the absence goes
#' back to being an error everywhere else.
#'
#' @keywords internal
#' @noRd
.io_aggregation_resolve <- function(aggregation, scope, area) {
  if (!is.null(aggregation) && nzchar(as.character(aggregation)[1]))
    return(core_name_cleaning(as.character(aggregation)[1]))

  if (identical(scope, "INSTANCE") && io_area_is_single_position(area))
    return("VALUE")

  stop("aggregation is required for scope '", scope, "' on area '", area,
       "': a block holding more than one position was reduced somehow, and the ",
       "name has to say how. Legal names: ",
       paste(util_aggregation_vocabulary(), collapse = ", "), ".", call. = FALSE)
}

#' @keywords internal
#' @noRd
.io_pivot_path <- function(marker, figure, area, subarea, aggregation, scope,
                           extension, add_gz)
{
  ssEnv <- core_get_session_info()
  scope <- io_scope_validate(scope)

  # No normalisation of an empty subarea. There is nothing to normalise: the
  # registry assigns every area an explicit subarea, WHOLE at minimum
  # (util_keys_create(): PROBE/WHOLE, POSITION/WHOLE, DMR/WHOLE, ...), so "" is not
  # a value the system produces. Where it used to appear it came from a caller
  # writing it by hand, and turning it into WHOLE here would accept a hole in a
  # coordinate - the same thing as filling a missing SUBAREA with "TOTAL", which
  # is how a defect hides inside a key. io_artefact_key() refuses it instead, and
  # the message names the caller's mistake.
  aggregation <- .io_aggregation_resolve(aggregation, scope, area)

  reportFolder <- io_dir_check_and_create(ssEnv$result_folderData, "Pivots")
  pivot_subfolder <- io_dir_check_and_create(reportFolder, marker)

  pivot_base <- io_artefact_key(marker = marker, figure = figure, scope = scope,
                                area = area, subarea = subarea,
                                aggregation = aggregation)

  io_file_path_build(baseFolder = pivot_subfolder, detailsFilename = pivot_base,
                     extension = extension, add_gz = add_gz)
}
