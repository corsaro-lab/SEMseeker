#' Check a selection against the space of a meta analysis and apply it
#'
#' The space says what can be pooled; the selection says what the caller wants
#' pooled. The selection has to fall inside the space, and the interesting part
#' is what "fall inside" means for four coordinates that multiply.
#'
#' The check is \strong{per coordinate}, not on the product. A selection of two
#' markers and two areas names four combinations, and a cross product almost
#' always has an empty cell for reasons of biology rather than of error:
#' refusing the whole request over one of them would make the argument
#' unusable. A value that no study measured is a different thing, a mistake in
#' the request, so it is refused and named.
#'
#' What survives is the product intersected with the space, and the number of
#' combinations that fell outside is returned so the caller can say it.
#'
#' @param space The intersection from [meta_key_space()].
#' @param markers,figures,areas,subareas Optional character vectors.
#'   \code{NULL} leaves that coordinate unrestricted.
#'
#' @return A list with \code{keys}, the region classes to analyse;
#'   \code{requested}, how many combinations the selection named; and
#'   \code{dropped}, how many of them the space does not contain.
#'
#' @keywords internal
#' @noRd
meta_keys_select <- function(space, markers = NULL, figures = NULL,
                             areas = NULL, subareas = NULL)
{
  if (!is.data.frame(space) || nrow(space) == 0)
    stop("space must be a non-empty table of region classes.")
  missing_cols <- setdiff(.META_KEY_COLS, colnames(space))
  if (length(missing_cols) > 0)
    stop("space is missing the columns ", paste(missing_cols, collapse = ", "), ".")

  selection <- list(MARKER = markers, FIGURE = figures,
                    AREA = areas, SUBAREA = subareas)

  for (coord in names(selection)) {
    wanted <- selection[[coord]]
    if (is.null(wanted)) next
    if (!is.character(wanted) || length(wanted) == 0)
      stop(tolower(coord), " must be a non-empty character vector, or NULL.")
    available <- unique(space[[coord]])
    unknown <- setdiff(wanted, available)
    if (length(unknown) > 0)
      stop("no study measured ", tolower(coord), " ", paste(unknown, collapse = ", "),
           ". Available in the space of this meta analysis: ",
           paste(sort(available), collapse = ", "), ".")
  }

  keys <- space
  for (coord in names(selection)) {
    wanted <- selection[[coord]]
    if (!is.null(wanted)) keys <- keys[keys[[coord]] %in% wanted, , drop = FALSE]
  }
  rownames(keys) <- NULL

  # What the caller named, counted only over the coordinates they actually
  # named. An unrestricted coordinate is not a request for its product: it is a
  # request for whatever the space holds, so it cannot contribute a combination
  # that the space fails to contain. Counting it as a product reported classes
  # as dropped when nothing had been asked for.
  given <- selection[!vapply(selection, is.null, logical(1))]
  if (length(given) == 0) {
    requested <- nrow(keys)
    dropped <- 0L
  } else {
    named <- expand.grid(given, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    requested <- nrow(named)
    present <- vapply(seq_len(nrow(named)), function(i) {
      hit <- rep(TRUE, nrow(space))
      for (coord in names(given)) hit <- hit & space[[coord]] == named[i, coord]
      any(hit)
    }, logical(1))
    dropped <- sum(!present)
  }


  if (nrow(keys) == 0)
    stop("the selection names ", requested, " combinations and the space contains ",
         "none of them. Each value exists on its own, but no combination of them ",
         "was measured by every study.")

  list(keys = keys, requested = requested, dropped = dropped)
}
