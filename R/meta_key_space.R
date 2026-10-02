# The space of a meta_ analysis: the region classes its studies have in common.
#
# Split in two on purpose. Reading a study's classes needs a session, so it is
# not testable without result folders; intersecting the tables is arithmetic and
# has to be, because the strictness of the intersection is a decision and a
# decision is worth a test.

.META_KEY_COLS <- c("MARKER", "FIGURE", "AREA", "SUBAREA")

#' Read the region classes each study declares
#'
#' No selection is applied here, and that is the point: the space has to be the
#' classes the studies actually share, not the classes the caller asked about.
#' Filtering the sessions first would compute the intersection over the narrowed
#' space, and a class missing from one study would vanish from the space instead
#' of falling outside it, leaving nothing to refuse.
#'
#' @param studies A table from [meta_studies_normalise()].
#' @param ... Passed to the per-study session setup.
#'
#' @return A named list of \code{data.frame}s, one per study, each with the four
#'   region-class columns.
#'
#' @keywords internal
#' @noRd
meta_key_sets_read <- function(studies, ...)
{
  key_sets <- lapply(seq_len(nrow(studies)), function(s) {
    ssEnv <- core_init_env(result_folder = studies$STUDY_FOLDER[s],
                           start_fresh = FALSE, ...)
    keys <- ssEnv$keys_areas_subareas_markers_figures
    if (is.null(keys))
      stop("study '", studies$STUDY[s], "' declares no region classes: ",
           "its session carries no key table.")
    unique(keys[, .META_KEY_COLS, drop = FALSE])
  })
  names(key_sets) <- studies$STUDY
  key_sets
}

#' Intersect the studies' region classes
#'
#' Strict: a class is in the space only when every study has it. This is
#' narrower than the pooling itself requires, which needs two studies rather
#' than all of them, and it is chosen deliberately so that every row of the
#' result rests on the same set of studies.
#'
#' What the strictness costs is reported rather than dropped. A class measured
#' by four studies out of five is excluded, and silently excluding it is the
#' difference between a decision and a disappearance.
#'
#' @param key_sets A named list of region-class tables, from
#'   [meta_key_sets_read()].
#'
#' @return A list with \code{space}, the intersection, and \code{excluded}, the
#'   classes at least one study had and the space does not, each with the
#'   studies that lacked it.
#'
#' @keywords internal
#' @noRd
meta_key_space <- function(key_sets)
{
  if (length(key_sets) < 2)
    stop("a meta analysis needs at least two studies; got ", length(key_sets), ".")
  if (is.null(names(key_sets)) || any(!nzchar(names(key_sets))))
    stop("key_sets must be named by study.")

  as_id <- function(df) do.call(paste, c(unname(as.list(df[, .META_KEY_COLS, drop = FALSE])), sep = "\r"))

  ids <- lapply(key_sets, as_id)
  common <- Reduce(intersect, ids)

  first <- key_sets[[1]]
  space <- first[as_id(first) %in% common, .META_KEY_COLS, drop = FALSE]
  space <- space[order(space$MARKER, space$FIGURE, space$AREA, space$SUBAREA), , drop = FALSE]
  rownames(space) <- NULL

  # Every class any study declared, minus the space, each with who lacked it.
  union_ids <- unique(unlist(ids, use.names = FALSE))
  outside <- setdiff(union_ids, common)
  if (length(outside) == 0) {
    excluded <- data.frame(MARKER = character(), FIGURE = character(),
                           AREA = character(), SUBAREA = character(),
                           STUDIES_MISSING = character(), stringsAsFactors = FALSE)
  } else {
    parts <- do.call(rbind, strsplit(outside, "\r", fixed = TRUE))
    excluded <- as.data.frame(parts, stringsAsFactors = FALSE)
    colnames(excluded) <- .META_KEY_COLS
    excluded$STUDIES_MISSING <- vapply(outside, function(id)
      paste(names(ids)[!vapply(ids, function(x) id %in% x, logical(1))], collapse = " "),
      character(1), USE.NAMES = FALSE)
    excluded <- excluded[order(excluded$MARKER, excluded$FIGURE,
                               excluded$AREA, excluded$SUBAREA), , drop = FALSE]
    rownames(excluded) <- NULL
  }

  if (nrow(space) == 0)
    stop("the studies share no region class: there is nothing to meta analyse. ",
         nrow(excluded), " classes were declared by some study but not by all.")

  list(space = space, excluded = excluded)
}
