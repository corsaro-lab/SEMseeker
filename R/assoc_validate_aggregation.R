#' Validate the requested aggregations against the markers of the run
#'
#' AI-248. Semantic validation of `inference_details$aggregation`, run **at the
#' door** of `association_analysis()` — before the session does any work, before
#' a single result row is written. Checking it deeper, inside the per-marker
#' loop, would mean discovering the mistake after part of the output already
#' exists.
#'
#' It is the counterpart of `assoc_validate_inference_schema()`, which checks
#' the *shape* of the request (which columns are legal). This one checks the
#' *meaning*: that the aggregation exists in the taxonomy, and that at least one
#' marker of this run admits it. Asking for the median of a 0/1 marker is
#' refused here rather than answered with an empty result.
#'
#' The check lives here and not inside [io_feature_colname()] on purpose: the
#' registry ([util_aggregations_allowed()]) says what is legal, the compositor
#' says how it is named. Making the compositor validate would couple naming to
#' the domain vocabulary and would break the callers that legitimately compose
#' hypothetical names — the scope probe does exactly that to find out which
#' columns a sibling actually carries.
#'
#' Two failure modes, treated differently on purpose:
#' \itemize{
#'   \item a **malformed** request — no aggregation named, or a name that is not
#'     in the taxonomy — stops the run. It is a mistake in how the request was
#'     written, and no rewriting of it can be trusted to mean what the author
#'     intended;
#'   \item an **impossible** request — a legal aggregation that none of this
#'     run's markers admits, e.g. the median of a 0/1 marker — does not stop
#'     anything: the row is dropped with a loud warning and the remaining rows
#'     run. Losing one row of a batch is better than losing the batch, as long
#'     as the loss is said out loud.
#' }
#' A request that is admissible for *some* markers keeps all its rows; the
#' combinations that cannot be computed are named in a warning and skipped
#' downstream.
#'
#' AI-308: it runs **after** [assoc_validate_scope()], and needs to. Which
#' aggregations a marker admits depends on the scope: the two peaks of a
#' bimodal density need one big group, so they exist at `SCOPE = SAMPLE` and
#' nowhere else, so a request whose scope has not been validated yet cannot be
#' judged here. `inference_details$scope` is therefore a precondition, not an
#' optional column.
#'
#' @param inference_details validated data.frame of requests, with `scope`
#'   already validated by [assoc_validate_scope()].
#' @param keys marker/figure keys of the run, defaults to the session's.
#' @return `inference_details`, possibly with impossible rows removed. Callers
#'   must use the returned value.
#' @keywords internal
#' @noRd
assoc_validate_aggregation <- function(inference_details, keys = NULL) {

  if (is.null(keys)) {
    ssEnv <- core_get_session_info()
    keys  <- ssEnv$keys_markers_figures
  }
  if (is.null(keys) || nrow(keys) == 0)
    return(inference_details)

  impossible <- logical(nrow(inference_details))

  # Spell-check only: the semantic check against the keys of the run is below.
  legal <- util_aggregation_vocabulary()

  for (z in seq_len(nrow(inference_details))) {
    detail <- inference_details[z, ]

    # AI-255: every artefact, whatever its scope. A per-area value is an
    # aggregate of the positions of that area exactly as a per-sample value is
    # an aggregate of the positions it masks. If the request does not name the
    # operator, it does not identify what it is asking for.
    requested <- detail$aggregation
    if (is.null(requested) || length(requested) == 0 || all(is.na(requested)) ||
        !any(nzchar(as.character(requested))))
      stop("inference_details row ", z, ": 'aggregation' is required. ",
           "Name which aggregation of the feature to test (",
           paste(legal, collapse = ", "), "). It used to be implicit because ",
           "every marker admitted exactly one; a scope now carries several.",
           call. = FALSE)

    requested <- core_name_cleaning(as.character(requested)[1])
    if (!(requested %in% legal))
      stop("inference_details row ", z, ": aggregation = '", requested,
           "' is not an aggregation of the taxonomy. Legal names: ",
           paste(legal, collapse = ", "), ".", call. = FALSE)

    # AI-308: the scope of the request reaches the registry. Which aggregations
    # a marker admits is not a property of the marker alone: MODELOW/MODEHIGH
    # estimate two peaks of a density and need the whole distribution, so they
    # exist at SCOPE = SAMPLE and nowhere else. Without the scope this check
    # could not see that, and a request for them per instance was answered here
    # and refused later, inside io_pivot_build(), with the run already going.
    #
    # The AREA is deliberately NOT passed. Its restriction, a single-position
    # block admits only VALUE, is not a refusal but a renaming: the request is
    # honoured and the artefact takes the name that says nothing was reduced
    # (see .assoc_aggregation_get()). Refusing it here would reject a legal
    # request on the grounds that one of the run's classes needs it spelled
    # differently.
    scope <- io_scope_validate(detail$scope)
    admitted <- vapply(seq_len(nrow(keys)), function(i)
      requested %in% util_aggregations_allowed(keys$MARKER[i], keys$FIGURE[i],
                                               discrete = isTRUE(keys$DISCRETE[i]),
                                               default  = FALSE,
                                               scope    = scope),
      logical(1))
    if (!any(admitted)) {
      impossible[z] <- TRUE
      message_text <- paste0(
        "inference_details row ", z, ": aggregation '", requested,
        "' is admissible for none of the markers of this run (",
        paste(unique(keys$MARKER), collapse = ", "),
        "), so the row is dropped. A count marker carries SUM and MEAN — its ",
        "median, variance and IQR are degenerate on a vector of zeros and ",
        "ones; the two modes exist only for SIGNAL on the BETA scale.")
      core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " ",
                     message_text)
      warning(message_text, call. = FALSE)
      next
    }

    # Some admit it, some do not: the run goes ahead on the combinations that
    # are possible, and says out loud which ones it is leaving out. Silence here
    # would look exactly like a run that tested everything.
    if (!all(admitted)) {
      excluded <- unique(paste(keys$MARKER[!admitted], keys$FIGURE[!admitted],
                               sep = "/"))
      message_text <- paste0(
        "aggregation '", requested, "' is not admissible for ",
        paste(excluded, collapse = ", "),
        ": those combinations are skipped, the run continues on the others.")
      core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " ",
                     message_text)
      warning(message_text, call. = FALSE)
    }
  }

  if (any(impossible)) {
    inference_details <- inference_details[!impossible, , drop = FALSE]
    if (nrow(inference_details) == 0)
      stop("every row of inference_details asks for an aggregation that no ",
           "marker of this run admits: there is nothing left to analyse.",
           call. = FALSE)
  }

  inference_details
}
