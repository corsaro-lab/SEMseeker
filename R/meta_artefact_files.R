# Which file of which study a meta analysis is going to read, and whether it is
# there. Both kinds of artefact are named from the inference specification, so
# the question is answerable before any reading starts.
#
# The two kinds are not shaped alike, and the difference is not cosmetic:
#
#   inference   one file per MARKER. The region-class coordinates live inside
#               the file, in its AREA and SUBAREA columns.
#   enrichment  one file per KEY. The coordinates are in the file NAME, because
#               each class was enriched on its own.
#
# So the inference side asks for as many files as there are markers, and the
# enrichment side for as many as there are classes.

#' Name and check the per-study artefacts a meta analysis will read
#'
#' A study with no file for the request contributes nothing, and a missing
#' contribution does not announce itself: it lowers the number of studies behind
#' a result while the output still looks complete. Hence the refusal here rather
#' than a shorter answer later.
#'
#' @param studies A table from [meta_studies_normalise()].
#' @param keys The region classes to read, from [meta_keys_select()].
#' @param inference_detail One row of the inference specification. It is both
#'   the request and the name of the files that answer it.
#' @param artefact Which kind to look for, \code{"inference"} or
#'   \code{"enrichment"}.
#' @param enrichment_package Required for \code{artefact = "enrichment"}: the
#'   enricher whose results are compared, one of \code{"WebGestalt"},
#'   \code{"STRINGdb"}, \code{"ctdR"} or \code{"pathfindR"}.
#' @param pvalue_column,alpha,significance Passed to the enrichment name
#'   builder, which puts them in the file name.
#'
#' @return A \code{data.frame} of \code{STUDY}, \code{MARKER}, the class
#'   coordinates where the artefact has them in its name, and \code{PATH}. Stops
#'   instead when any of them is absent.
#'
#' @keywords internal
#' @noRd
meta_artefact_files <- function(studies, keys, inference_detail,
                                artefact = c("inference", "enrichment"),
                                enrichment_package = NULL,
                                pvalue_column = "PVALUE_ADJ_ALL_BH",
                                alpha = NULL, significance = TRUE)
{
  artefact <- match.arg(artefact)
  if (!is.data.frame(keys) || nrow(keys) == 0)
    stop("keys is empty: there is no artefact to look for.")

  if (artefact == "enrichment") {
    # The enrichers are not a list to keep here: the session declares them, and
    # a hand-written list was already missing two of them.
    ssEnv <- core_get_session_info()
    known <- as.character(ssEnv$key_enrichment_format$label)
    if (is.null(enrichment_package) || length(enrichment_package) != 1)
      stop("artefact = \"enrichment\" needs a single enrichment_package, one of ",
           paste(known, collapse = ", "), ".")
    if (!enrichment_package %in% known)
      stop("unknown enrichment_package '", enrichment_package, "'; the session declares: ",
           paste(known, collapse = ", "), ".")
  }

  rows <- list()
  for (s in seq_len(nrow(studies))) {
    folder <- studies$STUDY_FOLDER[s]
    if (artefact == "inference") {
      # Per marker: the class coordinates are columns inside the file.
      for (marker in unique(as.character(keys$MARKER))) {
        path <- io_inference_file_name(inference_detail, marker,
                                       file.path(folder, "Inference"),
                                       skip_dir_create = TRUE)
        rows[[length(rows) + 1]] <- data.frame(
          STUDY = studies$STUDY[s], MARKER = marker,
          FIGURE = NA_character_, AREA = NA_character_, SUBAREA = NA_character_,
          PATH = path, stringsAsFactors = FALSE)
      }
    } else {
      # Per class: the coordinates are in the file name, so one file each. Both
      # the folder and the name come from the helpers the writers use, hung off
      # this study's enrichment folder rather than the current session's.
      base <- io_enrichment_folder(inference_detail, enrichment_package,
                                   base = file.path(folder, "Enrichment"),
                                   skip_dir_create = TRUE)
      for (k in seq_len(nrow(keys))) {
        project <- enrich_phenotype_analysis_name(
          inference_detail = inference_detail, key = keys[k, , drop = FALSE],
          prefix = "", suffix = "", pvalue_column = pvalue_column,
          alpha = alpha, significance = significance)
        rows[[length(rows) + 1]] <- data.frame(
          STUDY = studies$STUDY[s], MARKER = keys[k, "MARKER"],
          FIGURE = keys[k, "FIGURE"], AREA = keys[k, "AREA"],
          SUBAREA = keys[k, "SUBAREA"],
          PATH = io_file_path_build(base, project, "csv"),
          stringsAsFactors = FALSE)
      }
    }
  }

  out <- do.call(rbind, rows)
  absent <- !file.exists(out$PATH)
  if (any(absent)) {
    described <- if (artefact == "inference")
      paste0("  ", out$STUDY[absent], " / ", out$MARKER[absent])
    else
      paste0("  ", out$STUDY[absent], " / ", out$MARKER[absent], " ",
             out$FIGURE[absent], " ", out$AREA[absent], " ", out$SUBAREA[absent])
    stop("these studies hold no ", artefact, " results for the request being pooled:\n",
         paste(paste0(described, ": ", out$PATH[absent]), collapse = "\n"),
         "\nEither run that step in those studies first, or narrow the request ",
         "to what they have.")
  }
  rownames(out) <- NULL
  out
}
