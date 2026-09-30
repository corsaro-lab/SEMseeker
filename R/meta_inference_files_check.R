#' Check that every study holds the results of the request being pooled
#'
#' The inference specification is both the meta analysis's input and the name of
#' the files it will read: \code{io_inference_file_name()} derives one from the
#' other. So whether a study can answer this particular request is knowable
#' before any reading starts, and it is checked here rather than discovered
#' after a long pass over the ones that can.
#'
#' It matters more than it looks. A study with no file for the request
#' contributes nothing, and a missing contribution does not announce itself: it
#' lowers the number of studies behind a pooled estimate, or with a strict
#' intersection empties the result entirely, and either way the output still
#' looks complete.
#'
#' @param studies A table from [meta_studies_normalise()].
#' @param keys The region classes to analyse, from [meta_keys_select()]. Only
#'   their markers are used: the file name is per marker, the class coordinates
#'   live inside the file.
#' @param inference_detail One row of the inference specification.
#'
#' @return The studies table, unchanged, invisibly. Called for the refusal.
#'
#' @keywords internal
#' @noRd
meta_inference_files_check <- function(studies, keys, inference_detail)
{
  markers <- unique(as.character(keys$MARKER))
  if (length(markers) == 0)
    stop("keys names no marker: there is no file to look for.")

  missing <- character()
  for (s in seq_len(nrow(studies))) {
    inference_folder <- file.path(studies$STUDY_FOLDER[s], "Inference")
    for (marker in markers) {
      wanted <- io_inference_file_name(inference_detail, marker, inference_folder,
                                       skip_dir_create = TRUE)
      if (!file.exists(wanted))
        missing <- c(missing, paste0("  ", studies$STUDY[s], " / ", marker, ": ", wanted))
    }
  }

  if (length(missing) > 0)
    stop("these studies hold no results for the request being pooled:\n",
         paste(missing, collapse = "\n"),
         "\nEither run the inference in those studies first, or narrow the ",
         "request to the markers they have.")

  invisible(studies)
}
