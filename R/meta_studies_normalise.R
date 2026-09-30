#' Normalise the studies argument of the meta_ family
#'
#' Every function in the \code{meta_} family reads the results of several
#' studies, and each one used to ask for them differently. This turns the three
#' accepted forms into one table.
#'
#' A study needs two things, a label and a folder, and the label is what keys
#' the output: two studies that end up with the same label merge into one row
#' and nobody is told. That is why duplicate labels are refused here rather
#' than discovered later.
#'
#' The column is called \code{STUDY_FOLDER} and not \code{RESULT_FOLDER}
#' because the meta_ functions also take a \code{result_folder}, which is where
#' they write. The two point in opposite directions and were one letter apart.
#'
#' @param studies One of: a character vector of study paths, where the final
#'   path component labels the study; a named character vector, where the name
#'   labels it; or a \code{data.frame} with \code{STUDY} and
#'   \code{RESULT_FOLDER} columns, the form the overlap functions already take.
#'
#' @return A \code{data.frame} with columns \code{STUDY} and
#'   \code{STUDY_FOLDER}, in the order given.
#'
#' @keywords internal
#' @noRd
meta_studies_normalise <- function(studies)
{
  if (is.data.frame(studies)) {
    missing_cols <- setdiff(c("STUDY", "RESULT_FOLDER"), colnames(studies))
    if (length(missing_cols) > 0)
      stop("studies given as a data.frame must have columns STUDY and RESULT_FOLDER; missing: ",
           paste(missing_cols, collapse = ", "))
    out <- data.frame(STUDY = as.character(studies$STUDY),
                      STUDY_FOLDER = as.character(studies$RESULT_FOLDER),
                      stringsAsFactors = FALSE)
  } else {
    if (!is.character(studies))
      stop("studies must be a character vector of study paths, or a data.frame ",
           "with STUDY and RESULT_FOLDER columns; got ", class(studies)[1])
    folders <- as.character(studies)
    labels <- names(studies)
    if (is.null(labels)) labels <- rep(NA_character_, length(folders))
    # An unnamed element is labelled by the final component of its path. Trailing
    # separators are dropped first, otherwise basename() of "/a/b/" and of "/a/b"
    # would not agree.
    bare <- !nzchar(labels) | is.na(labels)
    labels[bare] <- basename(sub("[/\\\\]+$", "", folders[bare]))
    out <- data.frame(STUDY = labels, STUDY_FOLDER = folders, stringsAsFactors = FALSE)
  }

  if (nrow(out) == 0)
    stop("studies is empty: there is nothing to read.")
  if (any(!nzchar(out$STUDY_FOLDER)) || any(is.na(out$STUDY_FOLDER)))
    stop("studies contains an empty path.")
  if (any(!nzchar(out$STUDY)) || any(is.na(out$STUDY)))
    stop("studies contains a path with no usable label: name that element explicitly.")

  duplicated_labels <- unique(out$STUDY[duplicated(out$STUDY)])
  if (length(duplicated_labels) > 0)
    stop("two or more studies carry the same label (", paste(duplicated_labels, collapse = ", "),
         "). The label keys the result, so they would merge into one row. ",
         "Name the elements explicitly to tell them apart.")

  # A folder that exists is not a study. Two things make one, and both are
  # checked here rather than discovered halfway through a long read.
  #
  # The session written to disk, because that is where the region classes the
  # study analysed are recorded: meta_key_sets_read() takes them from there to
  # build the space of the meta analysis, so a folder without it is not a study
  # whose classes are unknown, it is a folder we cannot ask the question of.
  #
  # And the inference results, in the place the readers look for them. Whether a
  # study holds the results of a PARTICULAR request is a sharper question, and
  # meta_inference_files_check() answers it once the request is known; this only
  # establishes that it is a study at all.
  reasons <- vapply(out$STUDY_FOLDER, function(folder) {
    if (!dir.exists(folder)) return("the folder does not exist")
    if (!file.exists(file.path(folder, "Log", "session_info.rds")))
      return("no session was written here (Log/session_info.rds is missing)")
    inference <- file.path(folder, "Inference")
    if (!dir.exists(inference))
      return("there is no Inference folder")
    if (length(list.files(inference, recursive = TRUE)) == 0)
      return("the Inference folder is empty")
    ""
  }, character(1), USE.NAMES = FALSE)

  bad <- nzchar(reasons)
  if (any(bad))
    stop("these paths are not SEMseeker studies:\n",
         paste0("  ", out$STUDY_FOLDER[bad], ": ", reasons[bad], collapse = "\n"))

  out
}
