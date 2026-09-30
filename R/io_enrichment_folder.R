#' The folder an enrichment result of a given enricher lives in
#'
#' The same four-part convention was written out at eight call sites, differing
#' only in the label. It is one place now, which is what lets a reader in
#' another study's folder be named without rebuilding the convention: the defect
#' that pattern produces is a path that agrees with nothing and finds no files,
#' which looks like an empty result rather than a mistake.
#'
#' @param inference_detail One row of the inference specification. Its three sql
#'   conditions are the last three components of the path.
#' @param label The enricher's label, as it appears in
#'   \code{ssEnv$key_enrichment_format$label}.
#' @param base The folder the convention hangs off. Defaults to the current
#'   session's enrichment folder; pass another study's to name a file there, or
#'   the phenotype folder for the phenotype side.
#' @param skip_dir_create \code{TRUE} returns the path without creating
#'   anything, for callers that want the name and must not touch the filesystem.
#'
#' @return The folder path.
#'
#' @keywords internal
#' @noRd
io_enrichment_folder <- function(inference_detail, label, base = NULL,
                                 skip_dir_create = FALSE)
{
  if (is.null(base)) {
    ssEnv <- core_get_session_info()
    base <- ssEnv$result_folderEnrichment
  }
  parts <- c(as.character(label),
             core_name_cleaning(inference_detail$areas_sql_condition),
             core_name_cleaning(inference_detail$samples_sql_condition),
             core_name_cleaning(inference_detail$association_results_sql_condition))
  if (skip_dir_create) do.call(file.path, as.list(c(base, parts)))
  else io_dir_check_and_create(base, parts)
}
