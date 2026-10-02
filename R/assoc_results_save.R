# Symmetric counterpart to assoc_results_get(): writes an
# association inference CSV to the canonical path computed by
# io_inference_file_name().
#
# Kept deliberately thin (no p-adjust, no SQL filtering, no annotation
# enrichment) - those belong in upstream analyzers, not in the save
# step. Call sites that currently do
#   utils::write.csv2(results_inference,
#                     io_inference_file_name(inference_detail, marker, folder),
#                     row.names = FALSE)
# can switch to a one-line
#   assoc_results_save(results_inference, inference_detail, marker, folder)
# without behaviour change.
#
# @param results_inference data.frame to write
# @param inference_detail   inference metadata row (used to derive file name)
# @param marker             marker identifier (used in the file name)
# @param folder             output folder; if NULL falls back to
#                           ssEnv$result_folderInference
# @keywords internal
assoc_results_save <- function(results_inference, inference_detail, marker,
                                     folder = NULL) {
  if (is.null(folder)) {
    ssEnv <- core_get_session_info()
    folder <- ssEnv$result_folderInference
  }
  file_name <- io_inference_file_name(inference_detail, marker, folder)
  utils::write.csv2(results_inference, file_name, row.names = FALSE)
  invisible(file_name)
}
