# =============================================================================
# Rewrite the golden-master reference (tests/testthat/golden/).
#
# Run it only when a change of results is INTENDED, from the package root, on
# the INSTALLED package (install first, or it fingerprints the previous tree):
#
#   R CMD INSTALL . && NOT_CRAN=true Rscript dev/golden-update.R
#
# Then commit tests/testthat/golden/ together with the change: the diff of
# index.csv, and the failure the test printed before the update, are the record
# of which results moved and why.
# =============================================================================

suppressPackageStartupMessages(library(SEMseeker))
source(file.path("tests", "testthat", "helper-golden.R"))

result_folder <- file.path(tempdir(), "golden-update")
unlink(result_folder, recursive = TRUE)
.golden_run(result_folder)
fingerprint <- .golden_fingerprint(result_folder)

dir <- file.path("tests", "testthat", "golden")
if (file.exists(file.path(dir, "fingerprint.rds"))) {
  moved <- .golden_compare(readRDS(file.path(dir, "fingerprint.rds")), fingerprint)
  cat(length(moved), "difference(s) from the previous reference\n")
  if (length(moved)) cat(paste0("  ", utils::head(moved, 50), collapse = "\n"), "\n")
}
index <- .golden_write(fingerprint, dir)
cat("reference written:", nrow(index), "files in", dir, "\n")
