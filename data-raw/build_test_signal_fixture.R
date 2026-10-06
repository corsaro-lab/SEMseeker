#!/usr/bin/env Rscript
## ============================================================================
## data-raw/build_test_signal_fixture.R
##
## Build two package data fixtures from GSE133774 (EPIC 850k, BWS + MLID,
## GPL21145 — Infinium MethylationEPIC BeadChip):
##
##   data/test_signal_gse133774.rda      — beta matrix (≤20k probes × 10 samples)
##   data/test_samplesheet_gse133774.rda — sample sheet (16 rows: 3-class design)
##
## Source data:
##   semseeker_validationE2E/data/GSE13374/GSE133774_series_matrix.txt.gz (44 MB)
##
## Why GSE133774:
##   - EPIC 850k → full overlap with test_master_features (20k EPIC probe IDs)
##   - Contains real BWS + MLID epimutations (L1 = Beckwith-Wiedemann syndrome)
##   - Same dataset selected for the getting-started vignette
##   - Single consistent source for automated tests AND vignette
##
## Sample design (3-class SEMseeker pattern, Reference reuse):
##   Reference: CTRL01–CTRL06 (6 controls)
##   Control:   CTRL01–CTRL06 (same 6 controls, reused — SEMseeker canonical)
##   Case:      L1, L2, L3, L4 (BWS family — L1 is the BWS proband)
##
## Run from package root:
##     Rscript data-raw/build_test_signal_fixture.R
##
## Requirements: GEOquery, Biobase, usethis installed.
## ============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(usethis)
})

VAL_DIR    <- "~/Documents/Progetti/RICERCA/semseeker_validationE2E/data/GSE13374"
SERIES_MAT <- file.path(VAL_DIR, "GSE133774_series_matrix.txt.gz")
DERIVED    <- file.path(VAL_DIR, "derived")

stopifnot(file.exists(SERIES_MAT))
stopifnot(file.exists("data/test_master_features.rda"))

## ── 1. Parse beta matrix from series matrix (EPIC 850k) ───────────────────
beta_cache <- file.path(DERIVED, "beta.rds")
if (file.exists(beta_cache) &&
    file.info(beta_cache)$mtime > file.info(SERIES_MAT)$mtime) {
  cat("[fixture] cache OK:", beta_cache, "\n")
  beta_full <- readRDS(beta_cache)
} else {
  cat("[fixture] parsing GSE133774 series matrix (44 MB gz)...\n")
  dir.create(DERIVED, recursive = TRUE, showWarnings = FALSE)
  gset <- GEOquery::getGEO(filename = SERIES_MAT, GSEMatrix = TRUE, getGPL = FALSE)
  if (is.list(gset)) gset <- gset[[1]]
  beta_full <- Biobase::exprs(gset)
  if (max(beta_full, na.rm = TRUE) > 1.5 || min(beta_full, na.rm = TRUE) < -0.5) {
    cat("[fixture] M-values detected; converting to beta\n")
    beta_full <- (2 ^ beta_full) / (2 ^ beta_full + 1)
  }
  saveRDS(beta_full, beta_cache)
  cat(sprintf("[fixture] beta_full: %d probes × %d samples\n",
              nrow(beta_full), ncol(beta_full)))
}

## ── 2. Filter to 20k imprinting-aware EPIC probe IDs ──────────────────────
load("data/test_master_features.rda")
target_probes <- test_master_features$PROBE
overlap <- intersect(target_probes, rownames(beta_full))
cat(sprintf("[fixture] overlap: %d / %d target probes\n",
            length(overlap), length(target_probes)))
if (length(overlap) < 5000L)
  stop("[fixture] overlap too small — check probe IDs")

test_signal_gse133774 <- beta_full[overlap, , drop = FALSE]

## Rename columns from GEO accession IDs to readable labels matching the
## series matrix sample titles (CTRL01–CTRL06, L1–L4)
sample_titles <- c("CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06",
                   "L1","L2","L3","L4")
if (ncol(test_signal_gse133774) == length(sample_titles)) {
  colnames(test_signal_gse133774) <- sample_titles
} else {
  cat(sprintf("[fixture] WARNING: expected 10 samples, got %d — using GSM IDs\n",
              ncol(test_signal_gse133774)))
  sample_titles <- colnames(test_signal_gse133774)
}

cat(sprintf("[fixture] test_signal_gse133774: %d probes × %d samples\n",
            nrow(test_signal_gse133774), ncol(test_signal_gse133774)))

## ── 3. Build 3-class sample sheet (Reference reuse pattern) ───────────────
## Reference = Control (CTRL01-06 appear twice with different roles).
## This is the canonical SEMseeker design documented in the vignette.
test_samplesheet_gse133774 <- data.frame(
  Sample_ID = c(
    "CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06",   # Reference
    "CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06",   # Control (reuse)
    "L1","L2","L3","L4"                                       # Case (BWS family)
  ),
  Sample_Group = c(
    rep("Reference", 6L),
    rep("Control",   6L),
    rep("Case",      4L)
  ),
  stringsAsFactors = FALSE
)
cat(sprintf("[fixture] samplesheet: %d rows (%s)\n",
            nrow(test_samplesheet_gse133774),
            paste(names(table(test_samplesheet_gse133774$Sample_Group)),
                  table(test_samplesheet_gse133774$Sample_Group), sep="=",
                  collapse=", ")))

## ── 4. Save as package data ───────────────────────────────────────────────
usethis::use_data(test_signal_gse133774,      overwrite = TRUE, compress = "xz")
usethis::use_data(test_samplesheet_gse133774, overwrite = TRUE, compress = "xz")

sz1 <- file.info("data/test_signal_gse133774.rda")$size / 1024
sz2 <- file.info("data/test_samplesheet_gse133774.rda")$size / 1024
cat(sprintf("[fixture] test_signal_gse133774.rda      saved (%.1f KB)\n", sz1))
cat(sprintf("[fixture] test_samplesheet_gse133774.rda saved (%.1f KB)\n", sz2))
cat("[fixture] DONE\n")
cat("[fixture] Next: update setup.R + test-cross-format-convergence.R\n")
cat("[fixture]        + vignettes/getting-started.Rmd to load these fixtures\n")
cat("[fixture]        + add roxygen docs in R/data.R for both datasets\n")
