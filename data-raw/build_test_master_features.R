#!/usr/bin/env Rscript
## ============================================================================
## data-raw/build_test_master_features.R
##
## One-shot generator for data/test_master_features.rda — a 20k-probe master
## fixture sampled deterministically around 81 known human imprinting DMRs
## (KCNQ1OT1, H19/IGF2, MEG3, GNAS, PEG3, SNURF, PLAGL1, ...).
##
## The same fixture is consumed by:
##   - tests/testthat/setup.R  (automated tests)
##   - vignettes/getting-started.Rmd  (runnable example, subset to BWS loci)
##
## Re-run only when the curated DMR list (data/dmr_annotation.rda) or the EPIC
## annotation package changes.  Requires
## `IlluminaHumanMethylationEPICanno.ilm10b4.hg19` installed (needs XQuartz on
## macOS).
##
## Run from the package root:
##     Rscript data-raw/build_test_master_features.R
##
## set.seed(20210713) — date of v.0.1.9 Zenodo release (DOI 10.5281/zenodo.5095416).
## ============================================================================

suppressPackageStartupMessages({
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(usethis)
})

stopifnot(file.exists("data/dmr_annotation.rda"))
load("data/dmr_annotation.rda")

## ── 1. Load EPIC manifest ──────────────────────────────────────────────────
data("Locations", package = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
     envir = environment())
locs <- as.data.frame(Locations, stringsAsFactors = FALSE)
locs$PROBE <- rownames(locs)
locs$CHR_NO_PREFIX <- sub("^chr", "", locs$chr)
cat(sprintf("[build] EPIC manifest probes: %d\n", nrow(locs)))

## ── 2. Seed: imprinting-DMR probes from curated annotation ─────────────────
seed_df <- unique(dmr_annotation[, c("PROBE", "DMR_WHOLE")])
seed_df <- seed_df[seed_df$PROBE %in% locs$PROBE, ]
seed_set <- unique(seed_df$PROBE)
cat(sprintf("[build] Seed probes (unique DMR-internal): %d in %d DMRs\n",
            length(seed_set), length(unique(seed_df$DMR_WHOLE))))

## ── 3. Imprinting windows per DMR (±50 kb flank) ───────────────────────────
seed_locs <- merge(seed_df, locs[, c("PROBE", "chr", "pos", "CHR_NO_PREFIX")],
                   by = "PROBE", all.x = TRUE)

## ±500 kb captures most imprinted gene CLUSTERS (e.g. 11p15 ICR1+ICR2 span)
## without ballooning the fixture beyond the 1 MB budget.
FLANK <- 500000L
windows <- by(seed_locs, seed_locs$DMR_WHOLE, function(g) {
  data.frame(
    DMR_WHOLE = g$DMR_WHOLE[1L],
    chr       = g$chr[1L],
    win_start = as.integer(min(g$pos) - FLANK),
    win_end   = as.integer(max(g$pos) + FLANK),
    stringsAsFactors = FALSE
  )
})
windows <- do.call(rbind, windows)
rownames(windows) <- NULL
cat(sprintf("[build] Imprinting windows: %d (flank ±%d bp)\n",
            nrow(windows), FLANK))

## ── 4. Pull every EPIC probe inside any imprinting window ──────────────────
window_probes_list <- vector("list", nrow(windows))
for (i in seq_len(nrow(windows))) {
  w <- windows[i, ]
  hits <- locs[locs$chr == w$chr &
                 locs$pos >= w$win_start &
                 locs$pos <= w$win_end, ]
  if (nrow(hits) > 0L) {
    hits$DMR_WHOLE <- w$DMR_WHOLE
    window_probes_list[[i]] <- hits
  }
}
window_probes <- do.call(rbind, window_probes_list)
## A probe in overlapping windows would appear twice; dedup keeping first DMR
window_probes <- window_probes[!duplicated(window_probes$PROBE), ]
cat(sprintf("[build] Probes within imprinting windows (pre-cap): %d\n",
            nrow(window_probes)))

## ── 5. Cap at 20 000 — keep ALL seed (signal), sample flanking ─────────────
TARGET <- 20000L
set.seed(20210713)

is_seed       <- window_probes$PROBE %in% seed_set
signal_probes <- window_probes[is_seed, ]
flanking      <- window_probes[!is_seed, ]
cat(sprintf("[build] Signal probes kept: %d ; flanking pool: %d\n",
            nrow(signal_probes), nrow(flanking)))

if (nrow(signal_probes) >= TARGET) {
  out <- signal_probes[sample(nrow(signal_probes), TARGET), ]
} else {
  need <- TARGET - nrow(signal_probes)
  if (nrow(flanking) > need) {
    flanking <- flanking[sample(nrow(flanking), need), ]
  }
  out <- rbind(signal_probes, flanking)

  ## ── Top-up: if window-based pool is still below TARGET, fill with random
  ## EPIC probes drawn from outside any imprinting window. Keeps statistical
  ## validity (20k for IQR robustness) AND preserves biology-aware core.
  if (nrow(out) < TARGET) {
    used <- out$PROBE
    pool_off <- locs[!locs$PROBE %in% used, ]
    need_random <- TARGET - nrow(out)
    extra_idx <- sample(nrow(pool_off), need_random)
    extra <- pool_off[extra_idx, ]
    extra$DMR_WHOLE <- NA_character_
    extra <- extra[, c("PROBE", "chr", "pos", "DMR_WHOLE", "CHR_NO_PREFIX")]
    out <- rbind(out, extra)
    cat(sprintf("[build] Top-up random probes (background): %d\n", need_random))
  }
}

## ── 6. Final schema ────────────────────────────────────────────────────────
test_master_features <- data.frame(
  PROBE     = out$PROBE,
  CHR       = out$CHR_NO_PREFIX,
  START     = as.integer(out$pos),
  END       = as.integer(out$pos),
  ABSOLUTE  = paste(out$CHR_NO_PREFIX, out$pos, sep = "_"),
  DMR_LABEL = ifelse(out$PROBE %in% seed_set, out$DMR_WHOLE, NA_character_),
  stringsAsFactors = FALSE
)

## Sort by chr/start so downstream consumers get a stable ordering
## Chromosomes ordered as 1..22, X, Y for human readability
.chr_levels <- c(as.character(1:22), "X", "Y", "M", "MT")
test_master_features$CHR <- factor(test_master_features$CHR,
                                   levels = .chr_levels)
test_master_features <- test_master_features[
  order(test_master_features$CHR, test_master_features$START), ]
test_master_features$CHR <- as.character(test_master_features$CHR)
rownames(test_master_features) <- NULL

cat(sprintf("[build] Final dim: %d rows × %d cols\n",
            nrow(test_master_features), ncol(test_master_features)))
cat(sprintf("[build] Signal probes (with DMR_LABEL): %d\n",
            sum(!is.na(test_master_features$DMR_LABEL))))
cat("[build] Top 10 DMRs by probe count:\n")
.top <- sort(table(test_master_features$DMR_LABEL), decreasing = TRUE)
print(head(.top, 10L))

## ── 7. Save as package data ────────────────────────────────────────────────
usethis::use_data(test_master_features, overwrite = TRUE, compress = "xz")
.size_kb <- file.info("data/test_master_features.rda")$size / 1024
cat(sprintf("[build] data/test_master_features.rda saved (%.1f KB)\n", .size_kb))
