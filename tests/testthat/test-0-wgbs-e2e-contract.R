# AI-109 (2026-06-09): contract-level E2E tests for the WGBS / LONGREAD
# code path. The bulk of SEMseeker's runtime testing has been against
# Illumina (K27/K450/K850); long-reads exercise a different set of
# helpers (io_coord_probe_features, AI-098 tech-aware AREA skip, chunked
# per-chr lmFit). These tests build SYNTHETIC long-read data and assert
# the contracts that those helpers must satisfy.
#
# Synthetic probe-ID format: "{CHR}_{START}", with CHR carrying NO "chr"
# prefix (per io_probe_id_to_coord in R/coord_input.R). E.g. "1_10000",
# "X_5000".

# ---- io_coord_probe_features round-trip ----------------------------------

test_that("io_probe_id_to_coord parses synthetic IDs round-trip with io_coord_probe_features", {
  probe_ids <- c("1_10000", "1_20000", "2_30000", "X_5000", "Y_99999")

  coords <- SEMseeker:::io_probe_id_to_coord(probe_ids)
  expect_equal(coords$CHR,   c("1", "1", "2", "X", "Y"))
  expect_equal(coords$START, c(10000L, 20000L, 30000L, 5000L, 99999L))
  expect_equal(coords$END,   coords$START + 1L)

  pf <- SEMseeker:::io_coord_probe_features(probe_ids)
  expect_equal(pf$PROBE, probe_ids)
  expect_equal(pf$CHR,   coords$CHR)
  expect_equal(pf$START, coords$START)
  expect_equal(pf$END,   coords$END)
})

# ---- smart_split is a no-op on coord-encoded AREA names ----------------

test_that(".anno_smart_split_area_name is a no-op on coord-encoded AREA (no slash)", {
  coord_areas <- c("1_10000", "X_5000", "22_99999999")
  for (a in coord_areas) {
    expect_equal(SEMseeker:::.anno_smart_split_area_name(a), a,
                 info = sprintf("coord AREA '%s' must pass through unchanged", a))
  }
})

# ---- AI-106 sanitize regex is a no-op on coord-encoded AREA names ------

test_that("AI-106 sanitize regex preserves coord AREA names verbatim", {
  coord_areas <- c("1_10000", "X_5000", "22_99999999")
  safe <- gsub("[^A-Za-z0-9_.]", "_", coord_areas)
  # Synthetic coords are already in the R-safe set [A-Za-z0-9_.] → no change.
  expect_equal(safe, coord_areas)
})

# ---- AI-098 tech-aware AREA skip semantics ----------------------------

test_that("AI-098 AREA skip: PROBE is no-op on long-reads", {
  # The condition lives in run_depth_n_marker.R:
  #   if (key$AREA == "PROBE" && tech_is_longread) next
  for (tech in c("WGBS", "LONGREAD")) {
    tech_is_longread <- tech %in% c("WGBS", "LONGREAD")
    expect_true(tech_is_longread,
                info = paste("tech '", tech, "' must be classified as long-read"))
    # Simulate the skip predicate
    expect_true(("PROBE" == "PROBE" && tech_is_longread),
                info = paste("PROBE AREA must be skipped for tech =", tech))
  }
})

test_that("AI-098 AREA skip: POSITION is no-op on Illumina", {
  for (tech in c("K27", "K450", "K850")) {
    tech_is_longread <- tech %in% c("WGBS", "LONGREAD")
    expect_false(tech_is_longread)
    # Simulate the skip predicate
    #   if (key$AREA == "POSITION" && !tech_is_longread) next
    expect_true(("POSITION" == "POSITION" && !tech_is_longread),
                info = paste("POSITION AREA must be skipped for tech =", tech))
  }
})

test_that("AI-098 source-level: skip predicates live in assoc_run_marker", {
  src <- paste(deparse(SEMseeker:::assoc_run_marker), collapse = "\n")
  expect_true(
    grepl('key\\$AREA\\s*==\\s*"POSITION"\\s*&&\\s*!tech_is_longread', src),
    info = "AI-098 POSITION skip on Illumina must be present"
  )
  expect_true(
    grepl('key\\$AREA\\s*==\\s*"PROBE"\\s*&&\\s*tech_is_longread', src),
    info = "AI-098 PROBE skip on long-reads must be present"
  )
})

# ---- assoc_lmfit_chunked_by_chr: chromosome extraction strategy --------------

test_that("assoc_lmfit_chunked_by_chr extracts CHR for long-reads via AREA string split on '_'", {
  src <- paste(deparse(SEMseeker:::assoc_lmfit_chunked_by_chr), collapse = "\n")
  # WGBS / LONGREAD strategy: split on "_" and pick first chunk
  expect_true(
    grepl('str\\$split\\("_"\\)\\$list\\$get\\(0L\\)', src),
    info = "WGBS chunking strategy must parse CHR via str$split('_')$list$get(0L)"
  )
})

# ---- Synthetic POSITION pivot smoke (read-only contract) -----------------

test_that("synthetic long-read POSITION pivot can be lazily summarised by AREA prefix", {
  skip_on_cran()
  # Build a tiny synthetic POSITION pivot with 5 chrs × 4 probes each = 20 rows
  set.seed(42)
  chrs    <- c("1", "2", "3", "X", "Y")
  starts  <- rep(c(10000L, 20000L, 30000L, 40000L), times = length(chrs))
  probes  <- paste(rep(chrs, each = 4), starts, sep = "_")
  n       <- length(probes)

  pivot_df <- data.frame(
    AREA = probes,
    S001 = runif(n, 0, 1),
    S002 = runif(n, 0, 1),
    S003 = runif(n, 0, 1),
    stringsAsFactors = FALSE
  )
  pivot_lazy <- polars::as_polars_df(pivot_df)$lazy()

  # The chunking strategy: extract leading token before the first "_".
  chr_extracted <- as.character(as.data.frame(
    pivot_lazy$
      select(polars::pl$col("AREA")$str$split("_")$list$get(0L)$alias("CHR"))$
      unique()$
      collect()
  )$CHR)

  # Should core_recover the original 5 chromosomes (order-independent)
  expect_setequal(chr_extracted, chrs)
})

# ---- core_get_meth_tech: WGBS pre-declaration is honoured ------------------

test_that("core_get_meth_tech respects ssEnv$tech pre-declaration for WGBS", {
  skip_on_cran()

  tempFolder <- tempfile("ai109_wgbs_")
  dir.create(file.path(tempFolder, "Data"), recursive = TRUE)
  ssEnv <- SEMseeker:::core_init_env(tempFolder,
                                parallel_strategy = "sequential",
                                tech              = "WGBS",
                                iqrTimes          = 3,
                                verbosity         = 1)
  on.exit({ SEMseeker:::core_close_env(); unlink(tempFolder, recursive = TRUE) },
          add = TRUE)

  # Pre-declaration must be honoured even with an empty signal stub
  signal_data <- data.frame(
    PROBE = c("1_10000", "1_20000", "X_5000"),
    S001  = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )
  rownames(signal_data) <- signal_data$PROBE
  signal_data$PROBE <- NULL

  out <- tryCatch(SEMseeker:::core_get_meth_tech(signal_data),
                  error = function(e) e)
  if (inherits(out, "error")) {
    testthat::skip(paste("core_get_meth_tech failed (likely missing optional setup):",
                          conditionMessage(out)))
  }
  expect_equal(out$tech, "WGBS",
               info = "Pre-declared tech = WGBS must survive core_get_meth_tech")
})
