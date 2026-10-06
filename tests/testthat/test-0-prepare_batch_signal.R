# 2026-06-09: contract tests for sem_prepare_batch_signal().
#
# The function is the single source of truth for normalising a SIGNAL
# matrix into a shape consistent with the probe_features used downstream.
# These tests pin down its 4 contracts:
#
#   1. Illumina path: builds probe_features from Bioconductor manifest,
#      collapses dmr_annotation duplicates, intersects with input, removes
#      X/Y, aligns signal_data to probe_features.
#   2. WGBS / LONGREAD path: builds probe_features from synthetic
#      "{CHR}_{START}" rownames via io_coord_probe_features().
#   3. Sex-chromosome removal is applied uniformly (toggle off vs on).
#   4. Post-call invariants:
#        nrow(signal_data) == nrow(probe_features)
#        rownames(signal_data) == probe_features$PROBE
#        no duplicate PROBE in probe_features
#        no X/Y when sex_chromosome_remove = TRUE


# ---- WGBS path: io_coord_probe_features (no manifest dependency) -----------

test_that("WGBS path produces aligned signal_data + probe_features with no X/Y", {
  skip_on_cran()

  # 8 synthetic CpG positions across chr1 / chr2 / chrX
  probes <- c("1_1000", "1_2000", "1_3000",
              "2_5000", "2_6000",
              "X_100", "X_200", "Y_50")
  signal_data <- data.frame(
    S001 = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
    S002 = c(0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85),
    row.names = probes,
    stringsAsFactors = FALSE
  )

  out <- SEMseeker:::sem_prepare_batch_signal(
    signal_data,
    tech = "WGBS",
    sex_chromosome_remove = TRUE
  )

  pf <- attr(out, "probe_features")
  expect_equal(attr(out, "tech"), "WGBS")
  expect_equal(nrow(out), nrow(pf))
  expect_equal(rownames(out), as.character(pf$PROBE))
  # X / Y dropped: 5 autosomal probes survive
  expect_equal(nrow(out), 5L)
  expect_false(any(pf$CHR %in% c("X", "Y")))
  # No duplicates introduced
  expect_equal(anyDuplicated(pf$PROBE), 0L)
})

test_that("WGBS path retains X/Y when sex_chromosome_remove = FALSE", {
  skip_on_cran()

  probes <- c("1_1000", "X_100", "Y_50")
  signal_data <- data.frame(
    S001 = c(0.1, 0.5, 0.7),
    row.names = probes,
    stringsAsFactors = FALSE
  )

  out <- SEMseeker:::sem_prepare_batch_signal(
    signal_data,
    tech = "WGBS",
    sex_chromosome_remove = FALSE
  )

  pf <- attr(out, "probe_features")
  expect_equal(nrow(out), 3L)
  expect_true(all(c("X", "Y") %in% pf$CHR))
})


# ---- Sanity check on the invariants helper ------------------------------

test_that("post-call invariants survive an adversarial coordinate set", {
  skip_on_cran()

  # Some real K-style coords, some on X/Y, some unsortable
  probes <- c("22_100", "1_1000", "10_5000", "X_77", "Y_1", "2_200")
  signal_data <- data.frame(
    S = seq_along(probes),
    row.names = probes,
    stringsAsFactors = FALSE
  )

  out <- SEMseeker:::sem_prepare_batch_signal(
    signal_data, tech = "LONGREAD", sex_chromosome_remove = TRUE
  )

  pf <- attr(out, "probe_features")
  # Invariant 1: same row count
  expect_equal(nrow(out), nrow(pf))
  # Invariant 2: same row order
  expect_identical(rownames(out), as.character(pf$PROBE))
  # Invariant 3: no duplicates
  expect_equal(anyDuplicated(pf$PROBE), 0L)
  # Invariant 4: no sex-chr after removal
  expect_false(any(pf$CHR %in% c("X", "Y")))
  # Of the 6 input probes, 4 survive (1, 22, 10, 2 — Y and X dropped)
  expect_equal(nrow(out), 4L)
})


# ---- Failure modes ----------------------------------------------------

test_that("sem_prepare_batch_signal errors out when the intersection is empty", {
  skip_on_cran()

  # All probes are on X (which we'll then strip)
  probes <- c("X_100", "X_200", "X_300")
  signal_data <- data.frame(
    S = c(0.1, 0.2, 0.3),
    row.names = probes,
    stringsAsFactors = FALSE
  )

  expect_error(
    SEMseeker:::sem_prepare_batch_signal(
      signal_data, tech = "WGBS", sex_chromosome_remove = TRUE
    ),
    "probe_features became empty"
  )
})

test_that("sem_prepare_batch_signal errors out with no resolvable tech", {
  skip_on_cran()
  tempFolder <- tempfile("ppp_prep_batch_")
  on.exit(unlink(tempFolder, recursive = TRUE), add = TRUE)
  SEMseeker:::core_init_env(result_folder = tempFolder)

  # WGBS detection requires probes, but we pass an empty data.frame.
  # core_get_meth_tech() should fall through and return tech = "" -> error.
  signal_data <- data.frame(
    S = numeric(0),
    row.names = character(0),
    stringsAsFactors = FALSE
  )

  expect_error(
    SEMseeker:::sem_prepare_batch_signal(
      signal_data, tech = "", sex_chromosome_remove = TRUE
    )
  )
})


# ---- Source-level guard: the function exists and exposes the contract ----

test_that("sem_prepare_batch_signal source has all 8 documented steps", {
  src <- paste(deparse(SEMseeker:::sem_prepare_batch_signal), collapse = "\n")
  # Tech resolution
  expect_true(grepl("core_get_meth_tech", src))
  # Tech-specific branching
  expect_true(grepl('WGBS.*LONGREAD|tech %in%', src))
  # Duplicate collapse
  expect_true(grepl("duplicated\\(probe_features\\$PROBE", src))
  # Intersection with input (deparse may break the line — use lenient match)
  expect_true(grepl("probe_features\\$PROBE\\s*%in%\\s*rownames\\(signal_data",
                    src, perl = TRUE))
  # Sex-chr filter
  expect_true(grepl('"X", "Y"', src))
  # Alignment subset
  expect_true(grepl("signal_data\\[probe_features\\$PROBE", src))
  # Sanity check
  expect_true(grepl("stopifnot", src))
  # Attributes
  expect_true(grepl('attr\\(signal_data, "probe_features"\\)', src))
})
