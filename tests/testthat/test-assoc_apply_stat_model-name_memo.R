# 2026-06-09: guard the sanitize+memo+counter-rename
# pattern. The contract is:
#
#   1. The on-disk CSV writes AREA_OF_TEST as the RAW upstream-annotated
#      name (HLA-A, chr10:100028204-100028508, …).
#   2. Inside the formula machinery, the colname is a R-safe form
#      (HLA_A, chr10_100028204_100028508, …) so `lm(y ~ x)` parses.
#   3. The mapping is 1-to-1 even when sanitisation would collide
#      (HLA-A vs HLA_A both → HLA_A) — `make.unique()` disambiguates.
#
# These tests cover the pattern at the UNIT level. The downstream
# integration (full assoc_apply_stat_model run + result CSV inspection) is
# left to higher-level smoke tests because of the foreach + multisession
# setup overhead.

# ---- Pattern: sanitize ----------------------------------------------------

test_that("sanitize converts every non-R-safe character to underscore", {
  real <- c("HLA-A",
            "chr10:100028204-100028508",
            "TP53",
            "ANKHD1-EIF4EBP3",
            "gene/with/slash")
  safe <- gsub("[^A-Za-z0-9_.]", "_", real)

  # Every output character is in the R-identifier-safe set
  for (s in safe) {
    expect_true(grepl("^[A-Za-z0-9_.]+$", s),
                info = sprintf("sanitised name '%s' is not R-safe", s))
  }

  # Spot-checks
  expect_equal(safe[1], "HLA_A")
  expect_equal(safe[2], "chr10_100028204_100028508")
  expect_equal(safe[3], "TP53")                       # unchanged
  expect_equal(safe[4], "ANKHD1_EIF4EBP3")
  expect_equal(safe[5], "gene_with_slash")
})

# ---- Pattern: memo + reverse lookup --------------------------------------

test_that("safe_to_real reverse map recovers the raw name verbatim", {
  real <- c("HLA-A",
            "chr10:100028204-100028508",
            "TP53",
            "Sample_ID")
  safe <- gsub("[^A-Za-z0-9_.]", "_", real)
  safe_to_real <- setNames(real, safe)

  # Lookup: safe form → raw form
  expect_equal(safe_to_real[["HLA_A"]],                    "HLA-A")
  expect_equal(safe_to_real[["chr10_100028204_100028508"]],"chr10:100028204-100028508")
  expect_equal(safe_to_real[["TP53"]],                     "TP53")
  expect_equal(safe_to_real[["Sample_ID"]],                "Sample_ID")
})

# ---- Pattern: duplicate disambiguation ------------------------------------

test_that("make.unique disambiguates collisions after sanitisation", {
  # Adversarial: distinct raw names that collapse to the same safe form
  real <- c("HLA-A", "HLA_A", "HLA.A")
  safe <- gsub("[^A-Za-z0-9_.]", "_", real)
  # First two collapse to "HLA_A", third stays "HLA.A"
  # After make.unique: "HLA_A", "HLA_A_1", "HLA.A"
  if (anyDuplicated(safe)) {
    safe <- make.unique(safe, sep = "_")
  }
  expect_equal(length(unique(safe)), length(real))
  safe_to_real <- setNames(real, safe)

  # Each raw is still reachable through its (now unique) safe key
  expect_equal(safe_to_real[[safe[1]]], "HLA-A")
  expect_equal(safe_to_real[[safe[2]]], "HLA_A")
  expect_equal(safe_to_real[[safe[3]]], "HLA.A")
})

# ---- Source-level guards on the surface -------------------------
# These tests inspect deparse() of the installed functions to guarantee
# that the gsub("-","_") (or its polars equivalent) does NOT come back
# silently in a future refactor. They are intentionally NOT runtime tests
# — the full bulk path needs a non-trivial synthetic setup, which is
# brittle and unrelated to the contract under test.

test_that("assoc_apply_stat_model_batch_lazy does NOT '-'→'_' rewrite AREA (bulk path)", {
  src <- paste(deparse(SEMseeker:::assoc_apply_stat_model_batch_lazy), collapse = "\n")
  # Polars-side normalisation that was removed:
  expect_false(
    grepl('str\\$replace_all\\(\\s*"-"\\s*,\\s*"_"\\s*\\)', src),
    info = "the polars '-'→'_' on AREA was removed: keep it removed"
  )
  # R-side: should not gsub on AREA either:
  expect_false(
    grepl('gsub\\(\\s*"-"\\s*,\\s*"_"\\s*,\\s*[^)]*AREA[^)]*\\)', src),
    info = "the R-side '-'→'_' on AREA was removed: keep it removed"
  )
})

test_that("assoc_run_marker does NOT '-'→'_' rewrite AREA before resume match", {
  src <- paste(deparse(SEMseeker:::assoc_run_marker), collapse = "\n")
  expect_false(
    grepl('gsub\\(\\s*"-"\\s*,\\s*"_"\\s*,\\s*tempDataFrame\\$AREA\\s*\\)', src),
    info = "the AREA rewrite was removed: keep it removed"
  )
})

test_that("io_data_preparation does NOT '-'→'_' rewrite tempDataFrame colnames", {
  src <- paste(deparse(SEMseeker:::io_data_preparation), collapse = "\n")
  expect_false(
    grepl('gsub\\(\\s*"-"\\s*,\\s*"_"\\s*,\\s*colnames\\(tempDataFrame\\)\\s*\\)', src),
    info = "io_data_preparation's colname sanitisation was removed: keep it removed"
  )
})

test_that("assoc_apply_stat_model carries the sanitize+memo+counter-rename pattern (per-gene)", {
  src <- paste(deparse(SEMseeker:::assoc_apply_stat_model), collapse = "\n")
  # Sanitize step
  expect_true(
    grepl('safe_cols\\s*<-\\s*gsub\\(\\s*"\\[\\^A-Za-z0-9_\\.\\]"\\s*,\\s*"_"\\s*,\\s*real_cols\\s*\\)', src),
    info = "Per-gene path must build safe_cols from real_cols"
  )
  # Memo
  expect_true(
    grepl("safe_to_real\\s*<-\\s*(stats::)?setNames\\(real_cols,\\s*safe_cols\\)", src),
    info = "Per-gene path must memoise safe_to_real"
  )
  # Counter-rename when assigning AREA_OF_TEST
  expect_true(
    grepl("safe_to_real\\[\\[\\s*burdenValue\\s*\\]\\]", src),
    info = "Per-gene path must reverse-map AREA_OF_TEST via safe_to_real"
  )
})

test_that("diagnostic_performance carries the same sanitize+memo+counter-rename pattern", {
  src <- paste(deparse(SEMseeker:::diagnostic_performance), collapse = "\n")
  expect_true(
    grepl('safe_cols\\s*<-\\s*gsub\\(\\s*"\\[\\^A-Za-z0-9_\\.\\]"\\s*,\\s*"_"\\s*,\\s*real_cols\\s*\\)', src),
    info = "diagnostic_performance must build safe_cols from real_cols"
  )
  expect_true(
    grepl("safe_to_real\\s*<-\\s*(stats::)?setNames\\(real_cols,\\s*safe_cols\\)", src),
    info = "diagnostic_performance must memoise safe_to_real"
  )
  expect_true(
    grepl("safe_to_real\\[\\[\\s*area_of_test\\s*\\]\\]", src),
    info = "diagnostic_performance must reverse-map area_of_test via safe_to_real"
  )
})
