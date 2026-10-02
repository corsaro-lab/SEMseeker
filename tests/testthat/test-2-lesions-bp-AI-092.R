# AI-092: LESIONS clustering by genomic distance (bp) — canary tests
# verifying that LESIONS_BP behaves as the bp-distance threshold for cluster
# membership, with the expected monotonic relationship: smaller LESIONS_BP
# isolates the dense local cluster; larger LESIONS_BP merges spatially
# distant mutations into a single window.

.bp_make_mutations_df <- function(chr, starts, mutations) {
  data.frame(
    CHR       = rep(chr, length(starts)),
    START     = as.integer(starts),
    END       = as.integer(starts) + 1L,
    MUTATIONS = as.integer(mutations),
    stringsAsFactors = FALSE
  )
}

test_that("LESIONS_BP=500 finds 2 separate clusters at 500bp resolution", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(
    result_folder        = tf,
    start_fresh          = TRUE,
    LESIONS_BP           = 500L,
    bonferroni_threshold = 5     # loose so the tight clusters always fire
  )
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # Two tight clusters of mutated probes, 50kb apart (well outside LESIONS_BP=500):
  #   Cluster A: 5 probes at 100bp spacing → span 400bp ≤ 500bp window
  #   Cluster B: 5 probes at 100bp spacing → span 400bp ≤ 500bp window, +50kb
  # Background: 30 isolated unmutated probes scattered between the two.
  pos_A <- 1e6L + seq(0L, by = 100L, length.out = 5L)
  pos_B <- 1e6L + 50000L + seq(0L, by = 100L, length.out = 5L)
  pos_BG <- as.integer(seq(1e6L + 1000L, 1e6L + 49000L, length.out = 30L))
  starts <- c(pos_A, pos_B, pos_BG)
  muts   <- c(rep(1L, 5L), rep(1L, 5L), rep(0L, 30L))
  df <- .bp_make_mutations_df("chr1", starts, muts)

  result <- SEMseeker:::sem_lesions_get(grouping_column = "CHR",
                                    mutation_annotated_sorted = df)
  expect_s3_class(result, "data.frame")
  # With LESIONS_BP=500 the two clusters cannot merge (50kb apart). Each
  # cluster centre sees ~5 co-mutated probes within ±500bp, well above the
  # background, so we expect lesion probes from BOTH clusters.
  expect_gt(nrow(result), 0L)
  lesion_pos <- result$START
  in_A <- any(lesion_pos >= min(pos_A) & lesion_pos <= max(pos_A))
  in_B <- any(lesion_pos >= min(pos_B) & lesion_pos <= max(pos_B))
  expect_true(in_A, info = "expected at least one lesion in cluster A")
  expect_true(in_B, info = "expected at least one lesion in cluster B")
})

test_that("LESIONS_BP=0 reduces every probe to singleton window (no spatial leverage)", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(
    result_folder        = tf,
    start_fresh          = TRUE,
    LESIONS_BP           = 0L,
    bonferroni_threshold = 5
  )
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # 20 mutated + 20 unmutated probes, p0=0.5. With LESIONS_BP=0 every probe
  # is its own window (WINDOW_SIZE=1), so ENRICHMENT is either 0 or 1, and
  # pbinom(0, 1, 0.5, lower.tail=FALSE) = 0.5 — never beats any reasonable
  # Bonferroni threshold. No lesions expected.
  n <- 40L
  starts <- as.integer(seq(1e6L, by = 1000L, length.out = n))
  muts   <- c(rep(1L, 20L), rep(0L, 20L))
  df <- .bp_make_mutations_df("chr1", starts, muts)

  result <- SEMseeker:::sem_lesions_get(grouping_column = "CHR",
                                    mutation_annotated_sorted = df)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0L)
})

test_that("LESIONS_BP=100000 over-dilutes a sparse two-cluster signal — no lesions", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(
    result_folder        = tf,
    start_fresh          = TRUE,
    LESIONS_BP           = 100000L,
    bonferroni_threshold = 5
  )
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # Same two-cluster + background layout as the LESIONS_BP=500 test.
  # With LESIONS_BP=100000 every probe's window catches ALL 40 probes
  # (clusters + background). WINDOW_SIZE=40, ENRICHMENT=10, p0=10/40=0.25.
  # P(X >= 10 | Binom(40, 0.25)) ~= 0.5 — at the mean, never beats any
  # reasonable Bonferroni threshold. This is the CORRECT bp-distance
  # semantics: an over-wide window dilutes local enrichment with global
  # background and kills detection power. The test pins this property so a
  # future regression that re-introduces window-independent thresholding
  # would fail loudly.
  pos_A <- 1e6L + seq(0L, by = 100L, length.out = 5L)
  pos_B <- 1e6L + 50000L + seq(0L, by = 100L, length.out = 5L)
  pos_BG <- as.integer(seq(1e6L + 1000L, 1e6L + 49000L, length.out = 30L))
  starts <- c(pos_A, pos_B, pos_BG)
  muts   <- c(rep(1L, 5L), rep(1L, 5L), rep(0L, 30L))
  df <- .bp_make_mutations_df("chr1", starts, muts)

  result <- SEMseeker:::sem_lesions_get(grouping_column = "CHR",
                                    mutation_annotated_sorted = df)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0L)
})

test_that("sem_lesions_get errors on negative LESIONS_BP", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(
    result_folder = tf,
    start_fresh   = TRUE,
    LESIONS_BP    = -100L
  )
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df <- .bp_make_mutations_df("chr1", seq(1e6L, by = 1000L, length.out = 10L),
                              rep(0L, 10L))
  expect_error(
    SEMseeker:::sem_lesions_get(grouping_column = "CHR",
                            mutation_annotated_sorted = df),
    "LESIONS_BP must be a non-negative integer"
  )
})
