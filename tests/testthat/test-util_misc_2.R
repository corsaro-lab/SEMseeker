# Coverage for previously-untested PURE / UTILITY helpers (no session needed).
#
# These are deterministic, dependency-free functions (label formatting, order
# checks, gene-set set-algebra, string wrapping, BLAS introspection). Assertions
# are exact where the output is fully determined; the annotation-package helpers
# skip_if_not_installed() because they read Bioconductor Suggests data.
#
# Covered:
#   util_pretty_label()                     underscore->space display labels
#   util_test_match_order()                 same-order equality of two vectors
#   util_test_match_order_by_rownames()     same-rowname-order of two frames
#   enrich_find_unique_gene_sets()          per-group unique set-difference
#   enrich_wrap_it()                        fixed-width string wrapping
#   .core_blas_info()                       BLAS flavour introspection
#   .core_warn_blas_single_thread()         one-shot BLAS warning (no-op log)
#   .anno_detect_tech_from_anno()           array detection by probe overlap
#   .anno_pkg_load_table/_probe_ids/_to_df  Illumina anno-package readers

# ---------------------------------------------------------------------------
# util_pretty_label
# ---------------------------------------------------------------------------

test_that("util_pretty_label turns UPPER_SNAKE into spaced labels", {
  expect_equal(SEMseeker:::util_pretty_label("TUMOUR_STAGE_N"),   "TUMOUR STAGE N")
  expect_equal(SEMseeker:::util_pretty_label("PVALUE_ADJ_ALL_FDR"), "PVALUE ADJ ALL FDR")
})

test_that("util_pretty_label collapses repeated separators and trims", {
  expect_equal(SEMseeker:::util_pretty_label("A__B"), "A B")     # double underscore
  expect_equal(SEMseeker:::util_pretty_label("_A_"),  "A")       # leading/trailing
})

test_that("util_pretty_label passes NA and empty through, preserves length", {
  expect_equal(SEMseeker:::util_pretty_label(c("DELTARP_GENE_TSS200", NA, "")),
               c("DELTARP GENE TSS200", NA, ""))
  expect_equal(SEMseeker:::util_pretty_label(character(0)), character(0))
})

# ---------------------------------------------------------------------------
# util_test_match_order / util_test_match_order_by_rownames
# ---------------------------------------------------------------------------

test_that("util_test_match_order is TRUE only for identical same-order vectors", {
  expect_true (SEMseeker:::util_test_match_order(c(1, 2, 3), c(1, 2, 3)))
  expect_true (SEMseeker:::util_test_match_order(c("a", "b"), c("a", "b")))
  expect_false(SEMseeker:::util_test_match_order(c(1, 2, 3), c(3, 2, 1)))  # same set, wrong order
  expect_false(SEMseeker:::util_test_match_order(c(1, 2, 3), c(4, 5, 6)))  # disjoint
})

test_that("util_test_match_order returns FALSE on NULL input", {
  expect_false(SEMseeker:::util_test_match_order(NULL, c(1, 2, 3)))
  expect_false(SEMseeker:::util_test_match_order(c(1, 2, 3), NULL))
})

test_that("util_test_match_order_by_rownames compares rowname order", {
  a <- data.frame(v = 1:3, row.names = c("p1", "p2", "p3"))
  b <- data.frame(v = 4:6, row.names = c("p1", "p2", "p3"))
  c <- data.frame(v = 7:9, row.names = c("p3", "p2", "p1"))
  expect_true (SEMseeker:::util_test_match_order_by_rownames(a, b))
  expect_false(SEMseeker:::util_test_match_order_by_rownames(a, c))
  expect_false(SEMseeker:::util_test_match_order_by_rownames(NULL, a))
})

# ---------------------------------------------------------------------------
# enrich_find_unique_gene_sets
# ---------------------------------------------------------------------------

test_that("enrich_find_unique_gene_sets keeps only per-group unique members", {
  res <- SEMseeker:::enrich_find_unique_gene_sets(
    list(A = c("g1", "g2", "g3"), B = c("g3", "g4"), C = c("g5")))
  expect_equal(sort(res$A), c("g1", "g2"))
  expect_equal(res$B, "g4")
  expect_equal(res$C, "g5")
})

test_that("enrich_find_unique_gene_sets drops groups with no unique members", {
  res <- SEMseeker:::enrich_find_unique_gene_sets(
    list(A = c("g1"), B = c("g1")))
  expect_length(res, 0)
})

# ---------------------------------------------------------------------------
# enrich_wrap_it
# ---------------------------------------------------------------------------

test_that("enrich_wrap_it wraps long strings and leaves short ones intact", {
  expect_equal(SEMseeker:::enrich_wrap_it("hello", 20), "hello")           # fits
  wrapped <- SEMseeker:::enrich_wrap_it("aaaa bbbb cccc dddd", 9)
  expect_match(wrapped, "\n")                                              # broke into lines
  expect_equal(length(SEMseeker:::enrich_wrap_it(c("a", "b"), 5)), 2L)     # length preserved
  expect_null(names(SEMseeker:::enrich_wrap_it(c("a", "b"), 5)))           # USE.NAMES = FALSE
})

# ---------------------------------------------------------------------------
# .core_blas_info / .core_warn_blas_single_thread
# ---------------------------------------------------------------------------

test_that(".core_blas_info returns a coherent BLAS descriptor", {
  info <- SEMseeker:::.core_blas_info()
  expect_named(info, c("path", "multi_threaded", "flavor", "is_reference"),
               ignore.order = TRUE)
  expect_type(info$multi_threaded, "logical")
  expect_type(info$is_reference,   "logical")
  expect_true(info$flavor %in% c("Accelerate (vecLib)", "OpenBLAS", "MKL",
                                 "reference (single-thread)", "unknown"))
  # a reference BLAS is single-threaded by definition
  if (isTRUE(info$is_reference)) {
    expect_equal(info$flavor, "reference (single-thread)")
    expect_false(info$multi_threaded)
  }
})

test_that(".core_warn_blas_single_thread returns the BLAS info invisibly", {
  # core_log_event() no-ops without an initialised session, so this is safe.
  expect_no_error(res <- SEMseeker:::.core_warn_blas_single_thread())
  expect_equal(res$flavor, SEMseeker:::.core_blas_info()$flavor)
})

# ---------------------------------------------------------------------------
# .anno_detect_tech_from_anno  (deterministic no-overlap path)
# ---------------------------------------------------------------------------

test_that(".anno_detect_tech_from_anno returns '' when nothing overlaps", {
  # No real probe id matches -> every array count is 0 (or no anno pkg) -> "".
  expect_identical(SEMseeker:::.anno_detect_tech_from_anno("ZZZ_not_a_probe"), "")
  expect_identical(SEMseeker:::.anno_detect_tech_from_anno(character(0)),      "")
})

# ---------------------------------------------------------------------------
# .anno_pkg_* readers  (require an Illumina annotation Suggests package)
# ---------------------------------------------------------------------------

test_that(".anno_pkg_probe_ids / .anno_pkg_load_table expose Illumina probe ids", {
  pkg <- "IlluminaHumanMethylation27kanno.ilmn12.hg19"   # smallest array
  skip_if_not_installed(pkg)

  ids <- SEMseeker:::.anno_pkg_probe_ids(pkg)
  expect_type(ids, "character")
  expect_gt(length(ids), 0L)
  expect_true(any(grepl("^cg", ids)))

  locs <- SEMseeker:::.anno_pkg_load_table(pkg, "Locations")
  expect_s3_class(locs, "data.frame")
  expect_true("chr" %in% colnames(locs))
  expect_setequal(rownames(locs), ids)
})

test_that(".anno_pkg_to_df combines Locations + Islands + Other for one probe/row", {
  # Needs an array that ships the Islands.UCSC sub-table (450k / EPIC; the 27k
  # array lacks it and the current fallback in .anno_pkg_to_df cannot rebuild a
  # matching-length frame -- tracked separately, not exercised here).
  pkg <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  skip_if_not_installed(pkg)

  df  <- SEMseeker:::.anno_pkg_to_df(pkg)
  ids <- SEMseeker:::.anno_pkg_probe_ids(pkg)
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0L)
  expect_true(all(c("chr", "Relation_to_Island") %in% colnames(df)))
  expect_setequal(rownames(df), ids)
})
