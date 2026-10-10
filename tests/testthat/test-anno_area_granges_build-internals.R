## Coverage for previously-untested internal helpers (annotation + core).
## All deterministic; GRanges-based helpers build tiny synthetic ranges.
##
## Covered:
##   sem_metrics_name_collect()       currently a no-op (returns NULL)
##   .anno_normalize_label_set()      ;-split, chr-strip, upper, sort-unique
##   .anno_area_category()            AREA -> annotation backend switch
##   .anno_islands_gr_from_names()    "chrN:start-end" strings -> GRanges
##   .anno_opensea_gaps()             inter-island gaps as labelled GRanges
##   .anno_assign_opensea_labels()    probes -> containing-gap label / NA
##   .anno_island_columns()           Relation_to_Island -> ISLAND_* columns
##   .core_pop_arg()                  pop a named arg with a default

# ---------------------------------------------------------------------------
# sem_metrics_name_collect  (no-op placeholder)
# ---------------------------------------------------------------------------

test_that("sem_metrics_name_collect is a no-op returning NULL", {
  expect_null(SEMseeker:::sem_metrics_name_collect(data.frame(A = 1)))
})

# ---------------------------------------------------------------------------
# .anno_normalize_label_set
# ---------------------------------------------------------------------------

test_that(".anno_normalize_label_set splits, strips chr, uppercases, dedups, sorts", {
  expect_equal(SEMseeker:::.anno_normalize_label_set(character(0)), character(0))
  expect_equal(
    SEMseeker:::.anno_normalize_label_set(c("TP53;BRCA1", "chr1:100-200")),
    c("1:100-200", "BRCA1", "TP53"))
  # NA / empty dropped, case-insensitive de-duplication
  expect_equal(SEMseeker:::.anno_normalize_label_set(c(NA, "", "tp53", "TP53")), "TP53")
})

# ---------------------------------------------------------------------------
# .anno_area_category
# ---------------------------------------------------------------------------

test_that(".anno_area_category maps AREA to its annotation backend", {
  expect_equal(SEMseeker:::.anno_area_category("GENE_TSS200"),  "txdb")
  expect_equal(SEMseeker:::.anno_area_category("ISLAND_WHOLE"), "annotationhub")
  expect_equal(SEMseeker:::.anno_area_category("CHR_CYTOBAND"), "bundled")
  expect_equal(SEMseeker:::.anno_area_category("DMR_X"),        "bundled")
  expect_equal(SEMseeker:::.anno_area_category("FOO_BAR"),      "unknown")
})

# ---------------------------------------------------------------------------
# .anno_islands_gr_from_names
# ---------------------------------------------------------------------------

test_that(".anno_islands_gr_from_names parses coordinate strings into a GRanges", {
  gr <- SEMseeker:::.anno_islands_gr_from_names(
    c("chr1:100-200", "chrX:5-9", "garbage", NA))
  expect_s4_class(gr, "GRanges")
  expect_length(gr, 2L)
  expect_equal(as.character(GenomicRanges::seqnames(gr)), c("chr1", "chrX"))
  expect_equal(GenomicRanges::start(gr), c(100L, 5L))
  expect_equal(GenomicRanges::end(gr),   c(200L, 9L))

  expect_length(SEMseeker:::.anno_islands_gr_from_names(character(0)), 0L)
  expect_length(SEMseeker:::.anno_islands_gr_from_names(c("nope", "x")), 0L)
})

# ---------------------------------------------------------------------------
# .anno_opensea_gaps  /  .anno_assign_opensea_labels
# ---------------------------------------------------------------------------

test_that(".anno_opensea_gaps returns labelled inter-island gap ranges", {
  islands <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = c(1000L, 50000L), end = c(2000L, 51000L)))
  gaps <- SEMseeker:::.anno_opensea_gaps(islands)
  expect_s4_class(gaps, "GRanges")
  expect_true("label" %in% names(GenomicRanges::mcols(gaps)))
  expect_gte(length(gaps), 1L)   # at least the gap between the two islands
})

test_that(".anno_assign_opensea_labels labels open-sea probes and NA-s island probes", {
  islands <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000L, 2000L))
  out <- SEMseeker:::.anno_assign_opensea_labels(
    probe_chr = c("chr1", "chr1"),
    probe_pos = c(1500L, 30000L),   # inside island vs far away (open sea)
    island_gr = islands)
  expect_length(out, 2L)
  expect_true(is.na(out[1]))        # inside the island neighbourhood
  expect_false(is.na(out[2]))       # lands in an inter-island gap
})

# ---------------------------------------------------------------------------
# .anno_island_columns
# ---------------------------------------------------------------------------

test_that(".anno_island_columns fans Relation_to_Island out into ISLAND_* columns", {
  out <- SEMseeker:::.anno_island_columns(
    island_rel  = c("Island", "N_Shore", "OpenSea"),
    island_name = c("chr1:1000-2000", "chr1:1000-2000", "chr1:1000-2000"),
    chr         = c("1", "1", "1"),
    start       = c(1500L, 900L, 30000L))

  expect_equal(out$ISLAND_ISLAND,  c("chr1:1000-2000", NA, NA))
  expect_equal(out$ISLAND_N_SHORE, c(NA, "chr1:1000-2000", NA))
  expect_true(is.na(out$ISLAND_WHOLE[3]))              # open sea excluded from WHOLE
  expect_length(out$ISLAND_OPENSEA, 3L)
  expect_false(is.na(out$ISLAND_OPENSEA[3]))           # open-sea probe gets a gap label
})

# ---------------------------------------------------------------------------
# .core_pop_arg
# ---------------------------------------------------------------------------

test_that(".core_pop_arg pops a present arg and falls back to the default", {
  present <- SEMseeker:::.core_pop_arg(list(a = 1, b = 2), "a", 99)
  expect_equal(present$value, 1)
  expect_false("a" %in% names(present$args))
  expect_true("b" %in% names(present$args))

  missing <- SEMseeker:::.core_pop_arg(list(b = 2), "z", 99)
  expect_equal(missing$value, 99)
})
