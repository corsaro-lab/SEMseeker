# Unit tests for the shared CpG-island / OpenSea semantics (island_opensea.R).
# Pure GRanges logic — no Bioconductor annotation packages or downloads.

testthat::skip_if_not_installed("GenomicRanges")
testthat::skip_if_not_installed("IRanges")

test_that(".anno_islands_gr_from_names parses unique coordinate strings, ignores junk", {
  gr <- SEMseeker:::.anno_islands_gr_from_names(
    c("chr1:1000-2000", "chr1:1000-2000", "chr2:5000-6000", NA, "", "OpenSea"))
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2L)  # dedup + drop NA/""/non-coordinate
  expect_setequal(as.character(GenomicRanges::seqnames(gr)), c("chr1", "chr2"))
})

test_that(".anno_opensea_gaps returns inter-neighbourhood gaps labelled by coordinate", {
  islands <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = c(10000L, 30000L), end = c(11000L, 31000L)))

  gaps_gr <- SEMseeker:::.anno_opensea_gaps(islands, chrom_ends = list(chr1 = 50000L))

  # neighbourhoods: [6000,15000] and [26000,35000] -> gaps:
  #   [1,5999], [15001,25999], [35001,50000]
  expect_true(all(as.character(GenomicRanges::strand(gaps_gr)) == "*"))
  labs <- GenomicRanges::mcols(gaps_gr)$label
  expect_true("chr1:15001-25999" %in% labs)
  expect_true(all(grepl("^chr1:[0-9]+-[0-9]+$", labs)))
  # No gap crosses into a neighbourhood
  expect_false(any(GenomicRanges::start(gaps_gr) >= 6000 &
                   GenomicRanges::end(gaps_gr)   <= 15000))
})

test_that(".anno_assign_opensea_labels maps each open-sea CpG to its gap; NA inside islands", {
  islands <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = c(10000L, 30000L), end = c(11000L, 31000L)))

  labels <- SEMseeker:::.anno_assign_opensea_labels(
    probe_chr = rep("chr1", 4L),
    probe_pos = c(1000L, 20000L, 40000L, 10500L),
    island_gr = islands  # no seqlengths: helper extends universe to last CpG
  )

  expect_true(grepl("^chr1:[0-9]+-[0-9]+$", labels[1]))  # before first island
  expect_equal(labels[2], "chr1:15001-25999")            # between islands
  expect_true(grepl("^chr1:35001-", labels[3]))          # after last island
  expect_true(is.na(labels[4]))                          # inside island core
})

# --- Integration glue: the same paths anno_probe_annotation_build / anno_area_granges_build
# use, exercised without an Illumina annotation package or AnnotationHub. -------

test_that(".anno_island_columns recodes all 6 Relation_to_Island categories + OPENSEA", {
  rel   <- c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf",
             "OpenSea", "OpenSea")
  name  <- c(rep("chr1:10000-11000", 5L), "", "")  # OpenSea has empty Islands_Name
  chr   <- rep("1", 7L)                             # CHR is stored without "chr"
  start <- c(10500L, 9000L, 11500L, 7000L, 13500L, 1000L, 40000L)

  cols <- SEMseeker:::.anno_island_columns(rel, name, chr, start)

  expect_named(cols, c("ISLAND_WHOLE", "ISLAND_ISLAND", "ISLAND_N_SHORE",
    "ISLAND_S_SHORE", "ISLAND_N_SHELF", "ISLAND_S_SHELF", "ISLAND_OPENSEA"))

  # WHOLE = whole neighbourhood: all 5 island-context probes, NA for open-sea.
  expect_equal(cols$ISLAND_WHOLE,
    c(rep("chr1:10000-11000", 5L), NA_character_, NA_character_))
  # Core / shores / shelves: one probe each, NA elsewhere.
  expect_equal(cols$ISLAND_ISLAND[1], "chr1:10000-11000")
  expect_true(all(is.na(cols$ISLAND_ISLAND[-1])))
  expect_equal(cols$ISLAND_N_SHORE[2], "chr1:10000-11000")
  expect_equal(cols$ISLAND_S_SHELF[5], "chr1:10000-11000")
  # OPENSEA: only the two open-sea probes get a gap coordinate; core neighbourhood
  # is core +/- 4kb = [6000,15000], so probe@1000 -> gap [1,5999].
  expect_true(all(is.na(cols$ISLAND_OPENSEA[1:5])))
  expect_equal(cols$ISLAND_OPENSEA[6], "chr1:1-5999")
  expect_match(cols$ISLAND_OPENSEA[7], "^chr1:15001-")
})

test_that(".anno_build_island_area maps subareas with injected islands (no AnnotationHub)", {
  islands <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = c(10000L, 30000L), end = c(11000L, 31000L)))
  lab <- function(gr) GenomicRanges::mcols(gr)$label

  core <- SEMseeker:::.anno_build_island_area("ISLAND", "hg19", islands = islands)
  expect_equal(GenomicRanges::start(core), c(10000, 30000))
  expect_true(all(grepl("^chr1:", lab(core))))

  whole <- SEMseeker:::.anno_build_island_area("WHOLE", "hg19", islands = islands)
  expect_equal(GenomicRanges::start(whole), c(6000, 26000))   # core +/- 4kb
  expect_equal(GenomicRanges::end(whole),   c(15000, 35000))

  nsh <- SEMseeker:::.anno_build_island_area("N_SHORE", "hg19", islands = islands)
  expect_equal(length(nsh), 2L)

  expect_error(
    SEMseeker:::.anno_build_island_area("FOOBAR", "hg19", islands = islands),
    regexp = "Unknown ISLAND subarea")
})
# NOTE: .anno_build_island_area("OPENSEA") is intentionally NOT unit-tested here.
# Its only logic beyond the shared .anno_opensea_gaps() (already covered above) is a
# GenomeInfoDb::seqlengths() read, and GenomeInfoDb must not appear in a
# skip_if_not_installed() (test-pkg_suggests_installed.R forbids skip-guarding a
# package that is not in the CI install list). The OpenSea gap computation is
# fully exercised by the .anno_opensea_gaps / .anno_assign_opensea_labels tests.
