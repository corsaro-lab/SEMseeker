## Tests for anno_area_granges_build() — C-04
##
## All tests that require TxDb / AnnotationHub are guarded with
## skip_if_not_installed(). CI installs these packages explicitly, but they
## are optional for package users.

# ---------------------------------------------------------------------------
# Tests that do NOT require any external Bioconductor package
# ---------------------------------------------------------------------------

test_that("anno_area_granges_build errors on unknown area", {
  expect_error(
    SEMseeker:::anno_area_granges_build("FOOBAR_WHOLE", genome_build = "hg19"),
    regexp = "Unknown area"
  )
})

test_that("anno_area_granges_build errors on unknown GENE subarea", {
  # Bioc anno pkgs trigger requireNamespace -> minfi -> GEOquery -> tcltk segfault on R 4.6 arm64 macOS
  skip_on_os("mac")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomicFeatures")
  expect_error(
    SEMseeker:::anno_area_granges_build("GENE_FOOBAR", genome_build = "hg19"),
    regexp = "Unknown GENE subarea"
  )
})

test_that("anno_area_granges_build errors on unknown ISLAND subarea", {
  # Bioc anno pkgs trigger requireNamespace -> minfi -> GEOquery -> tcltk segfault on R 4.6 arm64 macOS
  skip_on_os("mac")
  skip_if_not_installed("AnnotationHub")
  skip_if_not_installed("GenomicRanges")
  expect_error(
    SEMseeker:::anno_area_granges_build("ISLAND_FOOBAR", genome_build = "hg19"),
    regexp = "Unknown ISLAND subarea"
  )
})

test_that("anno_area_granges_build appends _WHOLE when no underscore", {
  # Bioc anno pkgs trigger requireNamespace -> minfi -> GEOquery -> tcltk segfault on R 4.6 arm64 macOS
  skip_on_os("mac")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomicFeatures")
  # "GENE" alone should be treated as "GENE_WHOLE"
  # We just check it doesn't throw a "unknown area" error
  expect_error(
    SEMseeker:::anno_area_granges_build("GENE", genome_build = "hg19"),
    NA   # no error expected
  )
})

# ---------------------------------------------------------------------------
# DMR area — uses bundled data, no external packages needed
# ---------------------------------------------------------------------------

test_that("anno_area_granges_build DMR_WHOLE returns GRanges with label", {
  # Bioc anno pkgs trigger requireNamespace -> minfi -> GEOquery -> tcltk segfault on R 4.6 arm64 macOS
  skip_on_os("mac")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  gr <- SEMseeker:::anno_area_granges_build("DMR_WHOLE", genome_build = "hg19")
  expect_s4_class(gr, "GRanges")
  expect_true(length(gr) > 0)
  expect_true("label" %in% names(GenomicRanges::mcols(gr)))
  expect_false(any(is.na(GenomicRanges::mcols(gr)$label)))
})

# ---------------------------------------------------------------------------
# GENE areas — require TxDb
# ---------------------------------------------------------------------------

test_that("anno_area_granges_build GENE_BODY returns valid GRanges", {
  # Bioc anno pkgs trigger requireNamespace -> minfi -> GEOquery -> tcltk segfault on R 4.6 arm64 macOS
  skip_on_os("mac")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomicFeatures")
  gr <- SEMseeker:::anno_area_granges_build("GENE_BODY", genome_build = "hg19")
  expect_s4_class(gr, "GRanges")
  expect_true(length(gr) > 1000L)
  expect_true("label" %in% names(GenomicRanges::mcols(gr)))
})

test_that("anno_area_granges_build GENE_TSS200 ranges are ≤ 200bp wide", {
  # Bioc anno pkgs trigger requireNamespace -> minfi -> GEOquery -> tcltk segfault on R 4.6 arm64 macOS
  skip_on_os("mac")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomicFeatures")
  gr <- SEMseeker:::anno_area_granges_build("GENE_TSS200", genome_build = "hg19")
  expect_true(all(GenomicRanges::width(gr) <= 200L))
})

test_that("GENE_TSS1500 and GENE_TSS200 do not overlap within the same gene", {
  # Bioc anno pkgs trigger requireNamespace -> minfi -> GEOquery -> tcltk segfault on R 4.6 arm64 macOS
  skip_on_os("mac")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomicFeatures")
  gr200  <- SEMseeker:::anno_area_granges_build("GENE_TSS200",  genome_build = "hg19")
  gr1500 <- SEMseeker:::anno_area_granges_build("GENE_TSS1500", genome_build = "hg19")
  # Per-gene check: TSS200 and TSS1500 must not overlap within the same gene
  # (they are adjacent rings by design). Cross-gene overlaps between
  # neighbouring genes are expected and do not violate the ring structure.
  lbl200  <- GenomicRanges::mcols(gr200)$label
  lbl1500 <- GenomicRanges::mcols(gr1500)$label
  common  <- intersect(lbl200, lbl1500)
  skip_if(length(common) == 0L, "no gene label overlap between TSS200/TSS1500")
  per_gene_overlap <- vapply(utils::head(common, 500L), function(sym) {
    a <- gr200[lbl200 == sym]
    b <- gr1500[lbl1500 == sym]
    length(GenomicRanges::findOverlaps(a, b))
  }, integer(1L))
  expect_equal(sum(per_gene_overlap), 0L,
    info = "Per-gene: TSS200 and TSS1500 rings must not overlap within the same gene")
})

test_that("anno_area_granges_build result is cached on second call", {
  # Bioc anno pkgs trigger requireNamespace -> minfi -> GEOquery -> tcltk segfault on R 4.6 arm64 macOS
  skip_on_os("mac")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomicFeatures")
  gr1 <- SEMseeker:::anno_area_granges_build("GENE_BODY", genome_build = "hg19")
  gr2 <- SEMseeker:::anno_area_granges_build("GENE_BODY", genome_build = "hg19")
  expect_identical(gr1, gr2)  # exact same object from cache
})

# ---------------------------------------------------------------------------
# CHR_CYTOBAND — uses bundled data
# ---------------------------------------------------------------------------

test_that("anno_area_granges_build CHR_CYTOBAND returns GRanges with label", {
  # Bioc anno pkgs trigger requireNamespace -> minfi -> GEOquery -> tcltk segfault on R 4.6 arm64 macOS
  skip_on_os("mac")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
  gr <- SEMseeker:::anno_area_granges_build("CHR_CYTOBAND", genome_build = "hg19")
  expect_s4_class(gr, "GRanges")
  expect_true(length(gr) > 100L)
  expect_true("label" %in% names(GenomicRanges::mcols(gr)))
  # Labels should look like cytoband names (e.g. "p11.1", "q21.3")
  expect_true(any(grepl("[pq]", GenomicRanges::mcols(gr)$label)))
})

# ---------------------------------------------------------------------------
# GENE subareas whose label travels through names(gr)
#
# BODY, TSS200 and TSS1500 read the per-gene `symbols` vector that
# .anno_gene_granges() computes once from the gene ids, and the tests above
# cover them. 1STEXON, 5UTR, 3UTR and EXONBND take the other branch of the
# label block: they read names(gr) and map those through
# .anno_entrez_to_symbol(). Nothing exercised that branch, and it has two
# failure modes that are silent rather than loud:
#
#   * names(gr) is absent, and the label falls back to the positional
#     placeholder GENE_<SUBAREA>_<n>, which no consumer can join on;
#   * names(gr) holds identifiers of the wrong kind - a transcript id where an
#     Entrez id is expected - and the mapper, finding no match, hands the id
#     straight back, so the label is an identifier that names no gene.
#
# Either way the area is built, the run completes, and the result names genes
# that were never resolved.
#
# The assertion compares the label vocabulary against GENE_BODY, which is one
# range per gene of the same TxDb and is already covered above. Both sides go
# through the same mapper, so the comparison holds whether or not the symbol
# database is installed: what it measures is whether the two windows are
# naming genes of the same model, not how those genes are spelled.
# ---------------------------------------------------------------------------

for (.subarea in c("1STEXON", "5UTR", "3UTR", "EXONBND")) {

  test_that(paste0("anno_area_granges_build GENE_", .subarea,
                   " labels name genes of the same model as GENE_BODY"), {
    # Bioc anno pkgs trigger requireNamespace -> minfi -> GEOquery -> tcltk segfault on R 4.6 arm64 macOS
    skip_on_os("mac")
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("GenomicFeatures")

    area_subarea <- paste0("GENE_", .subarea)
    gr   <- SEMseeker:::anno_area_granges_build(area_subarea,
                                                genome_build = "hg19")
    body <- SEMseeker:::anno_area_granges_build("GENE_BODY",
                                                genome_build = "hg19")

    expect_s4_class(gr, "GRanges")
    expect_true(length(gr) > 0L)
    expect_true("label" %in% names(GenomicRanges::mcols(gr)))

    lbl <- unique(as.character(GenomicRanges::mcols(gr)$label))
    lbl <- lbl[!is.na(lbl) & nzchar(lbl)]
    expect_true(length(lbl) > 0L)

    expect_false(
      any(grepl(paste0("^", area_subarea, "_[0-9]+$"), lbl)),
      info = paste0(area_subarea, ": labels fell back to the positional ",
                    "placeholder, so names(gr) carried no parent gene"))

    ref    <- unique(as.character(GenomicRanges::mcols(body)$label))
    shared <- length(intersect(lbl, ref)) / length(lbl)
    expect_gt(shared, 0.5)
  })
}
