# Unit tests for the pure per-area column helpers used by anno_probe_annotation_build():
# .anno_gene_columns / .anno_chr_columns (ISLAND lives in test-0-island-opensea.R). These
# exercise the annotation recoding without an Illumina annotation package, so the
# glue is covered in CI.
#
# NOTE — DMR is deliberately NOT a column-helper: a probe can belong to several
# DMRs, so it is a row-EXPANDING join (merge), not a 1:1 per-probe column.
# anno_probe_features_get() relies on that expansion (+ distinct()) to preserve
# multi-DMR membership, so DMR stays as merge() inside anno_probe_annotation_build().

test_that(".anno_gene_columns recodes RefGene groups into the 9 GENE_* columns", {
  group <- c("Body", "TSS200;Body", "1stExon", "5'UTR;3'UTR", "")
  name  <- c("BRCA1", "GENEA;GENEB", "GENEC", "GENED;GENED", "")

  g <- SEMseeker:::.anno_gene_columns(group, name)

  expect_named(g, c("GENE_BODY", "GENE_TSS200", "GENE_TSS1500", "GENE_1STEXON",
                    "GENE_5UTR", "GENE_3UTR", "GENE_EXONBND", "GENE_WHOLE",
                    "GENE_PROMOTER"))
  expect_equal(g$GENE_BODY,    c("BRCA1", "GENEB", NA, NA, NA))
  expect_equal(g$GENE_TSS200,  c(NA, "GENEA", NA, NA, NA))
  expect_equal(g$GENE_1STEXON, c(NA, NA, "GENEC", NA, NA))
  expect_equal(g$GENE_5UTR,    c(NA, NA, NA, "GENED", NA))
  expect_equal(g$GENE_3UTR,    c(NA, NA, NA, "GENED", NA))
  # WHOLE = all genes overlapping the probe (deduplicated), like ISLAND_WHOLE
  expect_equal(g$GENE_WHOLE,   c("BRCA1", "GENEA;GENEB", "GENEC", "GENED", NA))
  # PROMOTER = the genes annotated to TSS200, TSS1500 or 1stExon. Body only is not
  # the promoter (probe 1), and neither are the UTRs (probe 4).
  expect_equal(g$GENE_PROMOTER, c(NA, "GENEA", "GENEC", NA, NA))
})

test_that("GENE_PROMOTER names a gene once however many of its windows are hit", {
  # The property that makes the class usable: it is a union over three windows, so
  # a gene reached by two of them is one gene here, not two. Summing the three
  # windows instead would count its positions twice.
  group <- c("TSS200;TSS1500", "TSS1500;1stExon", "TSS200;Body")
  name  <- c("GENEX;GENEX",    "GENEY;GENEZ",     "GENEW;GENEW")

  g <- SEMseeker:::.anno_gene_columns(group, name)

  # one gene through two promoter windows -> named once
  expect_equal(g$GENE_PROMOTER[1], "GENEX")
  # two genes through two different promoter windows -> both named
  expect_equal(g$GENE_PROMOTER[2], "GENEY;GENEZ")
  # a window outside the promoter does not bring its gene in on that account,
  # and a gene that is also in the promoter is still named once
  expect_equal(g$GENE_PROMOTER[3], "GENEW")
  expect_equal(g$GENE_BODY[3], "GENEW")

  # the single-window columns are unchanged by the generalisation of the helper
  expect_equal(g$GENE_TSS200,  c("GENEX", NA, "GENEW"))
  expect_equal(g$GENE_TSS1500, c("GENEX", "GENEY", NA))
})

test_that(".anno_chr_columns assigns cytoband by range overlap (injected table)", {
  cytoband <- data.frame(
    CHR      = c("1", "1"),
    START    = c(1L, 10000L),
    END      = c(9999L, 20000L),
    CYTOBAND = c("p36.33", "p36.32"),
    stringsAsFactors = FALSE)

  out <- SEMseeker:::.anno_chr_columns(
    chr = c("1", "1", "2"), start = c(5000L, 15000L, 500L), cytoband = cytoband)

  expect_named(out, "CHR_CYTOBAND")
  expect_equal(out$CHR_CYTOBAND, c("p36.33", "p36.32", NA))  # chr2 absent -> NA
})
