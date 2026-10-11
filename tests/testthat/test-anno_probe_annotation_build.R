## Coverage for previously-untested annotation column-builders.
## Both are pure: gene columns parse ";"-delimited manifest strings, and chr
## columns look CpGs up in a cytoband table (injected here for determinism).
##
## Covered:
##   .anno_gene_columns()   UCSC_RefGene_Group/Name -> GENE_<region> columns
##   .anno_chr_columns()    CHR + START -> CHR_CYTOBAND via a cytoband table

# ---------------------------------------------------------------------------
# .anno_gene_columns
# ---------------------------------------------------------------------------

test_that(".anno_gene_columns splits gene region/name pairs into GENE_* columns", {
  out <- SEMseeker:::.anno_gene_columns(
    group_str = c("Body;TSS200", "TSS1500"),
    name_str  = c("TP53;BRCA1",  "EGFR"))

  expect_equal(out$GENE_BODY,    c("TP53", NA))
  expect_equal(out$GENE_TSS200,  c("BRCA1", NA))
  expect_equal(out$GENE_TSS1500, c(NA, "EGFR"))
  # WHOLE is the union of all gene names on each probe
  expect_equal(out$GENE_WHOLE,   c("TP53;BRCA1", "EGFR"))
})

test_that(".anno_gene_columns yields NA for probes with no gene annotation", {
  out <- SEMseeker:::.anno_gene_columns(group_str = "", name_str = "")
  expect_true(is.na(out$GENE_BODY))
  expect_true(is.na(out$GENE_WHOLE))
})

# ---------------------------------------------------------------------------
# .anno_chr_columns
# ---------------------------------------------------------------------------

test_that(".anno_chr_columns maps positions to cytobands and NA outside any band", {
  cb <- data.frame(
    CHR      = c("1", "1"),
    START    = c(1L, 1000L),
    END      = c(999L, 2000L),
    CYTOBAND = c("p1", "p2"),
    stringsAsFactors = FALSE)

  res <- SEMseeker:::.anno_chr_columns(
    chr      = c("1", "1", "1"),
    start    = c(500L, 1500L, 3000L),   # in p1, in p2, past the last band
    cytoband = cb)

  expect_equal(res$CHR_CYTOBAND, c("p1", "p2", NA))
})

test_that(".anno_chr_columns returns NA for chromosomes absent from the table", {
  cb <- data.frame(CHR = "1", START = 1L, END = 999L, CYTOBAND = "p1",
                   stringsAsFactors = FALSE)
  res <- SEMseeker:::.anno_chr_columns(chr = "7", start = 500L, cytoband = cb)
  expect_true(is.na(res$CHR_CYTOBAND))
})
