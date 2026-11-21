## Coverage for previously-untested annotation column-builders.
## Both are pure: gene columns parse ";"-delimited manifest strings, and chr
## columns look CpGs up in a cytoband table (injected here for determinism).
##
## Covered:
##   .anno_gene_columns()   UCSC_RefGene_Group/Name -> GENE_<region> columns

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

