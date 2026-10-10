# Unit tests for the pure per-area column helpers used by anno_probe_annotation_build():
# .anno_gene_columns / .anno_chr_columns (ISLAND lives in test-anno_island_opensea.R). These
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

# ---- The cache has to know which columns it was written for ----------------
# A class added to the vocabulary adds a column to the annotation. A cache keyed
# only on technology and genome build is reused across that change, so the request
# asks for a column the stored table does not carry and the run stops on an
# undefined column, eleven files away from the change that caused it. Every
# machine that has run the package before the upgrade holds exactly such a file.
#
# These exercise the in-memory tier only: the key is one no file can exist for, so
# nothing is read from or written to the user's cache directory.

test_that("a cached annotation written for other columns is a miss, not an error", {
  env  <- SEMseeker:::.pkgglobalenv
  prev <- env$probe_annotation_memo
  on.exit(assign("probe_annotation_memo", prev, envir = env), add = TRUE)

  key    <- paste0("schema_probe_", as.integer(Sys.time()))
  before <- c("PROBE", "CHR", "START", "END", "K850", "GENE_BODY")
  after  <- c(before, "GENE_PROMOTER")
  table  <- data.frame(PROBE = "cg00000029", CHR = "1", START = 1L, END = 1L,
                       K850 = TRUE, GENE_BODY = NA_character_,
                       stringsAsFactors = FALSE)

  assign("probe_annotation_memo",
         list(key = key, schema = before, data = table), envir = env)

  # Same columns: the cache answers.
  expect_identical(SEMseeker:::.anno_probe_cache_get(key, before), table)

  # One column more: it is a miss, and a miss is NULL rather than a condition,
  # because the annotation can always be rebuilt.
  expect_null(SEMseeker:::.anno_probe_cache_get(key, after))
})

test_that("a cache with no schema recorded is a miss", {
  env  <- SEMseeker:::.pkgglobalenv
  prev <- env$probe_annotation_memo
  on.exit(assign("probe_annotation_memo", prev, envir = env), add = TRUE)

  # What every machine that ran an earlier build holds: the table on its own.
  key   <- paste0("schema_probe_legacy_", as.integer(Sys.time()))
  table <- data.frame(PROBE = "cg00000029", stringsAsFactors = FALSE)
  assign("probe_annotation_memo", list(key = key, data = table), envir = env)

  expect_null(SEMseeker:::.anno_probe_cache_get(key, c("PROBE", "GENE_PROMOTER")))
})

test_that(".anno_probe_schema names the columns the builder selects", {
  s <- SEMseeker:::.anno_probe_schema("K850")
  expect_true(all(c("PROBE", "CHR", "START", "END", "K850") %in% s))
  expect_true("GENE_PROMOTER" %in% s)
  # The technology is part of the schema, so a cache built for one array is not
  # offered to another even before the key is compared.
  expect_false("K450" %in% s)
  expect_true("K450" %in% SEMseeker:::.anno_probe_schema("K450"))
})

test_that("the cache file name carries the shape of what it holds", {
  # The schema inside the file settles one direction: a build reading a file written
  # for other columns rebuilds it. The other direction cannot be settled from
  # inside, because the build that gets it wrong is one that already shipped: it
  # reads the file and uses the result as a table, so a file holding anything else
  # stops it on an undefined column. The defence is the name, which is all it reads
  # before opening.
  f <- basename(SEMseeker:::.anno_probe_cache_file("K850_hg19"))

  expect_match(f, "^probe_annotation_v2_K850_hg19\\.rds$")

  # And explicitly not the name an earlier build goes looking for, which is the
  # whole point: the two never meet.
  expect_false(identical(f, "probe_annotation_K850_hg19.rds"))

  # The key still separates technology and genome build within the format.
  expect_false(identical(f, basename(SEMseeker:::.anno_probe_cache_file("K450_hg19"))))
  expect_false(identical(f, basename(SEMseeker:::.anno_probe_cache_file("K850_hg38"))))
})
