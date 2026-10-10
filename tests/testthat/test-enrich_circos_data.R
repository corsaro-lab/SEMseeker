# enrich_circos_plot() and the pure steps behind it. Through 0.99.6 the circos
# drawing existed twice: enrich_pathfindR_circlize() built its data and drew,
# and plot_area_plot_circlize() held the same thirteen drawing lines with its
# link frame `results` never assigned and never a parameter, so it raised on
# the line that chose the link colours. One drawing now, one door, and the
# enricher is a parameter.

.circos_session <- function() {
  folder <- sem_test_folder()
  dir.create(file.path(folder, "Data"), recursive = TRUE, showWarnings = FALSE)
  SEMseeker:::core_init_env(folder, parallel_strategy = "sequential",
                            iqrTimes = 3, verbosity = 1)
  folder
}

.pathfindr_rules <- function() {
  formats <- as.data.frame(SEMseeker:::core_get_session_info()$key_enrichment_format)
  formats[formats$label == "pathfindR", ]
}

.pathfindr_report <- function() data.frame(
  ID = c("hsa04110", "hsa04151"),
  Term_Description = c("Cell cycle", "PI3K-Akt"),
  Up_regulated = c("TP53, BRCA1", "AKT1"),
  Down_regulated = c("CDK2", ""),
  stringsAsFactors = FALSE)

.coords <- function() data.frame(
  GENE = c("TP53", "BRCA1", "CDK2", "AKT1"),
  CHR = c("chr17", "chr17", "chr12", "chr14"),
  START = c(7571720, 41196312, 55966769, 105235686),
  END = c(7590868, 41277500, 55972789, 105262088),
  stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# the rules table: the enricher is dispatched, and refused by name
# ---------------------------------------------------------------------------

test_that("the format table says where each enricher keeps its genes", {
  skip_on_cran()
  folder <- .circos_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  formats <- as.data.frame(SEMseeker:::core_get_session_info()$key_enrichment_format)
  expect_true("column_of_genes" %in% colnames(formats))
  # pathfindR splits its genes across two columns, which is why the field is
  # a list and not a single name.
  expect_equal(formats$column_of_genes[formats$label == "pathfindR"],
               "Up_regulated+Down_regulated")
  # the other six do not declare it yet, and that is a fact about this package
  # rather than about those enrichers.
  expect_equal(sum(nzchar(formats$column_of_genes)), 1L)
})

test_that("an unknown enricher and an undeclared one are both refused, by name", {
  skip_on_cran()
  folder <- .circos_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)
  session <- SEMseeker:::core_get_session_info()

  expect_error(SEMseeker:::.enrich_rules_for("not_an_enricher", session),
               "unknown enricher")
  # the message lists what it does know, so the caller does not have to guess.
  expect_error(SEMseeker:::.enrich_rules_for("not_an_enricher", session),
               "pathfindR")

  # An enricher with no gene column is a refusal and not an empty circle: an
  # empty circos looks exactly like a study whose enrichment found nothing.
  # The requirement lives in its own function and not in the lookup, because a
  # lollipop of enriched terms needs no genes and must serve all seven.
  web_rules <- SEMseeker:::.enrich_rules_for("WebGestalt", session)
  expect_equal(as.character(web_rules$label), "WebGestalt")
  expect_error(SEMseeker:::.enrich_gene_columns_of(web_rules, session),
               "does not declare which column")
  expect_error(SEMseeker:::.enrich_gene_columns_of(web_rules, session),
               "pathfindR")

  # and the enricher that does declare them gets both of its columns back.
  expect_setequal(
    SEMseeker:::.enrich_gene_columns_of(
      SEMseeker:::.enrich_rules_for("pathfindR", session), session),
    c("Up_regulated", "Down_regulated"))

  # and the lookup does not care about the case of the name.
  expect_equal(as.character(SEMseeker:::.enrich_rules_for("PATHFINDR",
                                                         session)$label),
               "pathfindR")
})

# ---------------------------------------------------------------------------
# the pairs
# ---------------------------------------------------------------------------

test_that("the gene columns of a report become one row per gene and term", {
  skip_on_cran()
  folder <- .circos_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  pairs <- SEMseeker:::.enrich_circos_pairs_build(.pathfindr_report(),
                                                 .pathfindr_rules())
  # both declared columns are read: TP53 and BRCA1 are up, CDK2 is down.
  expect_setequal(pairs$GENE[pairs$TERM == "hsa04110"],
                  c("TP53", "BRCA1", "CDK2"))
  # an empty cell contributes nothing rather than an empty gene name.
  expect_equal(pairs$GENE[pairs$TERM == "hsa04151"], "AKT1")
  expect_equal(nrow(pairs), 4L)
  # the description travels with the term, for the label on the chart.
  expect_equal(unique(pairs$TERM_LABEL[pairs$TERM == "hsa04151"]), "PI3K-Akt")
})

test_that("a term selection keeps only what it names, and nothing is NULL", {
  skip_on_cran()
  folder <- .circos_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  one <- SEMseeker:::.enrich_circos_pairs_build(
    .pathfindr_report(), .pathfindr_rules(), terms_selection = "hsa04151")
  expect_equal(nrow(one), 1L)
  expect_equal(one$TERM, "hsa04151")

  expect_null(SEMseeker:::.enrich_circos_pairs_build(
    .pathfindr_report(), .pathfindr_rules(), terms_selection = "nothing"))
})

test_that("a declared column absent from the report is a refusal, named", {
  skip_on_cran()
  folder <- .circos_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  report <- .pathfindr_report()
  report$Down_regulated <- NULL
  expect_error(
    SEMseeker:::.enrich_circos_pairs_build(report, .pathfindr_rules()),
    "Down_regulated")
})

# ---------------------------------------------------------------------------
# the layout
# ---------------------------------------------------------------------------

test_that("each term gets an equal slice and the last one ends at the sector", {
  skip_on_cran()
  folder <- .circos_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  pairs <- SEMseeker:::.enrich_circos_pairs_build(.pathfindr_report(),
                                                 .pathfindr_rules())
  layout <- SEMseeker:::.enrich_circos_layout_build(pairs, .coords(),
                                                   "term", 2e8)

  slices <- unique(layout$links_to[, c("START", "END")])
  expect_equal(nrow(slices), 2L)
  expect_equal(length(unique(slices$END - slices$START)), 1L)
  # starting at 1 and adding the width twice would end one past the sector,
  # which circlize answers with a warning rather than an error.
  expect_equal(max(layout$links_to$END), 2e8)
  expect_equal(min(layout$links_to$START), 1)
})

test_that("the terms are laid out in the order of the report, not of a merge", {
  skip_on_cran()
  folder <- .circos_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  pairs <- SEMseeker:::.enrich_circos_pairs_build(.pathfindr_report(),
                                                 .pathfindr_rules())
  layout <- SEMseeker:::.enrich_circos_layout_build(pairs, .coords(),
                                                   "term", 2e8)

  # hsa04110 is first in the report, so it holds the first slice. A report is
  # normally sorted by significance, so this is what makes the picture
  # reproducible instead of depending on merge order.
  first_slice <- layout$labels[layout$labels$CHR == "term" &
                               layout$labels$START == 1, "LABEL"]
  expect_equal(first_slice, "Cell cycle")
})

test_that("a link joins a gene to its own term, through two merges", {
  skip_on_cran()
  folder <- .circos_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  pairs <- SEMseeker:::.enrich_circos_pairs_build(.pathfindr_report(),
                                                 .pathfindr_rules())
  layout <- SEMseeker:::.enrich_circos_layout_build(pairs, .coords(),
                                                   "term", 2e8)

  # AKT1 is the only gene of hsa04151, which holds the second slice.
  akt_row <- which(layout$links_from$CHR == "chr14")
  expect_equal(length(akt_row), 1L)
  expect_equal(layout$links_to$START[akt_row], 1e8)

  # and the three genes of hsa04110 all land on the first slice.
  other_rows <- setdiff(seq_len(nrow(layout$links_from)), akt_row)
  expect_equal(unique(layout$links_to$START[other_rows]), 1)
})

test_that("a gene with no coordinates loses its link and is counted", {
  skip_on_cran()
  folder <- .circos_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  pairs <- SEMseeker:::.enrich_circos_pairs_build(.pathfindr_report(),
                                                 .pathfindr_rules())
  pairs <- rbind(pairs, data.frame(TERM = "hsa04151", TERM_LABEL = "PI3K-Akt",
                                   GENE = "NOT_ON_THE_ARRAY",
                                   stringsAsFactors = FALSE))

  layout <- SEMseeker:::.enrich_circos_layout_build(pairs, .coords(),
                                                   "term", 2e8)
  # a link needs two ends, so the gene is dropped rather than drawn at zero.
  expect_equal(nrow(layout$links_from), 4L)
  expect_false("NOT_ON_THE_ARRAY" %in% layout$labels$LABEL)

  # and no coordinates at all is nothing to draw, not an error.
  expect_null(SEMseeker:::.enrich_circos_layout_build(
    pairs, .coords()[0, ], "term", 2e8))
})

test_that("the labels carry both the genes and the terms", {
  skip_on_cran()
  folder <- .circos_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  pairs <- SEMseeker:::.enrich_circos_pairs_build(.pathfindr_report(),
                                                 .pathfindr_rules())
  layout <- SEMseeker:::.enrich_circos_layout_build(pairs, .coords(),
                                                   "term", 2e8)

  expect_setequal(layout$labels$LABEL[layout$labels$CHR != "term"],
                  c("TP53", "BRCA1", "CDK2", "AKT1"))
  expect_setequal(layout$labels$LABEL[layout$labels$CHR == "term"],
                  c("Cell cycle", "PI3K-Akt"))
})

# ---------------------------------------------------------------------------
# the drawing
# ---------------------------------------------------------------------------

test_that("the circos is drawn, and the genome build is the one asked for", {
  skip_on_cran()
  skip_if_not_installed("circlize")
  folder <- .circos_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  pairs <- SEMseeker:::.enrich_circos_pairs_build(.pathfindr_report(),
                                                 .pathfindr_rules())
  layout <- SEMseeker:::.enrich_circos_layout_build(pairs, .coords(),
                                                   "term", 2e8)

  for (build in c("hg19", "hg38")) {
    path <- tempfile(fileext = ".png")
    # The two charts this replaces read the cytoband with no species, which is
    # hg19: a study on hg38 was drawn against the wrong ideogram, with every
    # gene at the wrong cytoband and nothing saying so.
    expect_silent(SEMseeker:::.plot_circos_links_build(
      layout$labels, layout$links_from, layout$links_to,
      genome_build = build, plot_path = path, plot_format = "png",
      resolution_ppi = 150, extra_sector = "term"))
    expect_true(file.exists(path))
    expect_gt(file.size(path), 10000)
  }
})

test_that("the drawing works with no links at all", {
  skip_on_cran()
  skip_if_not_installed("circlize")
  folder <- .circos_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  labels <- data.frame(CHR = "chr1", START = 1e6, END = 1.1e6, LABEL = "BRCA1",
                       stringsAsFactors = FALSE)
  path <- tempfile(fileext = ".png")
  # an ideogram with labels is a usable chart, not a degenerate one.
  SEMseeker:::.plot_circos_links_build(labels, NULL, NULL, "hg19", path,
                                       "png", 150)
  expect_true(file.exists(path))
  expect_gt(file.size(path), 10000)
})

test_that("the drawing refuses a malformed request before opening a device", {
  skip_on_cran()
  skip_if_not_installed("circlize")

  labels <- data.frame(CHR = "chr1", START = 1, END = 2, LABEL = "G",
                       stringsAsFactors = FALSE)
  links <- data.frame(CHR = "chr1", START = 1, END = 2,
                      stringsAsFactors = FALSE)
  path <- tempfile(fileext = ".png")

  # a missing label column
  expect_error(SEMseeker:::.plot_circos_links_build(
    labels[, c("CHR", "START", "END")], NULL, NULL, "hg19", path, "png", 150),
    "LABEL")
  # one side of the links only
  expect_error(SEMseeker:::.plot_circos_links_build(
    labels, links, NULL, "hg19", path, "png", 150), "pass both or neither")
  # and a row count that cannot pair
  expect_error(SEMseeker:::.plot_circos_links_build(
    labels, rbind(links, links), links, "hg19", path, "png", 150),
    "same number of rows")

  # nothing was drawn, so no half-written file is left behind looking valid.
  expect_false(file.exists(path))
})

test_that("the door writes nothing when no report of that enricher exists", {
  skip_on_cran()
  folder <- .circos_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  inference_details <- data.frame(
    independent_variable = "Sample_Group", family_test = "wilcoxon",
    transformation_y = "none", aggregation = "MEAN", scope = "INSTANCE",
    stringsAsFactors = FALSE)

  written <- SEMseeker::enrich_circos_plot("pathfindR", inference_details)
  expect_equal(length(written), 0L)
})
