# enrich_term_lollipop_plot(), which through 0.99.6 was enrich_lollipop_plot()
# and could not run: performance_category was used six times and assigned
# nowhere, the ggplot object was never assigned so ggsave() saved last_plot(),
# and the descriptions were starred from a column produced only by a function
# with no callers. The drawing is asserted on through ggplot_build().

.lollipop_session <- function() {
  folder <- sem_test_folder()
  dir.create(file.path(folder, "Data"), recursive = TRUE, showWarnings = FALSE)
  SEMseeker:::core_init_env(folder, parallel_strategy = "sequential",
                            iqrTimes = 3, verbosity = 1)
  folder
}

.rules_of <- function(label) {
  formats <- as.data.frame(SEMseeker:::core_get_session_info()$key_enrichment_format)
  formats[formats$label == label, ]
}

# The shape the live pipeline writes: the enricher's OWN column names, plus the
# SS_* columns enrich_analysy_add_category() adds.
.raw_report <- function() data.frame(
  ID = c("hsa04110", "hsa04151", "hsa05200"),
  Term_Description = c("Cell cycle", "PI3K-Akt", "Cancer"),
  highest_p = c(1e-5, 1e-3, 0),
  Fold_Enrichment = c(3.2, 2.1, 1.8),
  SS_CATEGORY = c("KEGG", "KEGG", "KEGG"),
  SS_RANK = c(1, 2, 3),
  MARKER = c("DELTAS", "DELTAS", "DELTAR"),
  stringsAsFactors = FALSE)

.canonical <- function(report = .raw_report(), label = "pathfindR")
  SEMseeker:::.enrich_report_canonical(report, .rules_of(label))

.build <- function(data, ...) {
  session <- SEMseeker:::core_get_session_info()
  SEMseeker:::.enrich_lollipop_plot_build(
    data, palette = session$color_palette, alpha = 0.05, ...)
}

# ---------------------------------------------------------------------------
# the enricher's names become canonical ones
# ---------------------------------------------------------------------------

test_that("the four columns the chart needs are renamed from the rules row", {
  skip_on_cran()
  folder <- .lollipop_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  canonical <- .canonical()
  expect_true(all(c("ID", "Description", "P_Value", "Enrichment") %in%
                    colnames(canonical)))
  # the values travel with the name: highest_p is pathfindR's adjusted p-value.
  expect_equal(canonical$P_Value, c(1e-5, 1e-3, 0))
  expect_equal(canonical$Description[1], "Cell cycle")
})

test_that("another enricher's names are renamed from its own row", {
  skip_on_cran()
  folder <- .lollipop_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  # WebGestalt keeps its terms in geneSet/description/FDR/enrichmentRatio. The
  # drawing knows none of this, which is the point of renaming once.
  webgestalt <- data.frame(
    geneSet = "GO:0006915", description = "apoptosis", FDR = 0.001,
    enrichmentRatio = 2.5, SS_CATEGORY = "GO-BP", SS_RANK = 1,
    MARKER = "DELTAS", stringsAsFactors = FALSE)

  canonical <- .canonical(webgestalt, "WebGestalt")
  expect_equal(canonical$Description, "apoptosis")
  expect_equal(canonical$P_Value, 0.001)
  expect_equal(canonical$Enrichment, 2.5)
})

test_that("a declared column absent from the report is refused, and says which", {
  skip_on_cran()
  folder <- .lollipop_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  report <- .raw_report()
  report$Fold_Enrichment <- NULL
  expect_error(.canonical(report), "Fold_Enrichment")
  # the message says which canonical role the missing column was for.
  expect_error(.canonical(report), "Enrichment")
})

test_that("a report with no category is one category, not an undrawable chart", {
  skip_on_cran()
  folder <- .lollipop_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  report <- .raw_report()
  report$SS_CATEGORY <- NULL
  expect_equal(unique(.canonical(report)$SS_CATEGORY), "ALL")
})

# ---------------------------------------------------------------------------
# the drawing
# ---------------------------------------------------------------------------

test_that("the chart renders, one point per term", {
  skip_on_cran()
  folder <- .lollipop_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  built <- ggplot2::ggplot_build(
    .build(.canonical(), group_column = "MARKER", sort_column = "SS_RANK"))
  expect_equal(nrow(built$data[[1]]), 3L)
})

test_that("a p-value of zero gets a floor instead of leaving the chart", {
  skip_on_cran()
  folder <- .lollipop_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  built <- ggplot2::ggplot_build(
    .build(.canonical(), group_column = "MARKER", sort_column = "SS_RANK"))

  # -log10(0) is an infinity, which would silently drop the MOST significant
  # term. The floor is the smallest positive value in the same report, 1e-5.
  expect_true(all(is.finite(built$data[[1]]$x)))
  expect_equal(sum(built$data[[1]]$x == 5), 2L)
})

test_that("the group column gives both the colour and the shape", {
  skip_on_cran()
  folder <- .lollipop_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  plot_object <- .build(.canonical(), group_column = "MARKER",
                        sort_column = "SS_RANK")
  built <- ggplot2::ggplot_build(plot_object)

  # two markers in the fixture, so two colours and two shapes.
  expect_equal(length(unique(built$data[[1]]$colour)), 2L)
  expect_equal(length(unique(built$data[[1]]$shape)), 2L)
  # and the legend is titled with the column, so the reader knows what varies.
  expect_equal(plot_object$labels$colour, "MARKER")
  expect_equal(plot_object$labels$shape, "MARKER")
})

test_that("another column can group the points, which is why it is a parameter", {
  skip_on_cran()
  folder <- .lollipop_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  canonical <- .canonical()
  canonical$PHENOTYPE <- c(TRUE, FALSE, FALSE)
  plot_object <- .build(canonical, group_column = "PHENOTYPE",
                        sort_column = "SS_RANK")
  expect_equal(plot_object$labels$colour, "PHENOTYPE")
  expect_equal(length(unique(ggplot2::ggplot_build(plot_object)$data[[1]]$shape)),
               2L)
})

test_that("a column the report does not have is refused, with the list", {
  skip_on_cran()
  folder <- .lollipop_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  expect_error(.build(.canonical(), group_column = "performance_category",
                      sort_column = "SS_RANK"),
               "not a column of the report")
  expect_error(.build(.canonical(), group_column = "MARKER",
                      sort_column = "NOT_A_RANK"),
               "NOT_A_RANK")
})

test_that("top keeps that many terms and does not invent the missing ones", {
  skip_on_cran()
  folder <- .lollipop_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  two <- ggplot2::ggplot_build(
    .build(.canonical(), group_column = "MARKER", sort_column = "SS_RANK",
           top = 2L))
  expect_equal(nrow(two$data[[1]]), 2L)

  # asking for more terms than exist keeps the ones that do: indexing by
  # seq_len(top) would have produced NA descriptions instead.
  many <- ggplot2::ggplot_build(
    .build(.canonical(), group_column = "MARKER", sort_column = "SS_RANK",
           top = 500L))
  expect_equal(nrow(many$data[[1]]), 3L)
})

test_that("the keyword star is drawn when the column is there and skipped when not", {
  skip_on_cran()
  folder <- .lollipop_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  # by_keyword comes from a function with no callers, so the live report has no
  # such column. Assuming it assigned a zero-length vector to n rows.
  without <- .build(.canonical(), group_column = "MARKER",
                    sort_column = "SS_RANK")
  expect_false(any(grepl("\\*", levels(without$data$Description))))

  canonical <- .canonical()
  canonical$by_keyword <- c(TRUE, FALSE, FALSE)
  with_star <- .build(canonical, group_column = "MARKER",
                      sort_column = "SS_RANK")
  expect_true(any(grepl("\\*", levels(with_star$data$Description))))
  expect_equal(sum(grepl("\\*", levels(with_star$data$Description))), 1L)
})

test_that("the most significant term is drawn at the top", {
  skip_on_cran()
  folder <- .lollipop_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  plot_object <- .build(.canonical(), group_column = "MARKER",
                        sort_column = "SS_RANK")
  # ggplot2 draws the first level at the bottom, so the rank-1 term has to be
  # the LAST level for it to appear at the top.
  expect_equal(utils::tail(levels(plot_object$data$Description), 1),
               "Cell cycle")
})

test_that("the door writes nothing when no report of that enricher exists", {
  skip_on_cran()
  folder <- .lollipop_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  inference_details <- data.frame(
    independent_variable = "Sample_Group", family_test = "wilcoxon",
    transformation_y = "none", aggregation = "MEAN", scope = "INSTANCE",
    stringsAsFactors = FALSE)

  expect_equal(length(
    SEMseeker::enrich_term_lollipop_plot("pathfindR", inference_details)), 0L)
  # and an unknown enricher is refused before any file is looked for.
  expect_error(
    SEMseeker::enrich_term_lollipop_plot("not_an_enricher", inference_details),
    "unknown enricher")
})

test_that("the lollipop serves an enricher the circos has to refuse", {
  skip_on_cran()
  folder <- .lollipop_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  inference_details <- data.frame(
    independent_variable = "Sample_Group", family_test = "wilcoxon",
    transformation_y = "none", aggregation = "MEAN", scope = "INSTANCE",
    stringsAsFactors = FALSE)

  # A lollipop of enriched terms needs no genes, so requiring column_of_genes
  # in the shared lookup would have refused six enrichers a chart that draws
  # all seven. The circos refuses WebGestalt; this one does not.
  expect_error(SEMseeker::enrich_circos_plot("WebGestalt", inference_details),
               "does not declare which column")
  expect_equal(length(
    SEMseeker::enrich_term_lollipop_plot("WebGestalt", inference_details)), 0L)
})
