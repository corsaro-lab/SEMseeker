# AI-268. `filter_p_value` defaults to TRUE (assoc_analysis.R), and until AI-257
# the default path could not complete: assoc_analysis_save_results() subset on
# SIGNIFICATIVE_ADJ, a column it never creates — assoc_results_get() builds that
# one on the way out, not on the way in. Every fixture in the suite sets the flag
# to FALSE, which is why nothing ever caught it and why nothing exercised the
# path once it started working.
#
# These tests pin what the default now does, including the part that is a
# trade-off rather than a fix: the file a filtered run leaves behind is also the
# file the resume reads.

.filter_pv_session <- function(name) {
  folder <- file.path(tempdir(), name)
  unlink(folder, recursive = TRUE)
  SEMseeker:::core_init_env(result_folder = folder, start_fresh = TRUE,
                            multiple_test_adj = "BH", alpha = 0.05)
  folder
}

# One key, six instances. Only raw p-values are supplied: the adjusted column is
# what the code under test computes, and each test states its own expectation
# with an independent stats::p.adjust() rather than reading it back from the
# artefact it is meant to be checking.
.filter_pv_results <- function() {
  data.frame(
    MARKER       = "MUTATIONS",
    FIGURE       = "HYPER",
    SCOPE        = "INSTANCE",
    AREA         = "GENE",
    SUBAREA      = "WHOLE",
    AGGREGATION  = "SUM",
    AREA_OF_TEST = c("GENE-A", "GENE-B", "GENE-C", "GENE-D", "GENE-E", "GENE-F"),
    FAMILY_TEST  = "gaussian",
    TRANSFORMATION_Y = "none",
    R_MODEL      = "stats_lm",
    PVALUE       = c(0.0001, 0.0010, 0.0200, 0.2000, 0.5000, 0.9000),
    STATISTIC_PARAMETER = c(4.1, 3.6, 2.4, 1.1, 0.5, 0.1),
    stringsAsFactors = FALSE)
}

test_that("filter_p_value = TRUE writes only the rows significant at the ALL level", {
  folder <- .filter_pv_session("test_filter_pv_true")
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(folder, recursive = TRUE) }, add = TRUE)

  f <- file.path(folder, "filtered.csv")
  SEMseeker:::assoc_analysis_save_results(.filter_pv_results(), f,
                                          family_test = "gaussian",
                                          filter_p_value = TRUE)

  expect_true(file.exists(f))
  written <- utils::read.csv2(f, stringsAsFactors = FALSE)

  adj <- grep("PVALUE_ADJ_ALL_", colnames(written), value = TRUE)
  expect_length(adj, 1L)

  # Merit, not shape: which rows should survive is computed here, from the raw
  # p-values, with an independent p.adjust(). A test that only checked
  # "fewer rows than before" would pass on a filter that dropped the wrong ones.
  raw      <- .filter_pv_results()
  expected <- raw$AREA_OF_TEST[stats::p.adjust(raw$PVALUE, method = "BH") < 0.05]

  expect_setequal(written$AREA_OF_TEST, expected)
  expect_equal(nrow(written), length(expected))
  expect_true(all(written[[adj]] < 0.05))

  # And the filter has to be doing something: if every row survived, the three
  # assertions above would hold on a no-op.
  expect_lt(length(expected), nrow(raw))
})

test_that("filter_p_value = FALSE keeps the rows that did not reach significance", {
  folder <- .filter_pv_session("test_filter_pv_false")
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(folder, recursive = TRUE) }, add = TRUE)

  f <- file.path(folder, "unfiltered.csv")
  SEMseeker:::assoc_analysis_save_results(.filter_pv_results(), f,
                                          family_test = "gaussian",
                                          filter_p_value = FALSE)

  written <- utils::read.csv2(f, stringsAsFactors = FALSE)
  expect_equal(nrow(written), 6L)
  # The control that gives the previous test its meaning: the same input, the
  # same code, and the only difference is the flag.
  expect_setequal(written$AREA_OF_TEST,
                  c("GENE-A", "GENE-B", "GENE-C", "GENE-D", "GENE-E", "GENE-F"))
})

test_that("a filtered file tells the resume that the dropped instances were never tested", {
  # AI-268, the part that is a trade-off and not a fix. The inference CSV is
  # also what .assoc_resume_done() reads to decide what not to compute again.
  # Filtering it means the instances that failed to reach significance are no
  # longer on disk, so a resumed run cannot tell them apart from instances that
  # were never tested, and will test them again.
  #
  # This is pinned rather than corrected because the correction is a design
  # choice, not a bug fix: either the filter stops applying to the file the
  # resume reads, or the resume stops reading the file the filter wrote, or the
  # waste is accepted and documented. Whichever is chosen, this test is what
  # will notice it changing.
  folder <- .filter_pv_session("test_filter_pv_resume")
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(folder, recursive = TRUE) }, add = TRUE)

  f <- file.path(folder, "filtered.csv")
  SEMseeker:::assoc_analysis_save_results(.filter_pv_results(), f,
                                          family_test = "gaussian",
                                          filter_p_value = TRUE)
  written <- utils::read.csv2(f, stringsAsFactors = FALSE)

  key <- list(MARKER = "MUTATIONS", FIGURE = "HYPER", SCOPE = "INSTANCE",
              AREA = "GENE", SUBAREA = "WHOLE", AGGREGATION = "SUM")
  done <- SEMseeker:::.assoc_resume_done(written, key)

  # The resume reports as done exactly the rows that survived the filter.
  expect_setequal(done, written$AREA_OF_TEST)

  # And therefore NOT the ones the filter removed: a resumed run recomputes
  # them. Named explicitly, so that a change of behaviour reads as a decision.
  expect_false("GENE-F" %in% done)
  expect_lt(length(done), 6L)
})
