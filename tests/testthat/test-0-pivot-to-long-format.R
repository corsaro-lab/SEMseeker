# io_pivot_to_long_format(), which through 0.99.5 raised on its second line:
# it treated the lazy polars frame io_read_pivot() returns as a data.frame.
# Every test here reads a REAL parquet pivot from disk. Stubbing the reader is
# exactly what hid the defect: the chart test above it mocked this function, so
# the only thing never exercised was the thing that did not work.

.pivot_session <- function() {
  folder <- sem_test_folder()
  dir.create(file.path(folder, "Data"), recursive = TRUE, showWarnings = FALSE)
  SEMseeker:::core_init_env(folder, parallel_strategy = "sequential",
                            iqrTimes = 3, verbosity = 1)
  folder
}

# GENE blocks hold more than one position, so the pivot is identified by its
# aggregation too. Passing it is the point: the reader used to carry four of the
# six coordinates and could not name a GENE pivot at all.
.write_pivot <- function(areas, ..., aggregation = "SUM") {
  path <- SEMseeker:::io_pivot_file_name_parquet("MUTATIONS", "HYPER",
                                                "GENE", "WHOLE",
                                                aggregation = aggregation)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  polars::as_polars_df(data.frame(AREA = areas, ..., stringsAsFactors = FALSE)
                       )$write_parquet(path)
  path
}

.long <- function(sheet, ...) {
  SEMseeker:::io_pivot_to_long_format("MUTATIONS", "HYPER", "GENE", "WHOLE",
                                     "Sample_Group", sheet, ...,
                                     aggregation = "SUM")
}

.sheet <- function(...) data.frame(..., stringsAsFactors = FALSE)

test_that("a real pivot melts to one row per instance and sample", {
  skip_on_cran()
  folder <- .pivot_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  .write_pivot(c("BRCA1", "TP53"), S1 = c(1, 2), S2 = c(3, 4))
  sheet <- .sheet(Sample_ID = c("S1", "S2"),
                  Sample_Group = c("Case", "Reference"))

  long <- .long(sheet)

  # two instances x two samples, and the first instance is NOT dropped.
  expect_equal(nrow(long), 4L)
  expect_setequal(colnames(long), c("AREA", "VALUE", "SAMPLE", "phenotype"))
  expect_setequal(unique(long$AREA), c("BRCA1", "TP53"))
  expect_setequal(unique(long$SAMPLE), c("S1", "S2"))

  # the value of a cell survives the melt, in the right cell.
  expect_equal(long$VALUE[long$AREA == "TP53" & long$SAMPLE == "S2"], 4)
  expect_equal(long$VALUE[long$AREA == "BRCA1" & long$SAMPLE == "S1"], 1)
  # and the phenotype follows the sample, not the row order.
  expect_equal(unique(long$phenotype[long$SAMPLE == "S2"]), "Reference")
})

test_that("areas_selection = NULL keeps every area", {
  skip_on_cran()
  folder <- .pivot_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  .write_pivot(c("A", "B", "C"), S1 = c(1, 2, 3))
  sheet <- .sheet(Sample_ID = "S1", Sample_Group = "Case")

  # The documented default used to filter with `%in% NULL`, FALSE everywhere.
  long <- .long(sheet)
  expect_equal(nrow(long), 3L)
  expect_setequal(unique(long$AREA), c("A", "B", "C"))
})

test_that("areas_selection keeps only what it names", {
  skip_on_cran()
  folder <- .pivot_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  .write_pivot(c("A", "B", "C"), S1 = c(1, 2, 3))
  sheet <- .sheet(Sample_ID = "S1", Sample_Group = "Case")

  long <- .long(sheet, areas_selection = c("A", "C"))
  expect_setequal(unique(long$AREA), c("A", "C"))
  expect_equal(nrow(long), 2L)

  # a selection that matches nothing is nothing to draw, not an error.
  expect_null(.long(sheet, areas_selection = "NOT_AN_AREA"))
})

test_that("a numeric phenotype stays numeric", {
  skip_on_cran()
  folder <- .pivot_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  .write_pivot("A", S1 = 1, S2 = 2)
  sheet <- .sheet(Sample_ID = c("S1", "S2"), Age = c(31, 72))

  long <- SEMseeker:::io_pivot_to_long_format("MUTATIONS", "HYPER", "GENE",
                                             "WHOLE", "Age", sheet,
                                             aggregation = "SUM")
  # the type is what decides the fill scale downstream, so it is part of the
  # contract and not an implementation detail.
  expect_true(is.numeric(long$phenotype))
  expect_equal(long$phenotype[long$SAMPLE == "S2"], 72)
})

test_that("the first matching sheet row wins when a Sample_ID repeats", {
  skip_on_cran()
  folder <- .pivot_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  .write_pivot(c("A", "B"), REF_1 = c(1, 2))
  # A reference reused across groups is the normal case, not a duplicate.
  sheet <- .sheet(Sample_ID = c("REF_1", "REF_1"),
                  Sample_Group = c("Reference", "Control"))

  long <- .long(sheet)
  # one row per instance, not one per (instance, sheet row).
  expect_equal(nrow(long), 2L)
  expect_equal(unique(long$phenotype), "Reference")
})

test_that("a sample absent from the sheet keeps its values with an NA phenotype", {
  skip_on_cran()
  folder <- .pivot_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  .write_pivot("A", S1 = 1, GHOST = 9)
  sheet <- .sheet(Sample_ID = "S1", Sample_Group = "Case")

  long <- .long(sheet)
  # dropping it would remove data from a chart without the chart saying so.
  expect_equal(nrow(long), 2L)
  expect_equal(long$VALUE[long$SAMPLE == "GHOST"], 9)
  expect_true(is.na(long$phenotype[long$SAMPLE == "GHOST"]))
})

test_that("a missing pivot is nothing to draw and a bad phenotype column is a refusal", {
  skip_on_cran()
  folder <- .pivot_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  sheet <- .sheet(Sample_ID = "S1", Sample_Group = "Case")

  # no pivot written: NULL, not an error.
  expect_null(.long(sheet))

  .write_pivot("A", S1 = 1)
  # a column the sheet does not have is a malformed request. Answering with a
  # column of NA would save a chart indistinguishable from a real one.
  expect_error(
    SEMseeker:::io_pivot_to_long_format("MUTATIONS", "HYPER", "GENE", "WHOLE",
                                       "NoSuchColumn", sheet,
                                       aggregation = "SUM"),
    "not a column of the sample sheet")
})

test_that("the aggregation is part of the identity, not an ornament", {
  skip_on_cran()
  folder <- .pivot_session()
  on.exit({ SEMseeker:::core_close_env(); unlink(folder, recursive = TRUE) },
          add = TRUE)

  # Same marker, figure, area and subarea; two different reductions of the same
  # blocks. They are two pivots and must not be confused for one.
  .write_pivot("A", S1 = 10, aggregation = "SUM")
  .write_pivot("A", S1 = 99, aggregation = "MEAN")
  sheet <- .sheet(Sample_ID = "S1", Sample_Group = "Case")

  as_sum  <- .long(sheet)
  as_mean <- SEMseeker:::io_pivot_to_long_format("MUTATIONS", "HYPER", "GENE",
                                                "WHOLE", "Sample_Group", sheet,
                                                aggregation = "MEAN")
  expect_equal(as_sum$VALUE, 10)
  expect_equal(as_mean$VALUE, 99)

  # And an area whose blocks hold more than one position cannot be read without
  # it: the name of the file is not even constructible.
  expect_error(
    SEMseeker:::io_pivot_to_long_format("MUTATIONS", "HYPER", "GENE", "WHOLE",
                                       "Sample_Group", sheet),
    "aggregation is required")
})
