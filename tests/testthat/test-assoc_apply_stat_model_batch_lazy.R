# assoc_apply_stat_model_batch_lazy() on a lazy table of the data engine.
#
# The resume filter (area_to_remove) once crashed a pipeline because the
# engine of the time read the AREA values as column names (2026-06-05). The
# tests below pin the contract on the function itself, not on the engine:
# the areas already done are removed, the others are fitted, and the fit
# equals limma run directly on the same matrix. The expected values come from
# limma on a base-R matrix, never from the code under test.

# Build a small pivot and sample sheet shared by the tests.
# Returns a list with `pivot` (AREA + one column per sample), `values` (the
# same numbers as an areas x samples matrix) and `sample_sheet`.
batch_lazy_fixture <- function() {
  set.seed(20261120L)
  sample_count <- 10L
  area_count <- 25L
  values <- matrix(stats::rnorm(area_count * sample_count, mean = 3), nrow = area_count)
  colnames(values) <- sprintf("S%02d", seq_len(sample_count))
  rownames(values) <- c("A1BG", "B2M", "C3", sprintf("GENE_%02d", seq_len(area_count - 3L)))
  list(
    pivot = data.frame(AREA = rownames(values), values, check.names = FALSE, stringsAsFactors = FALSE),
    values = values,
    sample_sheet = data.frame(Sample_ID = colnames(values), AGE = seq(20, 65, length.out = sample_count),
                              stringsAsFactors = FALSE)
  )
}

batch_lazy_key <- data.frame(MARKER = "DELTAR", FIGURE = "HYPER", AREA = "GENE", SUBAREA = "WHOLE",
                             stringsAsFactors = FALSE)

test_that("area_to_remove drops the areas already done and keeps the others in order", {
  skip_if_not_installed("limma")
  fixture <- batch_lazy_fixture()
  area_to_remove <- c("A1BG", "B2M")
  result <- SEMseeker:::assoc_apply_stat_model_batch_lazy(
    pivot_lazy = SEMseeker:::io_register_frame(fixture$pivot), sample_sheet = fixture$sample_sheet,
    family_test = "limma_1", key = batch_lazy_key, transformation_y = "none",
    independent_variable = "AGE", area_to_remove = area_to_remove)
  expect_equal(result$AREA_OF_TEST, setdiff(rownames(fixture$values), area_to_remove))
})

test_that("the batch fit equals limma on the same matrix", {
  skip_if_not_installed("limma")
  fixture <- batch_lazy_fixture()
  result <- SEMseeker:::assoc_apply_stat_model_batch_lazy(
    pivot_lazy = SEMseeker:::io_register_frame(fixture$pivot), sample_sheet = fixture$sample_sheet,
    family_test = "limma_1", key = batch_lazy_key, transformation_y = "none",
    independent_variable = "AGE")
  design <- cbind(`(Intercept)` = 1, AGE = fixture$sample_sheet$AGE)
  expected <- limma::eBayes(limma::lmFit(fixture$values, design))
  expect_equal(result$AREA_OF_TEST, rownames(fixture$values))
  expect_equal(result$PVALUE, unname(expected$p.value[, "AGE"]), tolerance = 1e-10)
})

test_that("nothing left after the resume filter returns NULL", {
  skip_if_not_installed("limma")
  fixture <- batch_lazy_fixture()
  result <- SEMseeker:::assoc_apply_stat_model_batch_lazy(
    pivot_lazy = SEMseeker:::io_register_frame(fixture$pivot), sample_sheet = fixture$sample_sheet,
    family_test = "limma_1", key = batch_lazy_key, transformation_y = "none",
    independent_variable = "AGE", area_to_remove = rownames(fixture$values))
  expect_null(result)
})
