# The aggregations of io_pivot_build(), checked against values worked out by
# hand on a pivot small enough to read. The expected numbers are written out,
# not recomputed by code: a recomputation would share the definition under
# test (the IQR rank, the variance denominator, the sum of nothing).
#
# Sample S1 holds 1, 2, 3, 4, 10 and one missing value:
#   SUM 20; MEAN 4; MEDIAN 3; VARIANCE (9 + 4 + 1 + 0 + 36) / 4 = 12.5;
#   IQR by the "nearest" rank floor((n - 1) q + 0.5) on the sorted values:
#   q = 0.75 -> rank 3 -> fourth value 4; q = 0.25 -> rank 1 -> second value 2;
#   IQR 4 - 2 = 2.
# Sample S2 has no value at all: SUM 0, the others missing.
# Sample S3 holds 1, 5, 20, where the "nearest" rank and DuckDB's quantile_disc
# part company: q = 0.75 -> rank floor(1.5 + 0.5) = 2 -> 20; q = 0.25 -> rank
# floor(0.5 + 0.5) = 1 -> 5; IQR 15 (polars 1.16 gives 15; quantile_disc 4).

#' Open a session and write the SIGNAL position pivot of the fixture.
#'
#' @return the session list returned by core_init_env()
aggregation_fixture <- function() {
  ssEnv <- SEMseeker:::core_init_env(sem_test_folder(), parallel_strategy = "sequential", start_fresh = TRUE)
  pivot <- data.frame(CHR = c("chr1", "1", "1", "2", "2", "2"),
                      START = c(10L, 20L, 30L, 5L, 15L, 25L),
                      END = c(11L, 21L, 31L, 6L, 16L, 26L),
                      S1 = c(1, 2, 3, 4, 10, NA),
                      S2 = rep(NA_real_, 6),
                      S3 = c(1, 5, 20, NA, NA, NA))
  pivot_path <- SEMseeker:::io_pivot_file_name_parquet("SIGNAL", "BETA", "POSITION", "WHOLE")
  dir.create(dirname(pivot_path), recursive = TRUE, showWarnings = FALSE)
  nanoparquet::write_parquet(pivot, pivot_path)
  ssEnv
}

test_that("each algebraic aggregation at scope SAMPLE gives the value worked out by hand", {
  skip_if_not_installed("nanoparquet")
  aggregation_fixture()
  on.exit(SEMseeker:::core_close_env(), add = TRUE)

  written <- SEMseeker:::io_pivot_build("SIGNAL", "BETA", scope = "SAMPLE", area = "POSITION", subarea = "WHOLE",
                                        aggregations = c("SUM", "MEAN", "MEDIAN", "VARIANCE", "IQR"))
  # The artefacts are read with nanoparquet, a reader that is neither the
  # engine under test nor polars. aggregation: one of the names of `written`.
  read_one <- function(aggregation) as.data.frame(nanoparquet::read_parquet(written[[aggregation]]))

  expected_s1 <- c(SUM = 20, MEAN = 4, MEDIAN = 3, VARIANCE = 12.5, IQR = 2)
  for (aggregation in names(expected_s1)) {
    artefact <- read_one(aggregation)
    # One row, the samples first and the key last, the key naming the operator.
    expect_identical(colnames(artefact), c("S1", "S2", "S3", "AREA"), info = aggregation)
    expect_identical(artefact$AREA, aggregation, info = aggregation)
    expect_equal(artefact$S1, expected_s1[[aggregation]], info = aggregation)
  }
  expect_equal(read_one("SUM")$S2, 0)
  expect_equal(read_one("IQR")$S3, 15)
  for (aggregation in c("MEAN", "MEDIAN", "VARIANCE", "IQR"))
    expect_true(is.na(read_one(aggregation)$S2), info = aggregation)
})

test_that("the position mask at scope INSTANCE keys each row by CHR_START", {
  aggregation_fixture()
  on.exit(SEMseeker:::core_close_env(), add = TRUE)

  masked <- SEMseeker:::.io_pivot_masked_lazy("SIGNAL", "BETA", "INSTANCE", "POSITION", "WHOLE")
  got <- as.data.frame(dplyr::collect(masked$lazy))

  expect_identical(colnames(got), c("S1", "S2", "S3", "AREA"))
  expect_setequal(got$AREA, c("1_10", "1_20", "1_30", "2_5", "2_15", "2_25"))
})

test_that("composite area names explode into one trimmed row per area", {
  aggregation_fixture()
  on.exit(SEMseeker:::core_close_env(), add = TRUE)
  connection <- SEMseeker:::io_engine_connection()

  composite <- SEMseeker:::io_register_frame(data.frame(
    AREA  = c("GENEA; GENEB", "GENEC,GENED", "HLA-A/B", "GENEE;"),
    VALUE = c(1, 2, 3, 4), stringsAsFactors = FALSE))
  exploded_query <- SEMseeker:::.anno_area_explode(connection, dbplyr::sql_render(composite))
  got <- DBI::dbGetQuery(connection, paste("SELECT * FROM (", exploded_query, ") ORDER BY AREA"))

  expect_identical(colnames(got), c("AREA", "VALUE"))
  expect_identical(got$AREA, c("GENEA", "GENEB", "GENEC", "GENED", "GENEE", "HLA-A", "HLA-B"))
  expect_equal(got$VALUE, c(1, 1, 2, 2, 4, 3, 3))
})
