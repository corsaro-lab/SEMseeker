# Regression: assoc_apply_stat_model_batch_lazy's resume-cache filter must accept an
# R character vector as `area_to_remove` without polars interpreting the values
# as column references. polars 1.x is_in() parses bare R char vectors as a list
# of column expressions; the first AREA value (e.g. "A1BG") gets resolved as a
# nonexistent column and the whole $collect() fails with
#   "Column(s) not found: 'A1BG' not found"
# (ewas_cancer_stage pipeline crash, 2026-06-05 10:25). Fix wraps the vector
# via pl$lit()$implode() so it lands as a value-set literal.

test_that("is_in() with bare R char vector misinterprets values as column refs", {
  skip_on_cran()
  df <- polars::pl$DataFrame(AREA = c("A1BG", "B2M", "C3", "D4"), v = 1:4)
  expr <- !polars::pl$col("AREA")$is_in(c("A1BG", "B2M"))
  testthat::expect_error(
    df$lazy()$filter(expr)$collect(),
    "Column.*not found"
  )
})

test_that("is_in() with pl$lit()$implode() treats vector as a literal set", {
  skip_on_cran()
  df <- polars::pl$DataFrame(AREA = c("A1BG", "B2M", "C3", "D4"), v = 1:4)
  expr <- !polars::pl$col("AREA")$is_in(
    polars::pl$lit(c("A1BG", "B2M"))$implode()
  )
  out  <- df$lazy()$filter(expr)$collect()
  out_df <- as.data.frame(out)
  testthat::expect_equal(sort(out_df$AREA), c("C3", "D4"))
  testthat::expect_equal(sort(out_df$v), c(3L, 4L))
})

test_that("assoc_apply_stat_model_batch_lazy filter mirrors the patched chain", {
  skip_on_cran()
  # Same shape used inside assoc_apply_stat_model_batch_lazy: AREA column carries
  # hyphen-normalised gene identifiers and area_to_remove holds resume-cache
  # entries (callers pre-apply gsub("-","_") on both sides).
  df <- polars::pl$DataFrame(
    AREA = c("A1BG", "GENE_X", "B2M", "C3"),
    v    = 1:4
  )
  area_to_remove <- c("A1BG", "B2M")
  filtered <- df$lazy()$filter(
    !polars::pl$col("AREA")$str$replace_all("-", "_")$is_in(
      polars::pl$lit(area_to_remove)$implode()
    )
  )$collect()
  testthat::expect_equal(filtered$height, 2L)
  testthat::expect_setequal(
    as.character(filtered$select("AREA")$to_series()),
    c("GENE_X", "C3")
  )
})
