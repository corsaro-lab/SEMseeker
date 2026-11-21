# Tests for the public semseeker() dispatcher and its helpers.
# The dispatcher is covered by unit tests; the full end-to-end pipeline is
# exercised by test-e2e_semseeker.R.

# ---- .sem_detect_input_type -------------------------------------------------

test_that(".sem_detect_input_type recognises bedmethyl file paths", {
  expect_equal(SEMseeker:::.sem_detect_input_type("sample.bed"),       "bedmethyl")
  expect_equal(SEMseeker:::.sem_detect_input_type("sample.tsv"),       "bedmethyl")
  expect_equal(SEMseeker:::.sem_detect_input_type("sample.bedmethyl"), "bedmethyl")
  expect_equal(SEMseeker:::.sem_detect_input_type(c("a.bed","b.bed")), "bedmethyl")
})

test_that(".sem_detect_input_type recognises coord data frames", {
  df <- data.frame(CHR = "chr1", START = 100L, s1 = 0.5, s2 = 0.7)
  expect_equal(SEMseeker:::.sem_detect_input_type(df), "coord_df")
})

test_that(".sem_detect_input_type recognises probe-indexed matrices", {
  m <- matrix(runif(6), nrow = 3, dimnames = list(
    c("cg00000029","cg00000165","cg00000236"), c("s1","s2")))
  expect_equal(SEMseeker:::.sem_detect_input_type(m), "matrix")

  df <- data.frame(s1 = c(0.1, 0.2), s2 = c(0.3, 0.4),
                   row.names = c("cg00000029", "cg00000165"))
  expect_equal(SEMseeker:::.sem_detect_input_type(df), "matrix")
})

test_that(".sem_detect_input_type rejects unsupported extensions and classes", {
  expect_error(SEMseeker:::.sem_detect_input_type("data.csv"),  "unsupported extensions")
  expect_error(SEMseeker:::.sem_detect_input_type(list(1, 2)),  "cannot infer input_type")
})

# ---- sem_mvalue_to_beta -----------------------------------------------------

test_that("sem_mvalue_to_beta converts a matrix via 2^M / (1 + 2^M)", {
  m <- matrix(c(-2, 0, 2), ncol = 1)
  beta <- SEMseeker:::sem_mvalue_to_beta(m)
  expect_equal(as.numeric(beta), c(0.2, 0.5, 0.8), tolerance = 1e-6)
})

test_that("sem_mvalue_to_beta preserves CHR/START/END columns in coord data frames", {
  df <- data.frame(CHR = "chr1", START = 100L, END = 101L,
                   s1 = c(-2, 0, 2),
                   stringsAsFactors = FALSE)
  out <- SEMseeker:::sem_mvalue_to_beta(df, coord_cols = c("CHR","START","END"))
  expect_equal(out$CHR,   rep("chr1", 3))
  expect_equal(out$START, rep(100L, 3))
  expect_equal(out$s1,    c(0.2, 0.5, 0.8), tolerance = 1e-6)
})
