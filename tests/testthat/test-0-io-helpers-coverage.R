## Coverage for previously-untested IO / PIVOT helpers.
##
## Covered:
##   io_guess_decimal_separator()   comma-vs-period heuristic on first 5 lines
##   .io_find_coord_cols()          CHR/START/END column detection + aliases
##   .io_make_probe_id()            synthetic "{CHR}_{START}" probe id
##   io_inference_file_name()       inference-result path assembly (session)
##   io_save_latex_table()          data.frame -> .tex via xtable
##

# ---------------------------------------------------------------------------
# io_guess_decimal_separator
# ---------------------------------------------------------------------------

test_that("io_guess_decimal_separator picks the majority separator", {
  f_comma  <- tempfile(fileext = ".csv")
  f_period <- tempfile(fileext = ".csv")
  on.exit(unlink(c(f_comma, f_period)), add = TRUE)
  writeLines(c("1,5", "2,3", "4,1"), f_comma)
  writeLines(c("1.5", "2.3", "4.1"), f_period)

  expect_equal(SEMseeker:::io_guess_decimal_separator(f_comma),  ",")
  expect_equal(SEMseeker:::io_guess_decimal_separator(f_period), ".")
})

# ---------------------------------------------------------------------------
# .io_find_coord_cols
# ---------------------------------------------------------------------------

test_that(".io_find_coord_cols detects coordinate columns and aliases", {
  d1 <- data.frame(CHR = "chr1", START = 100L, s1 = 0.5)
  expect_equal(SEMseeker:::.io_find_coord_cols(d1),
               list(chr = "CHR", start = "START", end = NA_character_))

  d2 <- data.frame(CHR = "chr1", START = 100L, END = 101L, s1 = 0.5)
  expect_equal(SEMseeker:::.io_find_coord_cols(d2)$end, "END")

  d3 <- data.frame(seqnames = "chr1", pos = 100L, s1 = 0.5)  # bedmethyl-style aliases
  cols <- SEMseeker:::.io_find_coord_cols(d3)
  expect_equal(cols$chr,   "seqnames")
  expect_equal(cols$start, "pos")

  d_none <- data.frame(s1 = 0.5, s2 = 0.3)
  expect_null(SEMseeker:::.io_find_coord_cols(d_none))
})

# ---------------------------------------------------------------------------
# .io_make_probe_id
# ---------------------------------------------------------------------------

test_that(".io_make_probe_id strips the chr prefix and joins with START", {
  expect_equal(
    SEMseeker:::.io_make_probe_id(c("chr1", "chrX"), c(10000L, 543200L)),
    c("1_10000", "X_543200"))
})

# ---------------------------------------------------------------------------
# io_inference_file_name
# ---------------------------------------------------------------------------

test_that("io_inference_file_name assembles a marker path", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  detail <- list(
    covariates            = "",
    covariates_dummy      = "",
    covariates_pca        = FALSE,
    family_test           = "gaussian",
    transformation_y      = "",
    independent_variable  = "AGE",
    transformation_x      = "",
    samples_sql_condition = NULL
  )

  path <- SEMseeker:::io_inference_file_name(detail, marker = "MUTATIONS",
                                             folder = tf, file_extension = "csv")

  expect_type(path, "character")
  expect_length(path, 1L)
  # DEPTH left the path with the concept it named.
  expect_false(grepl("DEPTH", path, ignore.case = TRUE))
  expect_match(path, "MUTATIONS", ignore.case = TRUE)
  expect_match(path, "\\.csv$", ignore.case = TRUE)
})

# ---------------------------------------------------------------------------
# io_save_latex_table
# ---------------------------------------------------------------------------

test_that("io_save_latex_table writes a .tex table with the caption", {
  skip_if_not_installed("xtable")
  csv <- tempfile(fileext = ".csv")
  tex <- sub("\\.csv$", ".tex", csv)
  on.exit(unlink(c(csv, tex)), add = TRUE)

  df <- data.frame(GENE_NAME = c("TP53", "BRCA1"), P_VALUE = c("0.01", "0.20"))
  SEMseeker:::io_save_latex_table(df, csv, caption = "My caption")

  expect_true(file.exists(tex))
  contents <- paste(readLines(tex), collapse = "\n")
  expect_match(contents, "\\\\begin\\{table\\}")
  expect_match(contents, "My caption")
})
