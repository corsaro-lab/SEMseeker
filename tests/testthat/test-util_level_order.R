# util_level_order(): the order of the levels of an ordinal independent variable.
#
# The order is not presentation. It decides which comparisons are consecutive,
# the sign of a difference between two of them, and levels()[1], which is the
# reference category of every regression fitted on the variable. These assert
# the order itself where it is cheap, and its arrival in the data where it
# matters.

.stage_labels <- c("riferimento", "in situ", "I", "II", "III", "IV")
.stage_order  <- "riferimento + in situ + I + II + III + IV"

# ---- the fallback -----------------------------------------------------------

test_that("util_level_order: with no declaration the levels are mixedsorted", {
  # The property plain sorting got wrong: a number written as text reads as a
  # number, so 10 follows 2 instead of preceding it.
  got <- SEMseeker:::util_level_order(c("10", "2", "1", "3"))
  expect_equal(got, c("1", "2", "3", "10"))
  expect_false(identical(got, sort(c("10", "2", "1", "3"))))
})

test_that("util_level_order: NULL, NA and an empty string all mean no declaration", {
  v <- c("b", "a", "c")
  for (absent in list(NULL, NA, NA_character_, "", "   ")) {
    expect_equal(SEMseeker:::util_level_order(v, absent), c("a", "b", "c"),
                 info = paste("absent =", paste(format(absent), collapse = "")))
  }
})

test_that("util_level_order: the fallback does not rescue word labels", {
  # Stated as a test because it is the reason the declaration exists rather than
  # a better sort. mixedsort puts "in situ" between III and IV: a plausible
  # sequence that answers a different question from the one being asked.
  got <- SEMseeker:::util_level_order(.stage_labels)
  expect_false(identical(got, .stage_labels))
  expect_equal(which(got == "riferimento"), length(got))  # last, not first
})

# ---- the declaration --------------------------------------------------------

test_that("util_level_order: a declared order is the order", {
  expect_equal(SEMseeker:::util_level_order(.stage_labels, .stage_order), .stage_labels)
  # and it is the order whatever order the values arrive in
  expect_equal(SEMseeker:::util_level_order(rev(.stage_labels), .stage_order), .stage_labels)
})

test_that("util_level_order: a level the declaration does not place is refused", {
  expect_error(
    SEMseeker:::util_level_order(c(.stage_labels, "IVb"), .stage_order),
    "does not place every level present in the data")
  # the message names what was not placed, so the request can be corrected
  # without reading the data again
  expect_error(
    SEMseeker:::util_level_order(c(.stage_labels, "IVb"), .stage_order),
    "IVb")
})

test_that("util_level_order: a declared level absent from the subset is dropped, not refused", {
  # A samples filter legitimately removes a level, and the request should not
  # have to be rewritten per subset.
  observed <- c("riferimento", "I", "II")
  got <- SEMseeker:::util_level_order(observed, .stage_order)
  expect_equal(got, c("riferimento", "I", "II"))
})

# ---- the order reaches the data ---------------------------------------------

test_that("io_data_preparation: a grouping family gets the declared levels", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df <- data.frame(
    STAGE  = rep(.stage_labels, each = 3),
    BURDEN = as.numeric(seq_len(18)),
    stringsAsFactors = FALSE)

  out <- SEMseeker:::io_data_preparation(
    family_test = "kruskal.test", transformation_y = "none",
    tempDataFrame = df, independent_variable = "STAGE",
    g_start = 2L, g_end = 2L, covariates = character(0),
    key = list(AREA = "GENE", SUBAREA = "TSS200", MARKER = "MUTATIONS", FIGURE = "K850"),
    independent_variable_order = .stage_order)

  iv <- out$tempDataFrame[, "STAGE"]
  expect_true(is.factor(iv))
  expect_equal(levels(iv), .stage_labels)
  # levels()[1] is the reference category of any regression on this variable
  expect_equal(levels(iv)[1], "riferimento")
})

test_that("io_data_preparation: a regression family codes the labels by their declared position", {
  # Without a declared order the whole frame is coerced with as.numeric, so a
  # labelled variable became NA here and the caller had to add a second column
  # holding the same thing as a number. This is that column becoming unnecessary.
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df <- data.frame(
    STAGE  = rep(.stage_labels, each = 3),
    BURDEN = as.numeric(seq_len(18)),
    stringsAsFactors = FALSE)

  with_order <- SEMseeker:::io_data_preparation(
    family_test = "kendall", transformation_y = "none",
    tempDataFrame = df, independent_variable = "STAGE",
    g_start = 2L, g_end = 2L, covariates = character(0),
    key = list(AREA = "GENE", SUBAREA = "TSS200", MARKER = "MUTATIONS", FIGURE = "K850"),
    independent_variable_order = .stage_order)

  iv <- with_order$tempDataFrame[, "STAGE"]
  expect_true(is.numeric(iv))
  expect_equal(sort(unique(iv)), 1:6)
  # the code is the declared position, so the reference level is 1 and not
  # whatever the alphabet would have put there
  expect_equal(unique(iv[df$STAGE == "riferimento"]), 1)
  expect_equal(unique(iv[df$STAGE == "IV"]), 6)

  without_order <- SEMseeker:::io_data_preparation(
    family_test = "kendall", transformation_y = "none",
    tempDataFrame = df, independent_variable = "STAGE",
    g_start = 2L, g_end = 2L, covariates = character(0),
    key = list(AREA = "GENE", SUBAREA = "TSS200", MARKER = "MUTATIONS", FIGURE = "K850"))

  # Unchanged: with no declaration the labels still do not become numbers.
  expect_true(all(is.na(without_order$tempDataFrame[, "STAGE"])))
})
