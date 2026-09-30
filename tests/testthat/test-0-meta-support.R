# Support functions shared by the meta_ family: the studies table, the space of
# region classes the studies have in common, and the selection inside it.
#
# These are the three places where a wrong answer is silent rather than loud: a
# study that loses its label merges with another, a class that is not shared
# disappears instead of being reported, a selection that names something nobody
# measured returns fewer rows rather than an error. Each case below fixes one of
# those, and each refusal is asserted on the reason it gives, not just on the
# fact that it stopped.

# A folder is a study when a session was written to it and the inference results
# are where the readers look. Built by hand rather than by running a pipeline:
# these functions only look at the shape.
.fake_study <- function(parent, name, with_session = TRUE, with_inference = TRUE,
                        inference_empty = FALSE) {
  folder <- file.path(parent, name)
  dir.create(file.path(folder, "Log"), recursive = TRUE, showWarnings = FALSE)
  if (with_session) saveRDS(list(marker = name), file.path(folder, "Log", "session_info.rds"))
  if (with_inference) {
    dir.create(file.path(folder, "Inference"), recursive = TRUE, showWarnings = FALSE)
    if (!inference_empty)
      writeLines("AREA;SUBAREA", file.path(folder, "Inference", "INFERENCE_RESULT.csv"))
  }
  folder
}

.keys <- function(...) {
  rows <- list(...)
  do.call(rbind, lapply(rows, function(r)
    data.frame(MARKER = r[1], FIGURE = r[2], AREA = r[3], SUBAREA = r[4],
               stringsAsFactors = FALSE)))
}

# --- the studies table -------------------------------------------------------

test_that("an unnamed path is labelled by the last component of its path", {
  root <- sem_test_folder()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  a <- .fake_study(root, "BWS_2024")
  b <- .fake_study(root, "BWS_2025")

  out <- SEMseeker:::meta_studies_normalise(c(a, b))

  expect_identical(colnames(out), c("STUDY", "STUDY_FOLDER"))
  expect_identical(out$STUDY, c("BWS_2024", "BWS_2025"))
  expect_identical(out$STUDY_FOLDER, c(a, b))
})

test_that("a trailing separator does not change the label", {
  root <- sem_test_folder()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  a <- .fake_study(root, "BWS_2024")
  b <- .fake_study(root, "BWS_2025")

  out <- SEMseeker:::meta_studies_normalise(c(paste0(a, "/"), b))
  expect_identical(out$STUDY, c("BWS_2024", "BWS_2025"))
})

test_that("a name overrides the path component, which is the way out when the folder is called results", {
  root <- sem_test_folder()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  a <- .fake_study(file.path(root, "first"), "results")
  b <- .fake_study(file.path(root, "second"), "results")

  # Unnamed, both would be called "results" and merge into one row.
  expect_error(SEMseeker:::meta_studies_normalise(c(a, b)), "same label")

  out <- SEMseeker:::meta_studies_normalise(c(BWS_2024 = a, BWS_2025 = b))
  expect_identical(out$STUDY, c("BWS_2024", "BWS_2025"))
})

test_that("the data.frame form the overlap functions take is accepted, and RESULT_FOLDER is read as the study folder", {
  root <- sem_test_folder()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  a <- .fake_study(root, "one")
  b <- .fake_study(root, "two")

  out <- SEMseeker:::meta_studies_normalise(
    data.frame(STUDY = c("A", "B"), RESULT_FOLDER = c(a, b), stringsAsFactors = FALSE))

  expect_identical(out$STUDY, c("A", "B"))
  expect_identical(out$STUDY_FOLDER, c(a, b))
})

test_that("a folder that exists but is not a study is refused, and the reason names what is missing", {
  root <- sem_test_folder()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  good <- .fake_study(root, "good")

  no_session <- .fake_study(root, "no_session", with_session = FALSE)
  expect_error(SEMseeker:::meta_studies_normalise(c(good, no_session)),
               "session_info\\.rds")

  no_inference <- .fake_study(root, "no_inference", with_inference = FALSE)
  expect_error(SEMseeker:::meta_studies_normalise(c(good, no_inference)),
               "no Inference folder")

  empty_inference <- .fake_study(root, "empty_inference", inference_empty = TRUE)
  expect_error(SEMseeker:::meta_studies_normalise(c(good, empty_inference)),
               "Inference folder is empty")

  expect_error(SEMseeker:::meta_studies_normalise(c(good, file.path(root, "nowhere"))),
               "does not exist")
})

test_that("an empty studies argument is refused rather than producing an empty table", {
  expect_error(SEMseeker:::meta_studies_normalise(character(0)), "empty")
  expect_error(SEMseeker:::meta_studies_normalise(data.frame(NOPE = 1)), "STUDY")
  expect_error(SEMseeker:::meta_studies_normalise(list("a", "b")), "character vector")
})

# --- the space of region classes ---------------------------------------------

test_that("the space is the strict intersection, and what strictness costs is reported", {
  shared <- c("MUTATIONS", "HYPO", "GENE", "WHOLE")
  key_sets <- list(
    A = .keys(shared, c("MUTATIONS", "HYPO", "ISLAND", "WHOLE")),
    B = .keys(shared, c("DELTAS", "HYPO", "GENE", "WHOLE")),
    C = .keys(shared)
  )

  out <- SEMseeker:::meta_key_space(key_sets)

  # Only the class all three declared.
  expect_equal(nrow(out$space), 1L)
  expect_identical(as.character(unlist(out$space[1, ])), shared)

  # The two that are not shared are named, with who lacked them.
  expect_equal(nrow(out$excluded), 2L)
  island <- out$excluded[out$excluded$AREA == "ISLAND", ]
  expect_equal(nrow(island), 1L)
  expect_setequal(strsplit(island$STUDIES_MISSING, " ")[[1]], c("B", "C"))
  deltas <- out$excluded[out$excluded$MARKER == "DELTAS", ]
  expect_setequal(strsplit(deltas$STUDIES_MISSING, " ")[[1]], c("A", "C"))
})

test_that("an empty intersection stops, and one study is not a meta analysis", {
  key_sets <- list(
    A = .keys(c("MUTATIONS", "HYPO", "GENE", "WHOLE")),
    B = .keys(c("DELTAS", "HYPER", "ISLAND", "WHOLE"))
  )
  expect_error(SEMseeker:::meta_key_space(key_sets), "share no region class")

  expect_error(SEMseeker:::meta_key_space(list(A = .keys(c("M", "F", "A", "S")))),
               "at least two studies")
  expect_error(SEMseeker:::meta_key_space(list(.keys(c("M", "F", "A", "S")),
                                               .keys(c("M", "F", "A", "S")))),
               "named by study")
})

test_that("classes shared by every study all survive, so the intersection is not merely small", {
  same <- list(
    A = .keys(c("MUTATIONS", "HYPO", "GENE", "WHOLE"), c("MUTATIONS", "HYPER", "GENE", "TSS200")),
    B = .keys(c("MUTATIONS", "HYPER", "GENE", "TSS200"), c("MUTATIONS", "HYPO", "GENE", "WHOLE"))
  )
  out <- SEMseeker:::meta_key_space(same)
  expect_equal(nrow(out$space), 2L)
  expect_equal(nrow(out$excluded), 0L)
})

# --- the selection inside the space ------------------------------------------

test_that("a value no study measured is refused and named, and the message says what is available", {
  space <- .keys(c("MUTATIONS", "HYPO", "GENE", "WHOLE"),
                 c("DELTAS", "HYPO", "ISLAND", "WHOLE"))

  err <- tryCatch(SEMseeker:::meta_keys_select(space, markers = "MUTAZIONI"),
                  error = function(e) conditionMessage(e))
  expect_match(err, "MUTAZIONI")
  expect_match(err, "MUTATIONS")
  expect_match(err, "DELTAS")

  expect_error(SEMseeker:::meta_keys_select(space, areas = "PROBE"), "no study measured area PROBE")
})

test_that("the selection is checked per coordinate and applied as a product intersected with the space", {
  space <- .keys(c("MUTATIONS", "HYPO", "GENE",   "WHOLE"),
                 c("DELTAS",    "HYPO", "ISLAND", "WHOLE"))

  # Both markers and both areas exist, so nothing is refused, but two of the
  # four combinations were never measured: they are dropped and counted, not
  # turned into an error.
  out <- SEMseeker:::meta_keys_select(space,
                                      markers = c("MUTATIONS", "DELTAS"),
                                      areas = c("GENE", "ISLAND"))
  expect_equal(nrow(out$keys), 2L)
  expect_equal(out$requested, 4L)
  expect_equal(out$dropped, 2L)
})

test_that("no selection takes the whole space and drops nothing", {
  space <- .keys(c("MUTATIONS", "HYPO", "GENE", "WHOLE"),
                 c("DELTAS", "HYPO", "ISLAND", "WHOLE"))
  out <- SEMseeker:::meta_keys_select(space)
  expect_equal(nrow(out$keys), 2L)
  expect_equal(out$dropped, 0L)
})

test_that("a selection whose every value exists but whose every combination is absent stops", {
  space <- .keys(c("MUTATIONS", "HYPO", "GENE",   "WHOLE"),
                 c("DELTAS",    "HYPO", "ISLAND", "WHOLE"))
  expect_error(SEMseeker:::meta_keys_select(space, markers = "MUTATIONS", areas = "ISLAND"),
               "none of them")
})

# --- the request's files, study by study -------------------------------------
#
# The file name is derived by the producer here rather than written by hand,
# because the subject of these tests is the check and not the naming. What keeps
# them honest is the second end: a file created for one marker must NOT satisfy
# the check for another, which only holds if the name really varies with the
# marker.

.inference_detail <- function() {
  data.frame(independent_variable = "Sample_Group",
             family_test = "wilcoxon",
             transformation_y = "none",
             transformation_x = "none",
             covariates = "",
             covariates_dummy = "",
             covariates_pca = FALSE,
             samples_sql_condition = "",
             aggregation = "SUM",
             scope = "INSTANCE",
             stringsAsFactors = FALSE)
}

test_that("a study holding the request's file passes, and one that does not is named with its marker", {
  root <- sem_test_folder()
  session <- sem_test_folder()
  on.exit({
    SEMseeker:::core_close_env()
    unlink(c(root, session), recursive = TRUE)
  }, add = TRUE)
  SEMseeker:::core_init_env(session, parallel_strategy = "sequential")

  detail <- .inference_detail()
  keys <- .keys(c("MUTATIONS", "HYPO", "GENE", "WHOLE"))

  a <- .fake_study(root, "with_it")
  b <- .fake_study(root, "without_it")
  studies <- SEMseeker:::meta_studies_normalise(c(a, b))

  # Only the first study gets the file the request names.
  wanted <- SEMseeker:::io_inference_file_name(detail, "MUTATIONS",
                                               file.path(a, "Inference"),
                                               skip_dir_create = TRUE)
  dir.create(dirname(wanted), recursive = TRUE, showWarnings = FALSE)
  writeLines("AREA;SUBAREA", wanted)

  err <- tryCatch(SEMseeker:::meta_inference_files_check(studies, keys, detail),
                  error = function(e) conditionMessage(e))
  expect_match(err, "without_it")
  expect_match(err, "MUTATIONS")
  expect_false(grepl("with_it /", err))

  # Give the second one its file too and the check passes.
  wanted_b <- SEMseeker:::io_inference_file_name(detail, "MUTATIONS",
                                                 file.path(b, "Inference"),
                                                 skip_dir_create = TRUE)
  dir.create(dirname(wanted_b), recursive = TRUE, showWarnings = FALSE)
  writeLines("AREA;SUBAREA", wanted_b)
  expect_silent(SEMseeker:::meta_inference_files_check(studies, keys, detail))
})

test_that("a file present for one marker does not satisfy the check for another", {
  root <- sem_test_folder()
  session <- sem_test_folder()
  on.exit({
    SEMseeker:::core_close_env()
    unlink(c(root, session), recursive = TRUE)
  }, add = TRUE)
  SEMseeker:::core_init_env(session, parallel_strategy = "sequential")

  detail <- .inference_detail()
  a <- .fake_study(root, "one_marker_only")
  studies <- SEMseeker:::meta_studies_normalise(a)

  wanted <- SEMseeker:::io_inference_file_name(detail, "MUTATIONS",
                                               file.path(a, "Inference"),
                                               skip_dir_create = TRUE)
  dir.create(dirname(wanted), recursive = TRUE, showWarnings = FALSE)
  writeLines("AREA;SUBAREA", wanted)

  expect_silent(SEMseeker:::meta_inference_files_check(
    studies, .keys(c("MUTATIONS", "HYPO", "GENE", "WHOLE")), detail))
  expect_error(SEMseeker:::meta_inference_files_check(
    studies, .keys(c("DELTAS", "HYPO", "GENE", "WHOLE")), detail), "DELTAS")
})

test_that("keys with no marker is refused rather than passing vacuously", {
  empty <- data.frame(MARKER = character(), FIGURE = character(),
                      AREA = character(), SUBAREA = character(),
                      stringsAsFactors = FALSE)
  expect_error(SEMseeker:::meta_inference_files_check(
    data.frame(STUDY = "A", STUDY_FOLDER = ".", stringsAsFactors = FALSE),
    empty, .inference_detail()), "no marker")
})
