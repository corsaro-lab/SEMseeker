# start_fresh deletes only what SEMseeker created.
#
# The rule, fixed with the author before the code: a folder is deleted when it
# is missing, empty, or inside a result folder carrying Log/session_info.rds;
# a non-empty folder without that marker, and any folder that is or contains
# the working directory, the home directory or the root, are refused with
# their content untouched. Every folder here is built under tempdir(); the real
# home directory is never touched.

# Make a folder under tempdir() holding one file of the user.
# Returns the folder path and the md5 of the file, to prove it was not touched.
user_folder_with_file <- function() {
  folder <- file.path(tempfile("user_folder_"))
  dir.create(folder, recursive = TRUE)
  user_file <- file.path(folder, "thesis.txt")
  writeLines("months of work", user_file)
  list(folder = folder, user_file = user_file, md5 = unname(tools::md5sum(user_file)))
}

test_that("start_fresh on a SEMseeker result folder deletes it and recreates the session", {
  result_folder <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder, parallel_strategy = "sequential")
  SEMseeker:::core_close_env()
  stale_file <- file.path(result_folder, "Data", "stale.csv")
  dir.create(dirname(stale_file), recursive = TRUE, showWarnings = FALSE)
  writeLines("old", stale_file)

  SEMseeker:::core_init_env(result_folder, parallel_strategy = "sequential", start_fresh = TRUE)
  on.exit(SEMseeker:::core_close_env(), add = TRUE)

  expect_false(file.exists(stale_file))
  expect_true(file.exists(file.path(result_folder, "Log", "session_info.rds")))
})

test_that("start_fresh on a non-empty folder without the marker is refused and touches nothing", {
  user <- user_folder_with_file()
  on.exit(unlink(user$folder, recursive = TRUE), add = TRUE)

  expect_error(SEMseeker:::core_init_env(user$folder, parallel_strategy = "sequential", start_fresh = TRUE),
               "refusing to delete")
  expect_true(file.exists(user$user_file))
  expect_identical(unname(tools::md5sum(user$user_file)), user$md5)
})

test_that("an empty or missing folder is not an obstacle", {
  empty_folder <- tempfile("empty_folder_")
  dir.create(empty_folder)
  expect_true(SEMseeker:::core_folder_delete(empty_folder))
  expect_false(dir.exists(empty_folder))
  expect_false(SEMseeker:::core_folder_delete(tempfile("missing_folder_")))
})

test_that("the working directory and its parents are refused even with the marker", {
  parent_folder <- tempfile("parent_")
  working_folder <- file.path(parent_folder, "work")
  dir.create(file.path(working_folder, "Log"), recursive = TRUE)
  saveRDS(list(), file.path(working_folder, "Log", "session_info.rds"))
  saveRDS(list(), file.path(parent_folder, "session_marker_copy.rds"))
  dir.create(file.path(parent_folder, "Log"))
  saveRDS(list(), file.path(parent_folder, "Log", "session_info.rds"))
  on.exit(unlink(parent_folder, recursive = TRUE), add = TRUE)
  withr::local_dir(working_folder)

  expect_error(SEMseeker:::core_folder_delete(working_folder), "refusing to delete")
  expect_error(SEMseeker:::core_folder_delete(parent_folder), "refusing to delete")
  expect_true(file.exists(file.path(working_folder, "Log", "session_info.rds")))
  expect_true(file.exists(file.path(parent_folder, "session_marker_copy.rds")))
})

test_that("a folder outside the result folder is refused, as Inference pointing elsewhere would be", {
  result_folder <- tempfile("result_")
  dir.create(file.path(result_folder, "Log"), recursive = TRUE)
  saveRDS(list(), file.path(result_folder, "Log", "session_info.rds"))
  user <- user_folder_with_file()
  on.exit(unlink(c(result_folder, user$folder), recursive = TRUE), add = TRUE)

  expect_error(SEMseeker:::core_folder_delete(user$folder, result_folder = result_folder), "refusing to delete")
  expect_identical(unname(tools::md5sum(user$user_file)), user$md5)

  inference_folder <- file.path(result_folder, "Inference")
  dir.create(inference_folder)
  writeLines("previous run", file.path(inference_folder, "results.csv"))
  expect_true(SEMseeker:::core_folder_delete(inference_folder, result_folder = result_folder))
  expect_false(dir.exists(inference_folder))
})
