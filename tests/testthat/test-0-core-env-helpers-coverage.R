## Coverage for previously-untested core session/env helpers.
##
## Covered:
##   core_remove_empty_folders()      recursive prune of empty subdirectories
##   .core_log_info()                 INFO logging wrapper (no-op w/o session)
##   .core_init_env_check_kwargs()    explicit-args guard (pass-through)
##   .core_apply_defaults()           apply a defaults spec onto the session
##

# ---------------------------------------------------------------------------
# core_remove_empty_folders
# ---------------------------------------------------------------------------

test_that("core_remove_empty_folders prunes empty dirs recursively but keeps non-empty ones", {
  root <- tempfile("emptyprune")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  dir.create(file.path(root, "empty1"))
  dir.create(file.path(root, "nested"))
  dir.create(file.path(root, "nested", "empty2"))
  dir.create(file.path(root, "keep"))
  writeLines("x", file.path(root, "keep", "f.txt"))

  SEMseeker:::core_remove_empty_folders(root)

  expect_false(dir.exists(file.path(root, "empty1")))
  expect_false(dir.exists(file.path(root, "nested", "empty2")))
  expect_false(dir.exists(file.path(root, "nested")))       # emptied after child removed
  expect_true (dir.exists(file.path(root, "keep")))         # holds a file -> kept
  expect_true (file.exists(file.path(root, "keep", "f.txt")))
})

# ---------------------------------------------------------------------------
# .core_log_info  /  .core_init_env_check_kwargs
# ---------------------------------------------------------------------------

test_that(".core_log_info logs without error when no session is active", {
  expect_no_error(SEMseeker:::.core_log_info("coverage smoke message"))
})

test_that(".core_init_env_check_kwargs passes explicit args through unchanged", {
  expect_equal(SEMseeker:::.core_init_env_check_kwargs(list(a = 1, b = "x")),
               list(a = 1, b = "x"))
})

# ---------------------------------------------------------------------------
# .core_apply_defaults
# ---------------------------------------------------------------------------

test_that(".core_apply_defaults writes default values into the session", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  defaults <- list(cov_test_flag = list(value = TRUE))
  args <- SEMseeker:::.core_apply_defaults(list(), defaults)

  ssEnv <- SEMseeker:::core_get_session_info()
  expect_true(isTRUE(ssEnv$cov_test_flag) || ssEnv$cov_test_flag == "TRUE")
  # the consumed default is not left dangling in the returned argument list
  expect_null(args[["cov_test_flag"]])
})
