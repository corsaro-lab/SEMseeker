# What identifies a session: the folder it was opened on, or whatever happens to
# be in memory.
#
# core_get_session_info() reads the in-memory session first and consults the
# folder only when memory is empty, so for a while the answer was "whatever is
# in memory". Opening a new analysis on a new folder inherited the previous
# one's options, and an option nobody asked for arrives with no message. It
# showed up as a test failing several files later, which is how it stayed hidden:
# it needs two analyses in one process, in that order.
#
# Both ends are pinned in each test. Asserting only that the second session sees
# the default would pass with the option never applied in the first place, so the
# first assertion is that the non-default value took effect.

test_that("a session opened on a new folder takes the default, not the previous session's value", {
  first  <- sem_test_folder()
  second <- sem_test_folder()
  on.exit({
    SEMseeker:::core_close_env()
    unlink(c(first, second), recursive = TRUE)
  }, add = TRUE)

  SEMseeker:::core_init_env(first, parallel_strategy = "sequential",
                            start_fresh = TRUE, LESIONS_BP = 1234L)
  expect_equal(as.integer(SEMseeker:::core_get_session_info()$LESIONS_BP), 1234L)
  SEMseeker:::core_close_env()

  # A different folder is a different analysis. It never asked for 1234.
  SEMseeker:::core_init_env(second, parallel_strategy = "sequential")
  expect_equal(as.integer(SEMseeker:::core_get_session_info()$LESIONS_BP), 5000L)
})

test_that("reopening the SAME folder still recovers what was set there", {
  folder <- sem_test_folder()
  on.exit({
    SEMseeker:::core_close_env()
    unlink(folder, recursive = TRUE)
  }, add = TRUE)

  SEMseeker:::core_init_env(folder, parallel_strategy = "sequential",
                            start_fresh = TRUE, LESIONS_BP = 4321L)
  SEMseeker:::core_close_env()

  # Resuming is the case the folder check must not break: same folder, so the
  # value set there is still the answer.
  SEMseeker:::core_init_env(folder, parallel_strategy = "sequential")
  expect_equal(as.integer(SEMseeker:::core_get_session_info()$LESIONS_BP), 4321L)
})

test_that("asking for a folder that never held a session gives a session, not the last one", {
  first <- sem_test_folder()
  on.exit({
    SEMseeker:::core_close_env()
    unlink(first, recursive = TRUE)
  }, add = TRUE)

  SEMseeker:::core_init_env(first, parallel_strategy = "sequential",
                            start_fresh = TRUE, LESIONS_BP = 777L)

  # A folder with no Log/session_info.rds in it. The in-memory session belongs to
  # `first`, so it is not an answer to this question.
  elsewhere <- file.path(tempdir(), paste0("no_session_", as.integer(Sys.time())))
  dir.create(elsewhere, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(elsewhere, recursive = TRUE), add = TRUE)

  recovered <- SEMseeker:::core_get_session_info(elsewhere)
  expect_true(is.null(recovered$LESIONS_BP) || length(recovered) == 0L)
})

test_that("after a close there is no session, so reading one without a folder fails", {
  folder <- sem_test_folder()
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  SEMseeker:::core_init_env(folder, parallel_strategy = "sequential",
                            start_fresh = TRUE, LESIONS_BP = 2468L)
  # Both ends pinned, as above: the value has to have taken effect before the
  # close, or the assertion after it would pass on an option that never arrived.
  expect_equal(as.integer(SEMseeker:::core_get_session_info()$LESIONS_BP), 2468L)

  SEMseeker:::core_close_env()

  # No folder named and nothing left in memory: there is no session to answer
  # with. Answering with the closed one is how a caller that forgot to
  # initialise received a stale session instead of a failure.
  expect_error(SEMseeker:::core_get_session_info(),
               "called without result folder")
})

test_that("closing a session that is already closed is not an error", {
  folder <- sem_test_folder()
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  SEMseeker:::core_init_env(folder, parallel_strategy = "sequential",
                            start_fresh = TRUE)
  SEMseeker:::core_close_env()

  # Reachable only once a close empties the session, and reached often: the
  # on.exit() idiom in these tests closes again after the body already did.
  expect_silent(SEMseeker:::core_close_env())
})
