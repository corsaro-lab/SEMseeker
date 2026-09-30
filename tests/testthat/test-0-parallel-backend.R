# The statistical models are parallelised: assoc_apply_stat_model() runs one
# foreach %dorng% per area, and that is the loop the whole design is built
# around. Nothing ever asserted that the loop actually runs anywhere but here.
#
# The gap was known before this file existed; what was missing was the
# assertion, not the awareness of needing one.
#
# Every parallelism defect this package has had was silent in exactly that way.
# multicore on Windows was accepted and degraded to a single worker with no
# message. multicore on Linux forked a process holding a native thread pool and
# the children waited forever on a lock their parent's threads had taken. In
# both cases the suite stayed green: the results were correct, they were simply
# computed by one process while the caller believed otherwise.
#
# WHAT THIS ASSERTS, AND WHY NOT TIME.
# A timing assertion would be the obvious thing and the wrong one: two cores on
# a shared runner, a start-up cost per worker, and a margin small enough that
# the test would fail for weather. The mechanism, on the other hand, is
# deterministic. If the bodies of the loop report more than one process id, the
# work was distributed. If they all report this process, it was not. No
# threshold, no stopwatch, nothing to tune.
#
# WHY THE BACKEND AND NOT assoc_apply_stat_model() ITSELF.
# The model loop has no seam through which a body can report which worker ran
# it, and adding one would mean changing production code to satisfy a test. The
# backend under test here is the same one that loop uses: core_init_env()
# registers it once via core_parallel_session(), and every %dorng% in the
# package resolves through it. If it distributes, they distribute.

test_that("the backend behind every %dorng% distributes work across processes", {
  skip_on_cran()
  skip_if_not_installed("doRNG")
  skip_if_not_installed("foreach")
  skip_if_not_installed("future")
  skip_if(is.na(future::availableCores()) || future::availableCores() < 2L,
          "a single core has nothing to distribute")

  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  on.exit({
    SEMseeker:::core_close_env()
    unlink(tempFolder, recursive = TRUE)
  }, add = TRUE)

  SEMseeker:::core_init_env(tempFolder, parallel_strategy = "multisession")

  parent <- Sys.getpid()
  # More iterations than workers, so the result does not depend on how the
  # backend happens to chunk them.
  pids <- foreach::foreach(i = seq_len(16L), .combine = c) %dorng% Sys.getpid()

  expect_length(pids, 16L)
  expect_gt(length(unique(pids)), 1L)
  expect_false(all(pids == parent),
               info = "every iteration ran in the calling process: the request for parallelism was accepted and not honoured")
})

test_that("sequential runs in the calling process, so the check above can fail", {
  # A test that cannot fail proves nothing. This pins the other end: under
  # sequential every iteration must report this process, which is what makes
  # "more than one pid" evidence of something rather than an accident of how
  # foreach happens to behave.
  skip_on_cran()
  skip_if_not_installed("doRNG")
  skip_if_not_installed("foreach")

  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  on.exit({
    SEMseeker:::core_close_env()
    unlink(tempFolder, recursive = TRUE)
  }, add = TRUE)

  SEMseeker:::core_init_env(tempFolder, parallel_strategy = "sequential")

  parent <- Sys.getpid()
  pids <- foreach::foreach(i = seq_len(8L), .combine = c) %dorng% Sys.getpid()

  expect_length(pids, 8L)
  expect_true(all(pids == parent))
})

test_that("a multicore request is honoured as multisession, not as one worker", {
  # PR #46 fixed the Windows half of this: multicore was accepted there and
  # degraded to a single worker in silence. The Linux half was fixed later, when
  # forking a process that holds a native thread pool turned out to hang rather
  # than to crash. Both corrections say the same thing, and it is worth an
  # assertion: asking for multicore must leave a caller with real workers on
  # every platform, never with one process pretending.
  skip_on_cran()
  skip_if_not_installed("doRNG")
  skip_if_not_installed("foreach")
  skip_if(is.na(future::availableCores()) || future::availableCores() < 2L,
          "a single core has nothing to distribute")

  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  on.exit({
    SEMseeker:::core_close_env()
    unlink(tempFolder, recursive = TRUE)
  }, add = TRUE)

  SEMseeker:::core_init_env(tempFolder, parallel_strategy = "multicore")

  ssEnv <- SEMseeker:::core_get_session_info()
  expect_identical(ssEnv$parallel_strategy, "multisession")

  parent <- Sys.getpid()
  pids <- foreach::foreach(i = seq_len(16L), .combine = c) %dorng% Sys.getpid()
  expect_gt(length(unique(pids)), 1L)
  expect_false(all(pids == parent))
})
