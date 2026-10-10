# 2026-06-26: multisession/cluster workers are fresh R processes
# spawned via parallelly::makeClusterPSOCK. When the parent runs under `renv`,
# those workers must inherit the parent .libPaths() (which includes the renv
# project library) or they silently die at the first library() lookup — no
# R-visible error reaches the parent log. The fix in core_parallel_session() passes
# `rscript_libs = .libPaths()` to the multisession/cluster plans.
#
# This test asserts the contract end-to-end: after core_parallel_session() sets up a
# multisession plan, a worker future reports the SAME .libPaths() as the parent.

test_that("multisession workers inherit the parent .libPaths()", {
  skip_on_cran()

  # Minimal session env: just enough for core_parallel_session() to run. A real
  # (existing) session_folder lets core_log_event() write instead of no-op'ing,
  # but nothing here depends on the log contents.
  session_folder <- file.path(normalizePath(tempdir()),
                              paste0("ss_libpaths_", Sys.getpid()))
  dir.create(session_folder, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(session_folder, recursive = TRUE), add = TRUE)

  ssEnv <- list(
    parallel_strategy = "multisession",
    maxResources      = 50,
    session_folder    = session_folder,
    verbosity         = 1
  )
  # Populate .pkgglobalenv$ssEnv (in-memory only) so core_get_session_info() finds it.
  # Internal (non-exported) functions must be qualified with ::: so the test
  # resolves them under R CMD check (installed package), not just devtools::test().
  SEMseeker:::core_update_session_info(ssEnv, save_to_disk = FALSE)

  # Restore a clean sequential plan whatever happens.
  on.exit(future::plan(future::sequential), add = TRUE)

  SEMseeker:::core_parallel_session()

  parent_libs <- .libPaths()
  worker_libs <- future::value(future::future(.libPaths(), seed = TRUE))

  # The parent library paths must all be visible to the worker. We assert
  # set-inclusion (not strict equality): the worker may legitimately prepend
  # extra paths, but it must never LOSE the parent's — losing them is exactly
  # the renv-invisible-SEMseeker bug fixed by propagating the library paths.
  expect_true(all(parent_libs %in% worker_libs),
              info = paste0("worker .libPaths() is missing parent entries; ",
                            "rscript_libs propagation regressed."))
})
