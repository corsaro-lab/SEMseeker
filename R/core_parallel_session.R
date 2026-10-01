core_parallel_session <- function()
{
  #
  ssEnv <- core_get_session_info()
  parallel_strategy <- ssEnv$parallel_strategy

  # NOTE: `multicore` on macOS uses fork() and is known to be unsafe in
  # combination with Polars' C++ thread pool - forked children can be
  # killed by a Mach exception with no R-visible error. Tests on macOS
  # now default to `multisession` (see setup.R). End users on macOS
  # should use `multisession` or `sequential`.
  #
  # E-14: `multisession` workers are fresh R processes - .pkgglobalenv$ssEnv
  # starts empty. Every %dorng% foreach body must call
  # core_update_session_info(ssEnv) as its first statement to populate the
  # worker's namespace. See engineering-decisions.md §1.3.

  # macOS: OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES is recommended for any
  # non-sequential strategy (some packages still attempt fork internally).
  # Bioconductor disallows Sys.setenv() in package code; users on macOS
  # using a non-sequential strategy must set the env var themselves before
  # calling core_init_env(), e.g. in ~/.Renviron or Sys.setenv() at the top of
  # their script. We log a warning when the env var is missing so the
  # cause of any subsequent fork crash is obvious.
  if (Sys.info()["sysname"] == "Darwin" && parallel_strategy != "sequential") {
    env_var <- Sys.getenv("OBJC_DISABLE_INITIALIZE_FORK_SAFETY")
    if (env_var != "YES") {
      core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                " OBJC_DISABLE_INITIALIZE_FORK_SAFETY is not set to YES.",
                " On macOS with a non-sequential parallel strategy this can",
                " cause fork-related crashes. Set it in ~/.Renviron or call",
                " Sys.setenv(OBJC_DISABLE_INITIALIZE_FORK_SAFETY = 'YES')",
                " before core_init_env().")
    }
  }

  # E-13, widened: multicore means fork(), and fork() is unsafe for this package
  # on every platform that offers it. The hazard is not macOS, it is fork plus a
  # native thread pool: fork duplicates the memory of the parent's threads but
  # not the threads themselves, so a lock held at the moment of the fork stays
  # marked as taken in the child by a thread that will never exist to release
  # it. Polars builds a C++ thread pool on its first operation, so from then on
  # every child inherits the pool's queues and locks with no workers behind
  # them.
  #
  # The same cause shows two faces, which is why only one of them was guarded.
  # On macOS the system libraries detect a fork without exec and abort, so it
  # crashes loudly. On Linux nothing intervenes and the children wait forever:
  # measured on the installed package, nine children at 0.0% CPU for eleven
  # minutes on a test that takes 22 seconds under multisession.
  #
  # It is also order dependent. A fork before Polars has ever run is clean, so
  # the failure is intermittent and cannot be relied on to show up in testing.
  #
  # multicore stays accepted, because substituting it is safe and refusing it
  # would break existing scripts for nothing, but it is never honoured.
  if (parallel_strategy == "multicore") {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " multicore uses fork(), which is unsafe with this package's",
              " native thread pool: forked children crash on macOS and hang",
              " on Linux. Switching to multisession.")
    parallel_strategy <- "multisession"
    ssEnv$parallel_strategy <- parallel_strategy
  }

  # multicore requires fork(), which Windows does not provide. future does not
  # refuse the request: plan(multicore) is accepted and then degrades to a
  # single worker, so a caller who asked for parallelism gets none and is never
  # told. Convert to multisession, which every platform supports, and say so.
  # Measured on CI: the same suite took 285 minutes on Windows against 134 on
  # Linux, and the ratio was the number of workers, not the speed of the
  # platform.
  if (parallel_strategy == "multicore" && !parallelly::supportsMulticore()) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " multicore requires fork(), which this platform does not",
              " provide. Switching to multisession.")
    parallel_strategy <- "multisession"
    ssEnv$parallel_strategy <- parallel_strategy
  }

  # Nothing here may fork, so do not leave the door open for anything else in
  # the session either: enabling it only where it is available would let a
  # dependency fork on Linux for the same reasons this package must not.
  options(parallelly.fork.enable = FALSE)

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    nCore <- 2L
  } else {
    # use all cores in devtools::test()
    nCore <- future::availableCores() - 1
    nCore <- if(floor(future::availableCores() * ssEnv$maxResources/100 ) > nCore ) nCore else floor(future::availableCores() * ssEnv$maxResources/100 )
    nCore <- max(1L, nCore)  # guarantee at least 1 worker (e.g. covr subprocess with 1 core)
  }
  # permutation cluster
  outFile <- file.path(ssEnv$session_folder, "cluster_r.out")

  #
  ssEnv$parallel <- data.frame("parallel_strategy"="")
  ssEnv$parallel$parallel_strategy <- parallel_strategy
  ssEnv$parallel$nCore <- nCore

  options(readr.num_threads = nCore)  # Set to the desired number of threads

  options(doFuture.foreach.export = ".export-and-automatic-with-warning")
  # doFuture
  # is changed in onload see zzz.R
  # options(parallelly.fork.enable= FALSE)
  # options(future.globals.resolve = TRUE)
  # allow export of object of 32gb with future
  options(future.globals.maxSize= 64 * 1024^3)


  # check if future is registered
  # backend_name <- foreach::getDoParName()
  # if (!is.null(backend_name)) {
  #   core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " doFuture is registered as the %dopar% backend.")
  # } else {
  #   core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " doFuture is not registered as the %dopar% backend.")
  #   doFuture::registerDoFuture()
  # }
  doFuture::registerDoFuture()

  # get the future plan
  future_plan <- future::plan()
  core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Future Plan: ", future_plan)

  # AI-184: multisession/cluster workers are fresh R processes spawned via
  # parallelly::makeClusterPSOCK. When the parent runs under `renv`, those
  # workers start with the SYSTEM .libPaths() (or have renv reset them on
  # startup) and therefore cannot see SEMseeker or its dependencies, which
  # live in the renv project library. The result is a silent death right
  # after "I will work in multisession..." - workers fail at the first
  # library() lookup with no R-visible error in the parent log.
  #
  # Fix: capture the parent .libPaths() (which includes the renv project
  # library when renv is active) and propagate it to the workers via
  # `rscript_libs`. multicore (fork) inherits libPaths from the parent and
  # needs no special handling. See engineering-decisions.md.
  parent_libs <- .libPaths()

  # TODO: improve planning parallel management using also cluster
  if(parallel_strategy=="multisession")
  {
    future::plan( future::multisession, workers = nCore, rscript_libs = parent_libs)
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work in multisession with:", nCore, " Cores")
  }
  if(parallel_strategy=="multicore")
  {
    future::plan( future::multicore, workers = nCore)
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work in multicore with:", nCore," Cores")
  }
  if(parallel_strategy=="cluster")
  {
    if (!is.null(ssEnv$cluster_workers))
    {
      future::plan( future::cluster, workers = ssEnv$cluster_workers, rscript_libs = parent_libs)
      core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work with a cluster with:",ssEnv$cluster_workers)
    }
    else
    {
      future::plan( future::cluster, workers = nCore, rscript_libs = parent_libs)
      core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work with a cluster with:", nCore," Cores")
    }
  }

  if(parallel_strategy!="multisession" & parallel_strategy!="multicore"
    & parallel_strategy!="cluster")
  {
    options(parallelly.fork.enable= FALSE)
    future::plan(strategy = future::sequential)
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work in sequential mode")
  }

  core_update_session_info(ssEnv)

}
