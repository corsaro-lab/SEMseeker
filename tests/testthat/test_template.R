# =============================================================================
# TEMPLATE - copy this file, rename it test-<area>-<subject>.R, delete what you
# do not need. It is not run: no file called test_template.R matches testthat's
# ^test-.*\.R$ pattern, so it sits here as a skeleton rather than as a test.
#
# The rules below are in engineering-decisions.md 5.16, with the defect that
# produced each one. The traps are in 5.18. Read tests/testthat/setup.R for how
# the bench is sized. Everything here exists because its absence cost a day.
# =============================================================================


# --- 1. A SESSION, AND A FOLDER OF ITS OWN -----------------------------------
#
# sem_test_folder() returns a fresh directory per call. Never take one by index
# from a shared vector: that is how a test came to assert zero plots and find
# four, left by whoever else had been handed the same folder, and how adding a
# test file made a different file fail somewhere else.
#
# parallel_strategy is passed explicitly and is "sequential". The files run in
# parallel around this test; starting workers inside it as well would put
# ninety processes on ten cores. The only test that asks for anything else is
# the one whose subject is parallelism.

test_that("<what is true when this passes>", {
  tempFolder <- sem_test_folder()
  on.exit({
    SEMseeker:::core_close_env()
    unlink(tempFolder, recursive = TRUE)
  }, add = TRUE)

  SEMseeker:::core_init_env(tempFolder, parallel_strategy = "sequential")

  # ... arrange, act, assert ...
  expect_true(TRUE)
})


# --- 2. A SUGGESTS DEPENDENCY -------------------------------------------------
#
# skip_if_not_installed() is a promise, not an escape hatch. Whatever is named
# here MUST also be in .github/workflows/R-CMD-check.yml, in test-coverage.yml
# and in Dockerfile.ci. test-pkg_suggests_installed.R scans every call of this
# kind and fails on purpose when one of them would skip: a test that skips is a
# test that does not exist, and it still counts as green. Five defects lived
# behind silent skips for months.

test_that("<what needs the optional package>", {
  skip_if_not_installed("<package>")
  # skip_on_os("mac")   # only with a reason, and write the reason
  expect_true(TRUE)
})


# --- 3. ASSERT THE MECHANISM, NOT THE CLOCK -----------------------------------
#
# A timing assertion fails for the weather: shared cores, a start-up cost per
# worker, a margin too thin to survive a busy runner. Assert something
# deterministic instead. To show that work was distributed, count distinct
# process ids; to show it was not, count them too.
#
# And pin BOTH ends. "More than one pid" is evidence only because there is a
# case that yields exactly one; without it the assertion would pass with the
# backend broken. A test that cannot fail proves nothing.

test_that("<the mechanism happens>", {
  tempFolder <- sem_test_folder()
  on.exit({ SEMseeker:::core_close_env(); unlink(tempFolder, recursive = TRUE) }, add = TRUE)
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = "multisession")

  pids <- foreach::foreach(i = seq_len(16L), .combine = c) %dorng% Sys.getpid()
  expect_gt(length(unique(pids)), 1L)
})

test_that("<and the case that must NOT happen, so the check above can fail>", {
  tempFolder <- sem_test_folder()
  on.exit({ SEMseeker:::core_close_env(); unlink(tempFolder, recursive = TRUE) }, add = TRUE)
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = "sequential")

  pids <- foreach::foreach(i = seq_len(8L), .combine = c) %dorng% Sys.getpid()
  expect_true(all(pids == Sys.getpid()))
})


# --- 4. EXPECTED VALUES COME FROM OUTSIDE THE PRODUCER ------------------------
#
# Never compute what you expect by running the thing under test, or by
# reimplementing it beside the test. Such a value agrees with the code even
# when the code is wrong: it hid a mask that was not restricting anything,
# because both sides omitted the same filter.
#
# Write the expected value by hand, take it from an independent source, or
# assert a property that the result must have whatever the implementation does.

test_that("<the result has the property it must have>", {
  # expected <- 3L                      # by hand, from the fixture
  # expect_identical(result, expected)
  expect_true(TRUE)
})


# --- 5. A THRESHOLD NEEDS BOTH ENDS MEASURED ----------------------------------
#
# If a test compares against a number, that number comes from measuring the
# good case AND the bad case, not from extrapolating one of them. The symbol
# mapper's threshold holds because both ends were measured: 0.02 per cent of
# keys unmapped when the key type is right, 94 per cent when it is wrong, the
# line drawn at 20 in between. A memory coefficient taken from a single
# observation was wrong by more than a factor of two and killed three runs.


# --- 6. NOTHING THAT ASSUMES A TERMINAL ---------------------------------------
#
# There is no console in the subprocesses testthat starts for parallel files.
# flush.console() in setup.R killed every parallel run with "testthat
# subprocess failed to start", pointing at the line and saying nothing about
# why. The same goes for anything that expects interactive input, a display, or
# a working directory it did not create.


# --- 7. READ THE ALLOWANCE, NOT THE MACHINE -----------------------------------
#
# If a test cares how many cores or how much memory it has, ask for the limit
# the process was given, not for the hardware. Inside a container the two
# disagree: with --cpus=3, parallelly::availableCores() says 3 while
# detectCores() and nproc both say 10; with --memory=6g, the cgroup file says
# 6 GB while /proc/meminfo says 49.


# --- 8. INDEPENDENCE ----------------------------------------------------------
#
# No test may rely on another having run, on the order files are executed in,
# or on state surviving between them. The files run in parallel and in
# arbitrary order, and a filtered run executes a subset. If a test needs
# something built, it builds it.
