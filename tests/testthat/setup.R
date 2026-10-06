# =============================================================================
# Shared setup for the whole suite, and how to size the bench that runs it.
#
# Read this before changing anything here, and before concluding that the suite
# is slow, stuck, or broken. The rules for writing a test are in
# engineering-decisions.md 5.16, the traps that cost days are in 5.18, and
# tests/testthat/test_template.R is the executable skeleton.
#
# -----------------------------------------------------------------------------
# TWO PLANES OF PARALLELISM, AND THEY DO NOT COMPOSE
# -----------------------------------------------------------------------------
# The PRODUCT parallelises areas inside a model, through foreach %dorng% on a
# future plan set by core_parallel_session(). It serves the user.
#
# The BENCH parallelises test FILES, through Config/testthat/parallel in
# DESCRIPTION, one process per file. It serves whoever is waiting.
#
# Multiply them and they divide: nine package workers inside each of ten
# concurrent files is ninety processes on ten cores, contending rather than
# computing. So every test runs SEQUENTIAL inside, and the files run in
# parallel outside. The one exception is the test whose subject IS parallelism:
# it asks for a strategy explicitly and does not read the global below.
#
# -----------------------------------------------------------------------------
# HOW MUCH MEMORY THE BENCH NEEDS, WHICH IS THE BINDING CONSTRAINT
# -----------------------------------------------------------------------------
# Not the cores. Measured on this suite:
#
#   one operation in a fresh process          ~1.2 GB
#   a process that has been working            4.0 - 4.7 GB   <- it plateaus
#   cost assumed per worker (dev/run-suite.R)  4.6 GB
#
# It plateaus rather than growing: five processes doing eighteen files each
# reach the same figure as one process doing eighty-eight. It is STRUCTURE, not
# payload - the memoised probe annotation is 226 MB per process, the native
# query engine keeps its own pool, the loaded namespaces are the same whatever
# the fixture size. Do not expect it to shrink because the data is small.
#
# Therefore, for the Docker VM:
#
#   VM memory   workers that fit   outcome observed
#   ---------   ----------------   ----------------------------------------
#     8 GB            1            nine workers were killed by the OOM killer
#    24 GB            4            four workers also killed, at 65 of 88 files
#    48 GB           10            completes
#
# A VM sized for the cores instead of the memory is what killed three runs. And
# a host is not free either: 48 GB of 64 leaves macOS compressing. 32 GB is the
# sensible daily setting, at six or seven workers.
#
# -----------------------------------------------------------------------------
# READ THE ALLOWANCE, NOT THE MACHINE
# -----------------------------------------------------------------------------
# Inside a container these disagree, and only one pair is right. Measured with
# `docker run --memory=6g --cpus=3`:
#
#   parallelly::availableCores()          3    <- honours the cgroup
#   parallel::detectCores(), nproc       10    <- reports the host
#   /sys/fs/cgroup/memory.max            6 GB  <- honours the cgroup
#   /proc/meminfo MemTotal              49 GB  <- reports the host
#
# dev/run-suite.R derives the worker count from the first of each pair. Nothing
# here should use the second.
#
# -----------------------------------------------------------------------------
# WHY A LOCAL GREEN IS NOT A CI GREEN
# -----------------------------------------------------------------------------
# Under R CMD check, R sets _R_CHECK_LIMIT_CORES_ and everything sees 2 cores,
# so CI runs two package workers and two test processes whatever the runner has.
# That is why the OOM never appears there: not a healthier environment, a
# blocked knob. CI also runs vignettes and examples, which the local bench does
# not, and on amd64 plus Windows, which the local bench does not.
#
# The bench is therefore HARSHER on memory and BLINDER on portability. Neither
# result substitutes for the other, and their timings are not comparable.
# =============================================================================

# num_rows <- 3e^6
# num_cols <- 5200
# populationMatrix <- as.data.frame(matrix(runif(num_rows * num_cols), nrow = num_rows, ncol = num_cols))
# library(future)
# options(future.globals.maxSize = 10 * 1024^3)
Sys.setenv(OBJC_DISABLE_INITIALIZE_FORK_SAFETY='YES')

rm(list = ls())
# DEBUG: trace setup.R progress so we can see, in a run log, at which step the
# annotation and its dependency chain get loaded.
#
# flush.console() is guarded. There is no console in the subprocesses testthat
# starts to run test files in parallel, and calling it there raises: the whole
# run then dies with "testthat subprocess failed to start", pointing at this
# line and saying nothing about why. The flush only matters when a human is
# watching output arrive, so it is skipped where nobody is.
.trace_step <- function(msg) {
  cat(sprintf("[SETUP-TRACE %s] %s\n",
              format(Sys.time(), "%H:%M:%OS3"), msg))
  if (interactive()) try(flush.console(), silent = TRUE)
}
.trace_step("setup.R BEGIN")
loadNamespace("future")
.trace_step("after loadNamespace(future)")
loadNamespace("stats")
.trace_step("after loadNamespace(stats)")
Sys.setenv(OBJC_DISABLE_INITIALIZE_FORK_SAFETY = 'YES')

# ── Probe features: 20k real EPIC probe IDs (imprinting DMR-aware) ────────
.trace_step("loading bundled test_master_features fixture")
utils::data("test_master_features", package = "SEMseeker", envir = environment())
.probe_features_all <- as.data.frame(test_master_features, stringsAsFactors = FALSE)
.trace_step(sprintf("loaded test_master_features: %d probes", nrow(.probe_features_all)))

# ── Beta matrix: real GSE133774 (EPIC 850k, BWS + MLID, 10 samples) ───────
# Built once by data-raw/build_test_signal_fixture.R; ships in data/.
# Single source of truth for both automated tests and vignette.
.trace_step("loading bundled test_signal_gse133774 fixture")
utils::data("test_signal_gse133774",      package = "SEMseeker", envir = environment())
utils::data("test_samplesheet_gse133774", package = "SEMseeker", envir = environment())

# Align probe_features to probes present in the beta matrix
.common <- intersect(.probe_features_all$PROBE, rownames(test_signal_gse133774))
probe_features <- .probe_features_all[.probe_features_all$PROBE %in% .common, ]
signal_data    <- as.data.frame(test_signal_gse133774[.common, , drop = FALSE])

nprobes  <<- nrow(signal_data)
nsamples <<- ncol(signal_data)
.trace_step(sprintf("aligned fixture: %d probes × %d samples", nprobes, nsamples))

# ── Sample sheet ───────────────────────────────────────────────────────────
mySampleSheet <- test_samplesheet_gse133774
# Extra columns used by association / batch tests
set.seed(474693)
mySampleSheet$Phenotest   <- stats::rnorm(nrow(mySampleSheet), mean = 1000, sd = 567)
mySampleSheet$Group       <- c(rep(TRUE,  ceiling(nrow(mySampleSheet) / 2)),
                               rep(FALSE, floor(nrow(mySampleSheet)   / 2)))
mySampleSheet$Covariates1 <- stats::rnorm(nrow(mySampleSheet), mean = 567,  sd = 1000)
mySampleSheet$Covariates2 <- stats::rnorm(nrow(mySampleSheet), mean = 67,   sd = 100)

mySampleSheet_batch <<- list(mySampleSheet, mySampleSheet, mySampleSheet)
signal_data_batch   <<- list(signal_data,   signal_data,   signal_data)

# ── IQR thresholds from Reference samples only ────────────────────────────
.ref_ids     <- unique(mySampleSheet$Sample_ID[mySampleSheet$Sample_Group == "Reference"])
.ref_ids     <- .ref_ids[.ref_ids %in% colnames(signal_data)]
.ref_signal  <- signal_data[, .ref_ids, drop = FALSE]
q1           <- apply(.ref_signal, 1, function(x) stats::quantile(x, 0.25, na.rm = TRUE))
q3           <- apply(.ref_signal, 1, function(x) stats::quantile(x, 0.75, na.rm = TRUE))
signal_medians <- apply(.ref_signal, 1, stats::median)
iqr          <- data.frame(q3 - q1)

signal_superior_thresholds <- data.frame("HIGH" = q3 + 3 * iqr)
signal_inferior_thresholds <- data.frame("LOW"  = q1 - 3 * iqr)
colnames(signal_inferior_thresholds) <- "LOW"
colnames(signal_superior_thresholds) <- "HIGH"
row.names(signal_superior_thresholds) <- probe_features$PROBE
row.names(signal_inferior_thresholds) <- probe_features$PROBE

signal_thresholds <- data.frame(
  "signal_median_values"       = signal_medians,
  "signal_inferior_thresholds" = signal_inferior_thresholds,
  "signal_superior_thresholds" = signal_superior_thresholds,
  "iqr"   = iqr,
  "q1"    = q1,
  "q3"    = q3
)
colnames(signal_thresholds) <- c("signal_median_values", "signal_inferior_thresholds",
                                  "signal_superior_thresholds", "iqr", "q1", "q3")
signal_thresholds$CHR   <- probe_features$CHR
signal_thresholds$START <- probe_features$START
signal_thresholds$END   <- probe_features$END

mySampleSheet              <<- mySampleSheet
signal_data                <<- signal_data
signal_medians             <<- signal_medians
signal_inferior_thresholds <<- signal_inferior_thresholds
signal_superior_thresholds <<- signal_superior_thresholds
signal_thresholds          <<- signal_thresholds
nsamples                   <<- nsamples
nprobes                    <<- nprobes



LESIONS_BP <<- 5000L  # bp-based window, literature-aligned default 5 kbp.
bonferroni_threshold <<- 0.1
batch_id <<- 1
iqrTimes <<- 3
# Sequential inside a test, parallel across test files.
#
# The files run in parallel now (Config/testthat/parallel in DESCRIPTION), which
# is where the independent units are: sixty of them, none dominating except one.
# Inside a test there is nothing to gain and something to lose. The tests are
# small - a Wilcoxon over all 18,089 probes of the fixture takes 1.1 s - so
# starting nine workers costs more than the work, and nine workers per file
# times ten concurrent files would put ninety processes on ten cores, where they
# would contend rather than compute.
#
# Two tests ask for something else, explicitly, because parallelism is their
# subject: the backend test and the libPaths propagation test. They pass the
# strategy as an argument and are unaffected by this default. Twenty other files
# already passed "sequential" by hand before this was the default, which is why
# the change here moves nineteen callers and not sixty.
#
# fork() is not among the options at all: it is unsafe with this package's
# native thread pool on every platform that offers it, and core_parallel_session()
# converts any request for it.
parallel_strategy <<- "sequential"

# SEMSEEKER_TEST_PARALLEL overrides the choice above, for two uses.
#
# Forcing "sequential" when a failure inside a worker arrives without a usable
# stack trace, which is the ordinary way to debug a parallel path.
#
# And comparing the same tests under two strategies by wall clock, which is the
# only way to find out whether the parallel paths actually pay off: the obvious
# candidate, user/real from proc.time(), cannot answer it, because multisession
# workers are separate R sessions and their CPU is never attributed to this
# process. Measured on the whole suite: user 69m37s against real 70m3s, a ratio
# of 1.0, while ten processes were busy a tenth of the time.
.forced <- Sys.getenv("SEMSEEKER_TEST_PARALLEL", "")
if (nzchar(.forced)) {
  parallel_strategy <<- .forced
  message("setup.R: parallel_strategy forced to '", .forced,
          "' by SEMSEEKER_TEST_PARALLEL")
}
markers <<- c("MUTATIONS","DELTAQ","DELTARQ","DELTAP","DELTARP","LESIONS")


# TODO
# core_recover session stored
#
tmp <- normalizePath(tempdir())

# One way to get a temporary folder, and only one.
#
# There used to be a shared vector of 50 pre-generated folders, read in two
# incompatible ways at once: 53 places took the first element and 50 of them then
# shortened the vector from the head, while 84 other places indexed it by fixed
# position, from 1 to 53. Three consequences, all silent.
#
# Indices 51 to 53 did not exist, so those tests passed NA as a result folder.
# A fixed index meant a different folder depending on how many tests had
# already taken one, so which directory a test touched depended on the order
# the suite happened to run in, and on any filter applied to it. And 137
# accesses over 50 folders guaranteed that unrelated tests shared directories
# without anyone deciding they should.
#
# It surfaced as a test asserting that a function had written no plots and
# finding four, left there by whoever else had been handed the same folder.
# The function was correct; the fixture was not. Only 72 of the 156
# core_init_env() calls in the suite pass start_fresh = TRUE, so the other 84
# inherit whatever the directory already contains.
#
# A fresh directory per call removes the sharing rather than managing it. Tests
# that reuse a folder deliberately - resume paths, idempotency - keep doing so:
# they hold it in a local variable, and only where that variable comes from has
# changed.
sem_test_folder <- function() {
  d <- file.path(normalizePath(tempdir()), "semseeker",
                 stringi::stri_rand_strings(1, 12, pattern = "[A-Za-z0-9]"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}


check_execution_context <- function() {
  calls <- sys.calls()
  if (any(sapply(calls, function(x) "test_file" %in% names(x)))) {
    showprogress <<- FALSE
    verbosity <<- 1
    # message("Called from testthat")
  } else {
    showprogress <<- TRUE
    verbosity <<- 4
    # message("Called from source or directly")
  }
}

check_execution_context()

# TODO
# core_recover session stored
#
tmp <- tempdir()
