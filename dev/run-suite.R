# ---------------------------------------------------------------------------
# Run the test suite against the INSTALLED package and print explicit totals.
#
#   Rscript dev/run-suite.R                       # whole suite
#   Rscript dev/run-suite.R 7-association_analysis  # one file, to iterate
#
# Why a file and not a -e one-liner: in PowerShell a double-quoted string
# interpolates $passed, $failed, $error, so `sum(df$passed)` reaches R as
# `sum(df)` and the run ends in a parse error after the tests have already
# taken their time. A script has no quoting layer to get wrong.
#
# Why load_package = "installed" and not devtools::test(): under load_all a
# helper file's closure is package:SEMseeker, whose environment chain reaches
# neither the global environment nor setup.R's, so the two harnesses do not
# measure the same thing - and the one that skips tests silently is the one that
# hides defects. See AI-254.
#
# Why explicit totals and not the default reporter: the totals come from the
# returned data frame, so they cannot be lost to a truncated console.
# ---------------------------------------------------------------------------

filter <- commandArgs(trailingOnly = TRUE)
filter <- if (length(filter) && nzchar(filter[1])) filter[1] else NULL

## Peak memory of this process, read at the end from two sources that measure
## different things. Both are reported because either one alone misleads.
##
##  * gc() tracks the R heap only. It is blind to what a native engine allocates
##    outside it, which here is most of the interesting memory.
##  * VmHWM in /proc/self/status is the high-water mark of resident memory for
##    this process, which is the number the OOM killer acts on. Linux only.
##
## NEITHER counts the multisession workers: they are separate processes, so a
## suite that looks modest here can still take the machine down. The container
## total is the only figure that covers everything, and it has to be read from
## outside. Measured on this suite, the container sat at 5.7 GiB of a 7.7 GiB
## limit while the R heap was a fraction of that.
invisible(gc(reset = TRUE))

.peak_rss_mb <- function() {
  f <- "/proc/self/status"
  if (!file.exists(f)) return(NA_real_)
  l <- grep("^VmHWM:", readLines(f, warn = FALSE), value = TRUE)
  if (!length(l)) return(NA_real_)
  as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", l[1])) / 1024
}

## Which code produced these numbers. A commit alone does not answer that when
## the working tree has moved on, and a benchmark taken from a dirty tree that
## says only "abc1234" is a number nobody can reproduce, so the state is
## recorded next to the hash rather than assumed clean.
.git <- function(...) tryCatch(
  suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = FALSE)),
  error = function(e) character(0))
## Outside a repository, or with no git on PATH, system2() comes back either
## empty or with a zero-length vector depending on which of the two it is, so
## both have to be treated as "unknown" - a length-1 "" would otherwise print as
## a blank commit, which reads like a bug in the report rather than like a run
## taken outside version control.
.first <- function(x) if (length(x) && nzchar(x[1])) x[1] else "unknown"
.commit <- .first(.git("rev-parse", "--short", "HEAD"))
.branch <- .first(.git("rev-parse", "--abbrev-ref", "HEAD"))
.dirty  <- any(nzchar(.git("status", "--porcelain")))

.t0 <- Sys.time()
res <- testthat::test_dir("tests/testthat",
                          package        = "SEMseeker",
                          load_package   = "installed",
                          filter         = filter,
                          # "progress", not "silent". The return value is the
                          # same either way - the reporter only controls what is
                          # printed - so silence bought nothing and cost the
                          # ability to see where a run is. A suite that prints
                          # nothing for an hour cannot be told apart from one
                          # that has hung, and when a run is killed the silent
                          # reporter takes every result with it: the summary is
                          # assembled at the end and there is no end.
                          # Both failure modes happened here, on the same day.
                          reporter       = "progress",
                          stop_on_failure = FALSE)

df <- as.data.frame(res)

## The report goes to a FILE as well as to the console, under dev/test-reports/,
## named after the platform. On a VM the repository is a shared folder, so the
## file lands on the host and can be read there directly - no copying a console
## buffer between two machines, and nothing lost to a window that scrolled.
report_dir <- file.path("dev", "test-reports")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
report <- file.path(report_dir,
                    paste0("suite-", gsub("[^A-Za-z0-9]+", "-", R.version$platform),
                           if (!is.null(filter)) paste0("-", gsub("[^A-Za-z0-9]+", "-", filter)) else "",
                           ".txt"))

emit <- function(...) {
  line <- paste0(...)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = report, append = TRUE)
}

if (file.exists(report)) unlink(report)

.wall <- as.numeric(difftime(Sys.time(), .t0, units = "secs"))

emit("================================================")
emit(R.version.string, " / ", R.version$platform)
emit("run at        : ", format(.t0, "%Y-%m-%d %H:%M:%S"))
emit("commit        : ", .commit, if (.dirty) " (DIRTY: uncommitted changes)" else "",
     "   on ", .branch)
if (!is.null(filter)) emit("filter        : ", filter)
emit("tests         : ", nrow(df), " in ", length(unique(df$file)), " files")
emit("wall clock    : ", sprintf("%.0f s  (%.1f min)", .wall, .wall / 60))
emit("PASS          : ", sum(df$passed))
emit("FAIL          : ", sum(df$failed))
emit("ERROR         : ", sum(df$error))
emit("SKIP          : ", sum(df$skipped))

## Finding the peak in gc() is fiddlier than it looks, in two ways that both
## fail silently. Several columns are named "(Mb)", so indexing by that name
## picks the first and reports CURRENT use as the peak. And the column layout is
## not fixed: "limit (Mb)" is present only when a memory limit is set, so a
## position that is right on one machine is off by one on another. Observed
## here: used / (Mb) / gc trigger / (Mb) / limit (Mb) / max used / (Mb), where
## taking column 6 yields the peak in CELLS, about 2.5 million, which prints as
## a plausible and entirely wrong number of megabytes.
## Locate "max used" by name and take the column after it.
.gc_peak_mb <- function() {
  g <- gc()
  i <- which(colnames(g) == "max used")
  if (!length(i) || i[1] >= ncol(g)) return(NA_real_)
  sum(g[, i[1] + 1L])
}

.heap <- .gc_peak_mb()
if (is.finite(.heap))
  emit("R heap peak   : ", sprintf("%.0f MB", .heap),
       "   (this process, R objects only)")
## Peak resident memory is Linux-only, and it says so rather than vanishing.
## There is no portable R API for it: /proc/self/status carries VmHWM on Linux,
## macOS has no /proc and `ps` reports current RSS rather than the high-water
## mark, and utils::memory.size() is defunct. Reporting the current figure as a
## peak would be worse than reporting nothing, so the row declares the gap and
## the R heap line above stays as the number that means the same thing
## everywhere.
.rss <- .peak_rss_mb()
emit("process peak  : ",
     if (is.finite(.rss))
       paste0(sprintf("%.0f MB", .rss),
              "   (resident, this process; workers not included)")
     else
       "not available on this platform (Linux only; the R heap line is portable)")
emit("================================================")

## ---------------------------------------------------------------------------
## Where the time goes.
##
## testthat already times every test and as.data.frame() already carries the
## numbers; the report just never printed them, so the only way to find a slow
## test was to watch the run. That is how an initialisation costing 287 s per
## test file stayed invisible while it dominated the suite.
##
## Wall clock is the honest column here. Do NOT read user/real as a measure of
## parallelism: proc.time() counts the CPU of this process and of children it
## WAITS for, and multisession workers are separate R sessions talking over
## sockets, so their CPU is never attributed. The ratio stays near 1 whether the
## workers are saturating the machine or not; the whole suite reported user
## 69m37s against real 70m3s while 10 processes were busy a tenth of the time.
## To find out whether parallelism works, compare wall clock between strategies
## (dev/parallel-ab.R) instead.
## ---------------------------------------------------------------------------
if ("real" %in% names(df) && any(is.finite(df$real))) {
  total <- sum(df$real, na.rm = TRUE)
  emit("")
  emit("Slowest tests (wall clock; total ", sprintf("%.0f", total), " s)")

  ord <- order(df$real, decreasing = TRUE, na.last = NA)
  for (i in utils::head(ord, 15))
    emit(sprintf("  %7.1f s  %5.1f%%  %s :: %s",
                 df$real[i], 100 * df$real[i] / total, df$file[i], df$test[i]))

  by_file <- tapply(df$real, df$file, sum, na.rm = TRUE)
  by_file <- sort(by_file, decreasing = TRUE)
  emit("")
  emit("Slowest files")
  for (n in utils::head(names(by_file), 15))
    emit(sprintf("  %7.1f s  %5.1f%%  %s",
                 by_file[[n]], 100 * by_file[[n]] / total, n))
  emit("================================================")
}

bad_idx <- which(df$failed > 0 | df$error)
if (length(bad_idx)) {
  emit("")
  emit("Not green:")
  for (i in bad_idx)
    emit("  ", df$file[i], " :: ", df$test[i],
         "  (failed=", df$failed[i], " error=", df$error[i], ")")

  ## The counts say WHERE, the messages say WHAT. Without them a report from
  ## another machine only tells you that something is wrong, which is the
  ## position we were in with a CI job that died leaving no log.
  emit("")
  emit("================ messages ================")
  for (i in bad_idx) {
    emit("")
    emit("--- ", df$file[i], " :: ", df$test[i])
    for (r in res[[i]]$results) {
      if (inherits(r, c("expectation_failure", "expectation_error"))) {
        emit("  [", class(r)[1], "] ",
             gsub("\n", "\n      ", conditionMessage(r)))
        srcref <- r$srcref
        if (!is.null(srcref))
          emit("      at line ", srcref[1])
      }
    }
  }
  emit("")
  emit("==========================================")
} else {
  emit("")
  emit("Nothing failed.")
}

cat("\nReport written to: ", report, "\n", sep = "")
