# ---------------------------------------------------------------------------
# Windows VM - one-off setup, then the fast loop
#
# Why this exists. The Windows CI job takes 4-5 hours and, when it dies, can
# leave nothing behind: run 31852499666 failed with no failed step and no log
# (`BlobNotFound`). A local Windows VM does not make a full R CMD check faster -
# on Apple Silicon R runs under x64 emulation, so it is slower - but it gives
# the two things the CI cannot: the output, and the ability to re-run ONE test
# file in minutes instead of waiting for a whole matrix.
#
# Host: Apple Silicon -> Parallels or VMware Fusion with Windows 11 ARM.
# VirtualBox does not run Windows on ARM. R for Windows is x64 and runs under
# the OS emulation layer; for the defects Windows actually produces - path
# separators, case-insensitive file systems, file locking, encoding, `parallel`
# without fork - that is equivalent to the x64 runner. For a genuinely
# architecture-specific defect it is not, and only the CI can tell you.
#
# Inside the VM, with the repo reachable (Parallels shares the Mac home under
# \\Mac\Home), run this file once with:
#
#   Rscript dev\windows-vm-setup.R
#
# after dev\windows-vm-bootstrap.ps1 has put Rscript on PATH. Note that winget
# installs R under "C:\Program Files (x86)\R" on this host, which is why the
# bootstrap searches for Rscript.exe instead of predicting where it went.
#
# ---------------------------------------------------------------------------

options(repos = c(CRAN = "https://cloud.r-project.org"))

## Binaries, not sources. On Windows almost everything is available prebuilt,
## and compiling is where the host bites back: winget installs the toolchain
## that matches the OS (aarch64 on Apple Silicon) while R here is the x64 build
## running under emulation, so a source build fails with
##   clang-19: error: unsupported option '-msse2' for target 'aarch64-w64-mingw32'
## Preferring binaries sidesteps the mismatch for every package that has one; a
## package with no binary still needs an x64 Rtools (see the bootstrap script).
options(pkgType = "binary", install.packages.check.source = "no")

message("R: ", R.version.string, " on ", R.version$platform)
if (.Platform$OS.type != "windows")
  stop("This script is meant to be run inside the Windows VM.")

## 1. A writable library ----------------------------------------------------
## The default library lives under Program Files and needs administrator rights,
## so a plain `install.packages()` fails with "'lib = ...' is not writable".
## R's own answer to this is a personal library; create it rather than asking the
## user to run everything elevated, which would also install packages where the
## next non-elevated session cannot see them.
personal <- Sys.getenv("R_LIBS_USER")
if (!nzchar(personal) || personal == "NULL")
  personal <- file.path(Sys.getenv("USERPROFILE"), "R",
                        paste0("win-library-", getRversion()[, 1:2]))
dir.create(personal, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(personal, .libPaths()))
message("Library: ", personal)

## 2. Toolchain -------------------------------------------------------------
## Rtools builds the few packages with no Windows binary. It is installed but
## not necessarily on PATH: on Windows ARM winget lays down rtools45-aarch64
## under C:\rtools45. Find it and put its bin directories in front, for this
## session and for the user.
if (!nzchar(Sys.which("make"))) {
  roots <- Sys.glob(c("C:/rtools*", "C:/Program Files/Rtools*",
                      "C:/Program Files (x86)/Rtools*"))
  bins  <- unlist(lapply(roots, function(r)
    Filter(dir.exists, file.path(r, c("usr/bin", "x86_64-w64-mingw32.static.posix/bin",
                                      "aarch64-w64-mingw32.static.posix/bin", "mingw64/bin")))))
  if (length(bins)) {
    Sys.setenv(PATH = paste(c(bins, Sys.getenv("PATH")), collapse = .Platform$path.sep))
    message("Rtools added to PATH: ", paste(bins, collapse = "; "))
  }
}
if (!nzchar(Sys.which("make")))
  message("NOTE: still no `make` on PATH. Most dependencies are Windows binaries ",
          "and will install anyway; only a package with no binary would fail.")

## 3. Package managers ------------------------------------------------------
for (p in c("remotes", "BiocManager", "devtools", "testthat"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, lib = personal)

## 4. Dependencies ----------------------------------------------------------
## The long pole, but it is paid ONCE: on Windows these are mostly binaries.
## BiocManager is used for everything so the Bioconductor annotation packages
## resolve; limma goes first for the reason recorded in the CI workflow - the
## Illumina annotation packages load minfi while installing, and minfi needs it.
if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma", ask = FALSE, update = FALSE, lib = personal,
                       type = "binary")

## The same explicit list the CI installs, and for the same reason: the suite
## exercises Suggests that `dependencies = "Suggests"` would drag in along with
## heavy chains (pathfindR -> ggkegg, ctdR) that are not needed and that have
## repeatedly broken CI. Enumerating them is what keeps the VM comparable to the
## runner - and their absence is not loud: without the Illumina annotation the
## technology detection falls back to a row-count heuristic, the coverage gate
## answers "skipped" instead of pass/fail, and a dozen pipeline tests fail for
## what looks like a Windows defect and is not.
## type = "both", NOT "binary". These live in Bioconductor's annotation
## repository, which does not serve them the way `type = "binary"` looks for
## them, so forcing it installs nothing and says nothing - the first attempt here
## left all four missing while reporting success. They carry no compiled code, so
## a source install is just unpacking and needs no toolchain: the binary-only
## rule that limma needs is exactly wrong for them.
anno_pkgs <- c("IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
               "IlluminaHumanMethylation450kanno.ilmn12.hg19",
               "IlluminaHumanMethylation27kanno.ilmn12.hg19",
               "AnnotationHub",
               "TxDb.Hsapiens.UCSC.hg19.knownGene")
anno_missing <- anno_pkgs[!vapply(anno_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(anno_missing)) {
  message("Annotation packages to install: ", paste(anno_missing, collapse = ", "))
  BiocManager::install(anno_missing, ask = FALSE, update = FALSE, lib = personal,
                       type = "both")
}

cran_pkgs <- c("Rfast", "mockery", "ggpubr")
cran_missing <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(cran_missing)) install.packages(cran_missing, lib = personal)

deps <- remotes::local_package_deps(dependencies = TRUE)
missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing ", length(missing), " dependencies: ",
          paste(utils::head(missing, 10), collapse = ", "),
          if (length(missing) > 10) ", ..." else "")
  BiocManager::install(missing, ask = FALSE, update = FALSE, lib = personal,
                       type = "binary")
}

## polars is not on CRAN and needs NOT_CRAN, same as the workflow does.
Sys.setenv(NOT_CRAN = "true")
if (!requireNamespace("polars", quietly = TRUE))
  install.packages("polars", lib = personal,
                   repos = "https://community.r-multiverse.org")

message("\nSetup done. The loop from here on is:\n\n",
        "  devtools::install('.', quick = TRUE, upgrade = FALSE, build = FALSE)\n",
        "  testthat::test_dir('tests/testthat', package = 'SEMseeker',\n",
        "                     load_package = 'installed', stop_on_failure = FALSE)\n\n",
        "and to iterate on one file, add filter = '7-association_analysis'.\n\n",
        "Use load_package = 'installed' and not devtools::test(): under load_all\n",
        "a helper's closure is package:SEMseeker, whose environment chain reaches\n",
        "neither the global environment nor setup.R's, so the two harnesses do not\n",
        "measure the same thing (see AI-254).\n")

## 5. Verify, and say so in a file -------------------------------------------
## A setup that half-worked is worse than one that failed: the suite still runs,
## the technology detection falls back to a row-count heuristic, the coverage
## gate answers "skipped" instead of pass/fail, and a dozen tests fail in ways
## that look like defects of this platform. So check what is actually loadable
## and write the answer where the host can read it.
required <- c("limma",
              "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
              "IlluminaHumanMethylation450kanno.ilmn12.hg19",
              "IlluminaHumanMethylation27kanno.ilmn12.hg19",
              "AnnotationHub", "TxDb.Hsapiens.UCSC.hg19.knownGene",
              "Rfast", "mockery", "ggpubr", "polars", "devtools", "testthat")

ok <- vapply(required, requireNamespace, logical(1), quietly = TRUE)

report_dir <- file.path("dev", "test-reports")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
report <- file.path(report_dir, "setup-check.txt")
if (file.exists(report)) unlink(report)

emit <- function(...) {
  line <- paste0(...)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = report, append = TRUE)
}

emit("================================================")
emit(R.version.string, " / ", R.version$platform)
emit("run at   : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
emit("library  : ", personal)
emit("================================================")
for (p in names(ok))
  emit(if (ok[p]) "  OK      " else "  MISSING ", p)
emit("================================================")
if (all(ok)) {
  emit("All present. The suite will measure what CI measures.")
} else {
  emit("INCOMPLETE. Until these are installed the suite fails in ways that")
  emit("look like platform defects and are not - see the annotation notes above.")
}
