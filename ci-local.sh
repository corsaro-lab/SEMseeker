#!/usr/bin/env bash
# ─── Run CI check locally in Docker ──────────────────────────────────────────
# Usage:
#   ./ci-local.sh                 → full R CMD check, output to stdout
#   ./ci-local.sh check           → same as above
#   ./ci-local.sh coverage        → run covr::package_coverage() (mirrors test-coverage CI)
#   ./ci-local.sh suite           → install the package, then run the suite against the
#                                  INSTALLED build with explicit totals - the same
#                                  measurement macOS and the Windows VM take, so the
#                                  three are comparable (see dev/run-suite.R)
#   ./ci-local.sh suite <filter>  → same, one test file
#   ./ci-local.sh logs            → R CMD check + saves .Rcheck dir to ./ci-check-output/
#   ./ci-local.sh shell           → interactive container for debugging
#   ./ci-local.sh build           → (re)build image only, keep cache
#   ./ci-local.sh rebuild         → force full rebuild (no cache)
#
# ─── Native R (no Docker) — for macOS-specific bugs ──────────────────────────
#   ./ci-local.sh native          → smoke: load_all + source setup.R
#                                   Fastest check (~5s if deps cached). Use to
#                                   verify the package loads and setup.R runs
#                                   without crashing — catches load-time bugs
#                                   that are mac-arm64-specific (tcltk segfault,
#                                   onAttach chains, etc.).
#   ./ci-local.sh native test_check → full testthat::test_check("SEMseeker")
#                                     (slow, mirrors what R CMD check tests does)
#   ./ci-local.sh native <test-file>  → load_all + source setup.R + test_file
#                                     Example: ./ci-local.sh native tests/testthat/test-5-annotation-concordance.R
#   All native modes set NOT_CRAN=true so skip_on_cran() does not fire,
#   matching the GitHub Actions matrix.
#   Requires: devtools, testthat, polars, and the package's Imports
#   available in the local user library.
#
# Error collection:
#   - Plain stdout:   ./ci-local.sh 2>&1 | tee ci-local.log
#   - Full Rcheck:    ./ci-local.sh logs   → writes ./ci-check-output/
#     The directory contains semseeker.Rcheck/ with all log files:
#       00check.log        → full R CMD check output
#       00install.out      → package installation log
#       semseeker/tests/testthat.Rout → all test output (every test_that block)
#   - Coverage:       ./ci-local.sh coverage
#     On failure prints the full .Rout.fail content so you can see which test failed.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

IMAGE="semseeker-ci"
MODE="${1:-check}"
OUTDIR="$(pwd)/ci-check-output"

_build() {
  docker build -f Dockerfile.ci -t "$IMAGE" "$@" .
}

case "$MODE" in

  build)
    echo "==> Building $IMAGE (cached) ..."
    _build
    ;;

  rebuild)
    echo "==> Force-rebuilding $IMAGE (no cache) ..."
    _build --no-cache
    ;;

  shell)
    echo "==> Building $IMAGE (if needed) ..."
    _build -q
    echo "==> Starting interactive shell in /pkg ..."
    echo "    Useful commands inside:"
    echo "      R CMD check --no-manual --as-cran /pkg"
    echo "      Rscript -e 'devtools::load_all(); testthat::test_file(\"tests/testthat/test-7-association_analysis.R\")'"
    docker run --rm -it -w /pkg "$IMAGE" bash
    ;;

  logs)
    echo "==> Building $IMAGE (if needed) ..."
    _build -q
    echo ""
    echo "==> Running R CMD check, collecting artefacts in $OUTDIR ..."
    mkdir -p "$OUTDIR"

    # Write the R script to a temp file to avoid shell-escaping issues
    RSCRIPT_TMP="/tmp/ci-check-$$.R"
    cat > "$RSCRIPT_TMP" <<'REOF'
res <- rcmdcheck::rcmdcheck(
  args       = c("--no-manual", "--as-cran"),
  build_args = "--no-manual",
  error_on   = "never",
  check_dir  = "/check-output"
)
cat("\n=== SUMMARY ===\n")
print(res)
if (length(res$errors) > 0) quit(status = 1)
REOF

    docker run --rm \
      -v "$OUTDIR":/check-output \
      -v "$RSCRIPT_TMP":/tmp/ci-check.R:ro \
      "$IMAGE" \
      Rscript /tmp/ci-check.R

    EXIT_CODE=$?
    rm -f "$RSCRIPT_TMP"

    echo ""
    echo "==> Artefacts saved to: $OUTDIR"
    echo "    Key files:"
    echo "      $OUTDIR/semseeker.Rcheck/00check.log"
    echo "      $OUTDIR/semseeker.Rcheck/00install.out"
    echo "      $OUTDIR/semseeker.Rcheck/semseeker/tests/testthat.Rout"
    exit $EXIT_CODE
    ;;

  coverage)
    echo "==> Building $IMAGE (if needed) ..."
    _build -q
    echo ""
    echo "==> Running covr::package_coverage() (mirrors test-coverage CI) ..."
    echo ""

    RSCRIPT_TMP="/tmp/ci-coverage-$$.R"
    cat > "$RSCRIPT_TMP" <<'REOF'
Sys.setenv(NOT_CRAN = "true", RENV_CONFIG_SANDBOX_ENABLED = "false")
if (!requireNamespace("covr", quietly = TRUE))
  install.packages("covr", repos = "https://packagemanager.posit.co/cran/__linux__/noble/latest")
install_path <- "/tmp/covr-pkg"
tryCatch({
  cov <- covr::package_coverage(
    path        = "/pkg",
    quiet       = FALSE,
    clean       = FALSE,
    install_path = install_path
  )
  cat("\n=== COVERAGE SUMMARY ===\n")
  print(cov)
}, error = function(e) {
  cat("\n=== COVERAGE FAILED ===\n", conditionMessage(e), "\n")
  # Print .Rout.fail if present
  fail_files <- list.files(install_path, pattern = "\\.Rout\\.fail$",
                           recursive = TRUE, full.names = TRUE)
  for (f in fail_files) {
    cat("\n=== FAIL FILE:", f, "===\n")
    cat(readLines(f), sep = "\n")
  }
  quit(status = 1)
})
REOF

    docker run --rm \
      -v "$RSCRIPT_TMP":/tmp/ci-coverage.R:ro \
      -e NOT_CRAN=true \
      -e RENV_CONFIG_SANDBOX_ENABLED=false \
      -e RENV_ACTIVATE_PROJECT=0 \
      "$IMAGE" \
      Rscript /tmp/ci-coverage.R

    EXIT_CODE=$?
    rm -f "$RSCRIPT_TMP"
    exit $EXIT_CODE
    ;;

  native)
    # Native R (no Docker) — reproduces macOS-arm64 specific bugs.
    # Use this for fast iteration when the GitHub Actions macOS runner
    # fails in a way that Docker Linux cannot reproduce.
    TARGET="${2:-smoke}"
    RSCRIPT_TMP="/tmp/ci-native-$$.R"

    case "$TARGET" in
      smoke)
        echo "==> Native smoke: load_all + setup.R + data sanity ..."
        cat > "$RSCRIPT_TMP" <<'REOF'
cat("R:", R.version$version.string, "\n")
cat("sysname:", Sys.info()[["sysname"]], "\n\n")
Sys.setenv(NOT_CRAN = "true")
cat("=== devtools::load_all('.') ===\n")
suppressMessages(devtools::load_all(".", quiet = TRUE))
cat("OK\n\n")
cat("=== source('tests/testthat/setup.R') ===\n")
source("tests/testthat/setup.R")
cat("OK\n")
cat("probe_features dim:", paste(dim(probe_features), collapse = "x"), "\n")
cat("first 3 probes:", paste(head(probe_features$PROBE, 3), collapse = ", "), "\n\n")
cat("=== data files reachable through namespace ===\n")
# Each .rda in data/ must be auto-loadable when LazyData: true.
# This catches regressions where DESCRIPTION drops LazyData and breaks
# `package::dataset` access (silently passes load_all but fails CI).
for (obj in c("metrics_properties", "cytoband_hg19", "dmr_annotation")) {
  x <- tryCatch(get(obj, envir = asNamespace("SEMseeker")),
                error = function(e) NULL)
  cat(sprintf("  %-24s %s\n", obj,
              if (is.null(x)) "MISSING (LazyData broken?)"
              else sprintf("OK (%s)", paste(class(x), collapse = "/"))))
  if (is.null(x)) quit(status = 1)
}
REOF
        ;;

      test_check)
        echo "==> Native test_check (full testthat suite) ..."
        cat > "$RSCRIPT_TMP" <<'REOF'
Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(".", quiet = TRUE))
testthat::test_check("SEMseeker", reporter = "summary")
REOF
        ;;

      *)
        # Treat any other argument as a test file path
        if [ ! -f "$TARGET" ]; then
          echo "Error: '$TARGET' is not a file. Use:" >&2
          echo "  ./ci-local.sh native                    (smoke)" >&2
          echo "  ./ci-local.sh native test_check         (full suite)" >&2
          echo "  ./ci-local.sh native <path/to/test.R>   (single file)" >&2
          exit 2
        fi
        echo "==> Native test_file: $TARGET ..."
        cat > "$RSCRIPT_TMP" <<REOF
Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/setup.R")
testthat::test_file("$TARGET", reporter = "summary")
REOF
        ;;
    esac

    EXIT_CODE=0
    Rscript "$RSCRIPT_TMP" || EXIT_CODE=$?
    rm -f "$RSCRIPT_TMP"
    exit $EXIT_CODE
    ;;

  suite)
    # AI-255. The same measurement the other two systems take, so the three are
    # comparable: the package is INSTALLED and the suite runs against it, with
    # explicit totals from the returned data frame.
    #
    # This is deliberately not `native` (load_all) and not `check` (R CMD check).
    # Under load_all a helper file's closure is package:SEMseeker, whose
    # environment chain reaches neither the global environment nor setup.R's, so
    # helpers reading a fixture as a free variable die there and pass in CI - the
    # harness silently skips whole tests rather than failing them (AI-254).
    # Three benches that differ from each other is how a defect stays hidden.
    #
    #   ./ci-local.sh suite                       whole suite
    #   ./ci-local.sh suite 7-association_analysis one file, to iterate
    echo "==> Building $IMAGE (if needed) ..."
    _build -q
    echo ""
    echo "==> Installing the package and running the suite against it ..."
    docker run --rm -w /pkg "$IMAGE" bash -lc \
      "Rscript -e 'devtools::install(\".\", quick = TRUE, upgrade = FALSE, build = FALSE)' \
       && Rscript dev/run-suite.R ${2:-}"
    ;;

  check|*)
    echo "==> Building $IMAGE (if needed) ..."
    _build
    echo ""
    echo "==> Running R CMD check ..."
    echo "    Tip: pipe through 'tee' to save: ./ci-local.sh 2>&1 | tee ci-local.log"
    echo "    Tip: use './ci-local.sh logs' to also extract the full .Rcheck directory"
    echo "    Tip: for macOS-arm64-only bugs, use './ci-local.sh native' (no Docker)"
    echo ""
    docker run --rm "$IMAGE"
    ;;

esac
