# 2026-06-09: meta-test: every package referenced by
# `skip_if_not_installed("...")` (or `testthat::skip_if_not_installed`)
# MUST actually be installed in the test environment.
#
# Why this exists:
#   `skip_if_not_installed(pkg)` silently skips the enclosing test if
#   `pkg` is missing. This is the correct testthat idiom for optional
#   features (Suggests deps). BUT it has a failure mode: if the CI
#   pipeline forgets to install `pkg`, the test silently skips and the
#   CI goes green WITHOUT exercising the feature. Coverage looks fine,
#   but a whole engine (limma_2, voom_2, binomial_bulk, ...) goes
#   untested.
#
#   This meta-test scans every test-*.R file under tests/testthat/,
#   extracts the package names referenced inside skip_if_not_installed()
#   calls, and asserts they are all installed. If any is missing, the
#   test FAILS with a list of (package, test file) pairs so the operator
#   knows exactly which CI install step to add.
#
#   Convention enforced: every Suggests dep used in tests MUST be in the
#   CI install list (.github/workflows/R-CMD-check.yml +
#   .github/workflows/test-coverage.yml). The CI workflows have a
#   comment naming each test that exercises each Suggests; this meta-test
#   keeps that comment in sync with reality.

test_that("every package in skip_if_not_installed(...) is actually installed", {

  test_dir <- testthat::test_path()
  test_files <- list.files(test_dir, pattern = "^test-.*\\.R$",
                           full.names = TRUE)
  # Exclude the meta-test file itself from the scan: it documents the
  # skip_if_not_installed pattern in its comments, and the regex would
  # falsely match placeholder strings like "pkg" inside that documentation.
  test_files <- test_files[basename(test_files) != "test-pkg_suggests_installed.R"]
  expect_true(length(test_files) > 0L)

  # Scan each test file for skip_if_not_installed("pkg") patterns.
  # Permissive regex: matches single or double quotes, optional
  # testthat:: prefix, optional whitespace.
  pattern <- 'skip_if_not_installed\\s*\\(\\s*["\']([^"\']+)["\']'
  hits <- list()
  for (f in test_files) {
    lines <- readLines(f, warn = FALSE)
    for (i in seq_along(lines)) {
      m <- regmatches(lines[i], regexec(pattern, lines[i]))[[1]]
      if (length(m) >= 2L) {
        pkg <- m[2]
        hits[[length(hits) + 1L]] <- list(pkg = pkg,
                                          file = basename(f),
                                          line = i)
      }
    }
  }

  if (length(hits) == 0L) {
    succeed("No skip_if_not_installed(...) calls found — nothing to assert.")
    return()
  }

  # Deduplicate per pkg, but keep the locations so we can report them.
  pkg_to_locations <- split(
    vapply(hits, function(h) sprintf("%s:%d", h$file, h$line), character(1)),
    vapply(hits, `[[`, character(1), "pkg")
  )

  missing_pkgs <- character(0)
  missing_msg  <- character(0)
  for (pkg in names(pkg_to_locations)) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_pkgs <- c(missing_pkgs, pkg)
      missing_msg  <- c(missing_msg, sprintf(
        "  %-12s  referenced by: %s",
        pkg,
        paste(pkg_to_locations[[pkg]], collapse = ", ")
      ))
    }
  }

  if (length(missing_pkgs) > 0L) {
    msg <- paste0(
      "Suggests packages referenced via skip_if_not_installed() but NOT ",
      "installed in this test environment. The tests that need them will ",
      "silently SKIP, hiding real coverage gaps from CI:\n",
      paste(missing_msg, collapse = "\n"),
      "\n\nFix: add these to .github/workflows/R-CMD-check.yml and ",
      ".github/workflows/test-coverage.yml install lists. See the comments ",
      "in those files for the convention (CRAN vs BiocManager)."
    )
    # On CI we want a hard fail so the install lists stay in sync with
    # the test references. Locally the dev may not have every Bioc dep
    # installed (e.g. a 3 GB TxDb chain); skip with a clear pointer
    # instead of blocking iteration. Override with SEMSEEKER_FULL_CHECK=1
    # to force the hard check even locally.
    in_ci   <- nzchar(Sys.getenv("CI"))
    forced  <- nzchar(Sys.getenv("SEMSEEKER_FULL_CHECK"))
    if (in_ci || forced) {
      fail(msg)
    } else {
      skip(paste0(
        "Local dev: ", length(missing_pkgs),
        " Suggests pkg(s) missing — see message below. ",
        "Set SEMSEEKER_FULL_CHECK=1 to turn this into a hard failure ",
        "locally. CI will fail anyway if not fixed.\n\n", msg
      ))
    }
  } else {
    succeed(sprintf(
      "All %d Suggests package(s) referenced in tests are installed: %s",
      length(pkg_to_locations),
      paste(sort(names(pkg_to_locations)), collapse = ", ")
    ))
  }
})

test_that("known-required CI packages are explicitly checked", {
  # Belt-and-suspenders: even if a test file is removed or refactored,
  # these packages remain required by the binomial_bulk and limma feature paths. If
  # one regresses out of CI install, this test fails fast with the
  # specific package name (rather than relying on the scan above).
  required_for_features <- c(
    "Rfast"   ,   # binomial_bulk
    "limma"       # limma_2 / voom_2 batch path
  )
  missing <- required_for_features[
    !vapply(required_for_features, requireNamespace,
            FUN.VALUE = logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0L) {
    fail(paste0(
      "Required feature packages NOT installed: ",
      paste(missing, collapse = ", "),
      ". These are needed to exercise:\n",
      "  Rfast → tests/testthat/test-assoc_glm_model_bulk.R\n",
      "  limma → tests/testthat/test-assoc_apply_stat_model_batch.R\n",
      "Add them to both CI workflows install lists."
    ))
  }
  expect_length(missing, 0L)
})
