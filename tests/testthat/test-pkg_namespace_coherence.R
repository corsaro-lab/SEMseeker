# Meta-test: NAMESPACE coherence with R/ source.
#
# Catches the failure mode that broke develop @196e62b on three CI jobs
# (BiocCheck / R-CMD-check / test-coverage): `export(marker_distribution_info)`
# survived in NAMESPACE after the function file was moved to bench/ (which
# is excluded from the package build). The package built but failed at
# "testing if installed package can be loaded" with
#
#     undefined exports: marker_distribution_info
#
# Running this at test time fails *before* CI ever sees the install error.

.parse_namespace_exports <- function(ns_path) {
  txt <- readLines(ns_path, warn = FALSE)
  out <- regmatches(txt, regexec("^export\\(([^)]+)\\)\\s*$", txt))
  vapply(Filter(function(m) length(m) == 2L, out), `[[`, character(1L), 2L)
}

.parse_r_function_defs <- function(r_dir) {
  # Use R's own parser so multi-line definitions, comments, weird formatting,
  # or `name <-\n  function(...)` are all picked up correctly. A regex on
  # `^name <- function` misses e.g. `ss_analysis <-\n  function(...)`.
  files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  defs <- character(0L)
  for (f in files) {
    exprs <- tryCatch(parse(f), error = function(e) NULL)
    if (is.null(exprs)) next
    for (e in as.list(exprs)) {
      if (length(e) >= 3L &&
          is.symbol(e[[1L]]) &&
          as.character(e[[1L]]) %in% c("<-", "=", "<<-")) {
        nm <- e[[2L]]
        if (is.symbol(nm)) defs <- c(defs, as.character(nm))
      }
    }
  }
  unique(defs)
}

test_that("every NAMESPACE export(X) has X defined in R/ (no stale exports)", {
  # Dev-time canary. Requires the R/ SOURCE files (.R) to be visible, which
  # is true under devtools::load_all() but NOT under R CMD check / covr —
  # the installed package's R/ contains compiled .rdb/.rdx, not sources. So
  # gate on the presence of at least one .R source file, not just dir
  # existence.
  pkg_root <- testthat::test_path("..", "..")
  ns_path  <- file.path(pkg_root, "NAMESPACE")
  r_dir    <- file.path(pkg_root, "R")
  r_sources <- if (dir.exists(r_dir))
                 list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
               else character(0L)

  testthat::skip_if_not(
    file.exists(ns_path) && length(r_sources) > 0L,
    "R/ source files not visible (installed-package or non-source layout) — dev-time only"
  )

  exports <- .parse_namespace_exports(ns_path)
  defs    <- .parse_r_function_defs(r_dir)

  stale <- setdiff(exports, defs)
  testthat::expect_equal(
    length(stale), 0L,
    info = sprintf(
      "stale export(X) entries (X not defined in R/): %s",
      paste(stale, collapse = ", ")
    )
  )
})

test_that("every export(X) in NAMESPACE resolves via getExportedValue() on the loaded package", {
  # Works in BOTH dev (load_all) and installed contexts. In installed
  # context NAMESPACE lives inside the installed package tree; fall back to
  # system.file() when the source NAMESPACE path isn't visible.
  testthat::skip_if_not_installed("SEMseeker")
  ns_path <- testthat::test_path("..", "..", "NAMESPACE")
  if (!file.exists(ns_path)) {
    ns_path <- system.file("NAMESPACE", package = "SEMseeker")
  }
  testthat::skip_if_not(nzchar(ns_path) && file.exists(ns_path),
                        "NAMESPACE not found in source tree or installed package")

  exports <- .parse_namespace_exports(ns_path)
  unresolved <- character(0L)
  for (name in exports) {
    ok <- tryCatch({
      getExportedValue("SEMseeker", name)
      TRUE
    }, error = function(e) FALSE)
    if (!ok) unresolved <- c(unresolved, name)
  }
  testthat::expect_equal(
    length(unresolved), 0L,
    info = sprintf(
      "exports declared in NAMESPACE but not resolvable on the loaded namespace: %s",
      paste(unresolved, collapse = ", ")
    )
  )
})
