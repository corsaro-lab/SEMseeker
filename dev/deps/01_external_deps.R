# 01_external_deps.R: the external dependency graph.
#
# Answers, for every declared dependency: where does it come from (CRAN / Bioc /
# neither), how many packages does it drag in transitively, and is it actually
# used anywhere in R/, tests/ or vignettes/.
#
# Repository metadata is fetched once and cached under dev/deps/cache/. With no
# network the script falls back to the locally installed library, and records
# which mode produced the numbers.

deps_available_db <- function(paths, refresh = FALSE) {
  cache_file <- file.path(paths$cache, "available_packages.rds")
  if (!refresh && file.exists(cache_file)) {
    db <- readRDS(cache_file)
    deps_msg("  repository db: cache (%s, %d packages)",
             format(attr(db, "fetched_at"), "%Y-%m-%d"), nrow(db))
    return(db)
  }
  repos <- c(CRAN = "https://cloud.r-project.org")
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    repos <- c(repos, BiocManager::repositories()[c("BioCsoft", "BioCann", "BioCexp")])
  }
  repos <- c(repos, `r-multiverse` = "https://community.r-multiverse.org")
  db <- tryCatch(
    utils::available.packages(repos = repos, filters = list()),
    error = function(e) NULL, warning = function(w) NULL
  )
  if (is.null(db) || !nrow(db)) {
    deps_msg("  repository db: OFFLINE, falling back to the installed library")
    db <- utils::installed.packages()
    attr(db, "source") <- "installed"
  } else {
    attr(db, "source") <- "repositories"
    attr(db, "repos") <- repos
  }
  attr(db, "fetched_at") <- Sys.time()
  saveRDS(db, cache_file)
  deps_msg("  repository db: %s (%d packages)", attr(db, "source"), nrow(db))
  db
}

# Which repository does each package come from? Bioc vs CRAN matters: Bioc
# accepts dependencies from CRAN and Bioc only.
deps_origin_map <- function(db) {
  repos <- attr(db, "repos")
  if (is.null(repos) || !"Repository" %in% colnames(db)) {
    return(setNames(rep("installed", nrow(db)), rownames(db)))
  }
  repo_col <- db[, "Repository"]
  origin <- rep("other", length(repo_col))
  origin[grepl("bioconductor", repo_col, ignore.case = TRUE)] <- "Bioconductor"
  origin[grepl("cloud[.]r-project|cran", repo_col, ignore.case = TRUE)] <- "CRAN"
  origin[grepl("r-multiverse", repo_col, ignore.case = TRUE)] <- "r-multiverse"
  origin <- setNames(origin, rownames(db))
  c(origin[setdiff(names(origin), DEPS_BASE_PKGS)],
    setNames(rep("base", length(DEPS_BASE_PKGS)), DEPS_BASE_PKGS))
}

# Hard dependency closure of one package: what an install actually pulls in.
deps_closure <- function(pkg, db) {
  res <- tools::package_dependencies(
    pkg, db = db, which = c("Depends", "Imports", "LinkingTo"),
    recursive = TRUE
  )[[1]]
  if (is.null(res)) return(character())
  setdiff(res, DEPS_BASE_PKGS)
}

# Every way this package can reach for another package, in one table:
# `pkg::fn`, requireNamespace("pkg"), and bare symbols brought in by
# importFrom(). One row per call site, with the line, so 02 can attribute each
# site to the function that owns it.
#
# R/ is scanned non-recursively on purpose: R only sources files directly in
# R/, so anything in a subdirectory is dead weight, not usage.
deps_usage_scan <- function(dir, imports, pattern = "[.][Rr]$|[.]Rmd$",
                            recursive = FALSE, known_packages = character()) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE,
                      recursive = recursive)
  sym_map <- if (length(imports)) imports$symbol_to_pkg else NULL
  parts <- list()
  for (f in files) {
    if (grepl("[.]Rmd$", f)) {
      txt <- readLines(f, warn = FALSE)
      hits <- unlist(regmatches(
        txt, gregexpr("[A-Za-z][A-Za-z0-9.]*(?=:::?)", txt, perl = TRUE)))
      if (length(hits)) {
        parts[[length(parts) + 1L]] <- data.frame(
          file = basename(f), line = NA_integer_, package = hits,
          symbol = NA_character_, internal = FALSE, kind = "ns",
          stringsAsFactors = FALSE)
      }
      next
    }
    q <- deps_qualified_calls(f)
    if (!is.null(q)) { q$kind <- "ns"; parts[[length(parts) + 1L]] <- q }
    l <- deps_loader_calls(f)
    if (!is.null(l)) { l$kind <- "loader"; parts[[length(parts) + 1L]] <- l }
    b <- deps_import_symbol_calls(f, sym_map)
    if (!is.null(b)) { b$kind <- "import"; parts[[length(parts) + 1L]] <- b }
    st <- deps_string_references(f, known_packages)
    if (!is.null(st)) { st$kind <- "string"; parts[[length(parts) + 1L]] <- st }
  }
  qualified <- if (length(parts)) do.call(rbind, parts) else
    data.frame(file = character(), line = integer(), package = character(),
               symbol = character(), internal = logical(), kind = character(),
               stringsAsFactors = FALSE)
  bare <- if (any(qualified$kind == "import"))
    as.data.frame(table(package = qualified$package[qualified$kind == "import"]),
                  stringsAsFactors = FALSE) else
    data.frame(package = character(), Freq = integer(), stringsAsFactors = FALSE)
  names(bare)[2] <- "n_sites"
  list(qualified = qualified, bare = bare)
}

# Files R will never source: anything in R/ that is not a top-level .R file.
deps_r_dir_strays <- function(r_dir) {
  all_entries <- list.files(r_dir, full.names = TRUE, recursive = TRUE,
                            include.dirs = TRUE)
  top <- list.files(r_dir, pattern = "[.][Rr]$", full.names = TRUE)
  setdiff(all_entries, top)
}

run_external_deps <- function(pkg_root = ".", refresh = FALSE) {
  paths <- deps_paths(pkg_root)
  deps_msg("[01] external dependencies")

  declared <- deps_declared(paths$desc)
  imports  <- deps_namespace_imports(paths$namespace)
  db       <- deps_available_db(paths, refresh = refresh)
  origin   <- deps_origin_map(db)

  known <- declared$package
  usage_r    <- deps_usage_scan(paths$r_dir, imports, pattern = "[.][Rr]$",
                                recursive = FALSE, known_packages = known)
  usage_test <- deps_usage_scan(file.path(pkg_root, "tests"), list(),
                                pattern = "[.][Rr]$", recursive = TRUE,
                                known_packages = known)
  usage_vig  <- deps_usage_scan(file.path(pkg_root, "vignettes"), list(),
                                pattern = "[.]Rmd$", recursive = TRUE)
  strays     <- deps_r_dir_strays(paths$r_dir)

  count_files <- function(u, pkg) length(unique(u$qualified$file[u$qualified$package == pkg]))
  count_sites <- function(u, pkg) sum(u$qualified$package == pkg)
  # how we know a package is used, in decreasing order of strength:
  # ns (pkg::fn) > import (importFrom symbol) > loader (requireNamespace("x"))
  # > string (the name appears only as a string constant). A package whose only
  # evidence is "string" may well be dead: the reader has to look.
  evidence <- function(u, pkg) {
    k <- unique(u$qualified$kind[u$qualified$package == pkg])
    if (!length(k)) "" else paste(sort(k), collapse = "+")
  }

  rows <- lapply(seq_len(nrow(declared)), function(i) {
    pkg <- declared$package[i]
    closure <- deps_closure(pkg, db)
    bare_n <- if (nrow(usage_r$bare)) sum(usage_r$bare$n_sites[usage_r$bare$package == pkg]) else 0
    data.frame(
      package        = pkg,
      field          = declared$field[i],
      constraint     = declared$constraint[i],
      origin         = unname(ifelse(pkg %in% names(origin), origin[[pkg]], "unknown")),
      closure_size   = length(closure),
      r_files        = count_files(usage_r, pkg),
      r_call_sites   = count_sites(usage_r, pkg),
      ns_import_sites = bare_n,
      evidence       = evidence(usage_r, pkg),
      test_files     = count_files(usage_test, pkg),
      vignette_files = count_files(usage_vig, pkg),
      stringsAsFactors = FALSE
    )
  })
  ext <- do.call(rbind, rows)
  ext$used_in_R <- ext$r_call_sites > 0 | ext$ns_import_sites > 0
  # a name that only ever appears inside a string is not proof of a live call
  ext$used_strongly <- ext$evidence != "" & ext$evidence != "string"
  ext$used_anywhere <- ext$used_in_R | ext$test_files > 0 | ext$vignette_files > 0

  # packages used in R/ but never declared (BiocCheck / R CMD check failure)
  used_pkgs <- unique(usage_r$qualified$package)
  undeclared <- setdiff(used_pkgs, c(declared$package, DEPS_BASE_PKGS, "SEMseeker"))

  # the full install closure of Depends+Imports = what a user installs
  hard <- declared$package[declared$field %in% c("Depends", "Imports", "LinkingTo")]
  hard_closure <- unique(unlist(lapply(hard, deps_closure, db = db)))
  hard_closure <- union(hard_closure, hard)
  soft <- declared$package[declared$field == "Suggests"]
  soft_closure <- unique(unlist(lapply(soft, deps_closure, db = db)))
  full_closure <- union(hard_closure, union(soft, soft_closure))

  closure_origin <- unname(ifelse(hard_closure %in% names(origin),
                                  origin[hard_closure], "unknown"))

  deps_write_csv(ext, file.path(paths$out, "external_deps.csv"))
  deps_write_csv(
    data.frame(package = hard_closure, origin = closure_origin,
               direct = hard_closure %in% hard, stringsAsFactors = FALSE),
    file.path(paths$out, "external_closure.csv")
  )
  deps_write_csv(usage_r$qualified, file.path(paths$out, "external_call_sites.csv"))

  summary <- list(
    db_source        = attr(db, "source"),
    db_fetched_at    = format(attr(db, "fetched_at"), "%Y-%m-%d %H:%M:%S"),
    n_depends        = sum(declared$field == "Depends"),
    n_imports        = sum(declared$field == "Imports"),
    n_suggests       = sum(declared$field == "Suggests"),
    hard_closure     = length(hard_closure),
    full_closure     = length(full_closure),
    closure_by_origin = as.list(table(closure_origin)),
    imports_unused_in_R = sort(ext$package[ext$field == "Imports" & !ext$used_in_R]),
    imports_single_file = sort(setdiff(
      ext$package[ext$field == "Imports" & ext$r_files == 1], DEPS_BASE_PKGS)),
    suggests_unused     = sort(ext$package[ext$field == "Suggests" & !ext$used_anywhere]),
    non_cran_bioc       = sort(unique(hard_closure[!closure_origin %in% c("CRAN", "Bioconductor", "base")])),
    imports_string_only = sort(ext$package[ext$field == "Imports" &
                                             ext$evidence == "string"]),
    suggests_string_only = sort(ext$package[ext$field == "Suggests" &
                                              ext$evidence == "string"]),
    undeclared_used     = sort(undeclared),
    r_dir_strays        = sort(basename(strays))
  )
  jsonlite::write_json(summary, file.path(paths$out, "external_summary.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  deps_msg("  hard install closure: %d packages (%d direct)",
           length(hard_closure), length(hard))
  invisible(list(deps = ext, closure = hard_closure, summary = summary,
                 usage = usage_r))
}
