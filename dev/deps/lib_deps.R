# lib_deps.R: shared helpers for the dependency analysis pipeline.
#
# Development tooling only: `dev/` is listed in .Rbuildignore, so nothing here
# ships with the package. Pure base R + igraph; no package loading, no
# evaluation of top-level package code (only function definitions are
# materialised, which has no side effects).

DEPS_BASE_PKGS <- c(
  "base", "compiler", "datasets", "grDevices", "graphics", "grid", "methods",
  "parallel", "splines", "stats", "stats4", "tcltk", "tools", "translations",
  "utils"
)

DEPS_DOMAINS <- c("sem", "assoc", "enrich", "anno", "io", "plot", "core", "util")

deps_msg <- function(...) cat(sprintf(...), "\n", sep = "")

deps_paths <- function(pkg_root = ".") {
  list(
    root      = normalizePath(pkg_root, mustWork = TRUE),
    r_dir     = file.path(pkg_root, "R"),
    desc      = file.path(pkg_root, "DESCRIPTION"),
    namespace = file.path(pkg_root, "NAMESPACE"),
    out       = file.path(pkg_root, "dev", "deps", "output"),
    cache     = file.path(pkg_root, "dev", "deps", "cache"),
    snapshots = file.path(pkg_root, "dev", "deps", "snapshots")
  )
}

deps_write_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, na = "")
  deps_msg("  wrote %s (%d rows)", basename(path), nrow(df))
  invisible(path)
}

deps_read_csv <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, na.strings = "")
}

# ---------------------------------------------------------------- DESCRIPTION

# Parse one dependency field into a data.frame of (package, version_constraint).
deps_parse_field <- function(desc, field) {
  raw <- desc[[field]]
  if (is.na(raw) || !nzchar(raw)) {
    return(data.frame(package = character(), constraint = character(),
                      field = character(), stringsAsFactors = FALSE))
  }
  parts <- trimws(strsplit(raw, ",")[[1]])
  parts <- parts[nzchar(parts)]
  pkg <- trimws(sub("\\(.*", "", parts))
  con <- ifelse(grepl("\\(", parts), trimws(gsub(".*\\(|\\).*", "", parts)), NA_character_)
  data.frame(package = pkg, constraint = con, field = field,
             stringsAsFactors = FALSE)
}

deps_declared <- function(desc_path) {
  desc <- read.dcf(desc_path)[1, ]
  fields <- c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances")
  fields <- fields[fields %in% names(desc)]
  out <- do.call(rbind, lapply(fields, function(f) deps_parse_field(desc, f)))
  out[out$package != "R", , drop = FALSE]
}

# ------------------------------------------------------------------ NAMESPACE

# importFrom(pkg, symbol) / import(pkg) directives, as a (symbol -> package) map
# plus the set of wholesale-imported packages.
deps_namespace_imports <- function(ns_path) {
  exprs <- parse(ns_path)
  sym <- list()
  whole <- character()
  exported <- character()
  for (e in exprs) {
    if (!is.call(e)) next
    head_name <- as.character(e[[1]])
    args <- as.list(e)[-1]
    lit <- function(x) if (is.character(x)) x else as.character(x)
    if (head_name == "importFrom" && length(args) >= 2) {
      pkg <- lit(args[[1]])
      for (s in args[-1]) sym[[lit(s)]] <- pkg
    } else if (head_name == "import") {
      whole <- c(whole, vapply(args, lit, character(1)))
    } else if (head_name == "export") {
      exported <- c(exported, vapply(args, lit, character(1)))
    } else if (head_name == "exportPattern") {
      # not used by this package; recorded so the caller can notice if it appears
      whole <- whole
    }
  }
  list(
    symbol_to_pkg = unlist(sym),
    whole_import  = unique(whole),
    exported      = unique(exported)
  )
}

# ------------------------------------------------------------- source parsing

# Top-level `name <- function(...)` definitions in one file, with line spans.
deps_file_functions <- function(path) {
  exprs <- tryCatch(parse(path, keep.source = TRUE),
                    error = function(e) {
                      warning(sprintf("parse failed: %s (%s)", path, conditionMessage(e)))
                      NULL
                    })
  if (is.null(exprs)) return(NULL)
  refs <- utils::getSrcref(exprs)
  rows <- list()
  for (i in seq_along(exprs)) {
    e <- exprs[[i]]
    if (!is.call(e)) next
    op <- as.character(e[[1]])[1]
    if (!op %in% c("<-", "=", "assign")) next
    if (op == "assign") next
    lhs <- e[[2]]
    rhs <- e[[3]]
    if (!is.symbol(lhs) && !is.character(lhs)) next
    if (!(is.call(rhs) && identical(as.character(rhs[[1]])[1], "function"))) next
    sr <- refs[[i]]
    rows[[length(rows) + 1L]] <- data.frame(
      fun       = as.character(lhs),
      file      = basename(path),
      line_from = if (is.null(sr)) NA_integer_ else as.integer(sr[1]),
      line_to   = if (is.null(sr)) NA_integer_ else as.integer(sr[3]),
      n_lines   = if (is.null(sr)) NA_integer_ else as.integer(sr[3] - sr[1] + 1L),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(NULL)
  attr(rows, "exprs") <- NULL
  do.call(rbind, rows)
}

# Materialise every top-level function definition of the package into one env,
# so codetools can resolve the call graph statically.
deps_function_env <- function(r_dir) {
  env <- new.env(parent = globalenv())
  failed <- character()
  for (path in list.files(r_dir, pattern = "[.][Rr]$", full.names = TRUE)) {
    exprs <- tryCatch(parse(path, keep.source = FALSE), error = function(e) NULL)
    if (is.null(exprs)) { failed <- c(failed, basename(path)); next }
    for (e in exprs) {
      if (!is.call(e)) next
      op <- as.character(e[[1]])[1]
      if (!op %in% c("<-", "=")) next
      rhs <- e[[3]]
      if (!(is.call(rhs) && identical(as.character(rhs[[1]])[1], "function"))) next
      tryCatch(eval(e, envir = env),
               error = function(err) failed <<- c(failed, basename(path)))
    }
  }
  attr(env, "failed") <- unique(failed)
  env
}

# `pkg::fn` / `pkg:::fn` occurrences, one row per call site, with the line so the
# caller can attribute them to the enclosing function.
deps_qualified_calls <- function(path) {
  pd <- tryCatch(utils::getParseData(parse(path, keep.source = TRUE)),
                 error = function(e) NULL)
  if (is.null(pd) || !nrow(pd)) return(NULL)
  idx <- which(pd$token == "SYMBOL_PACKAGE")
  if (!length(idx)) return(NULL)
  rows <- lapply(idx, function(i) {
    # the token right after `pkg` is `::` or `:::`, then the symbol
    tail_tokens <- pd[pd$line1 == pd$line1[i] & pd$col1 > pd$col1[i], , drop = FALSE]
    tail_tokens <- tail_tokens[order(tail_tokens$col1), , drop = FALSE]
    op <- tail_tokens$token[tail_tokens$token %in% c("NS_GET", "NS_GET_INT")][1]
    fn <- tail_tokens$text[tail_tokens$token %in% c("SYMBOL", "SYMBOL_FUNCTION_CALL")][1]
    data.frame(
      file     = basename(path),
      line     = pd$line1[i],
      package  = pd$text[i],
      symbol   = if (is.na(fn)) NA_character_ else fn,
      internal = identical(op, "NS_GET_INT"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# --------------------------------------------------------------------- domain

deps_domain_of <- function(name) {
  # `.anno_x` is an anno-domain private helper, not an unprefixed name
  pref <- sub("_.*$", "", sub("^[.]+", "", name))
  ifelse(pref %in% DEPS_DOMAINS, pref, "unprefixed")
}

# ------------------------------------------------------------------- reporting

deps_pct <- function(x, total) if (total == 0) 0 else round(100 * x / total, 1)

# Packages named inside requireNamespace() / loadNamespace() / library() /
# require(). Optional backends are wired this way throughout the package, so a
# `pkg::`-only scan badly understates usage. Driven by the parser, not by a
# regex over lines: commented-out code must not count as usage, and a variable
# passed to requireNamespace() is not a package name.
DEPS_LOADER_CALLS <- c("requireNamespace", "loadNamespace", "library",
                       "require", "attachNamespace")
DEPS_LOADER_SYMBOL_OK <- c("library", "require")

deps_loader_calls <- function(path) {
  pd <- tryCatch(utils::getParseData(parse(path, keep.source = TRUE)),
                 error = function(e) NULL)
  if (is.null(pd) || !nrow(pd)) return(NULL)
  pd <- pd[pd$terminal, , drop = FALSE]
  pd <- pd[order(pd$line1, pd$col1), , drop = FALSE]
  idx <- which(pd$token == "SYMBOL_FUNCTION_CALL" & pd$text %in% DEPS_LOADER_CALLS)
  if (!length(idx)) return(NULL)
  rows <- lapply(idx, function(i) {
    nxt <- pd[seq.int(i + 1L, min(i + 3L, nrow(pd))), , drop = FALSE]
    nxt <- nxt[nxt$token %in% c("STR_CONST", "SYMBOL"), , drop = FALSE]
    if (!nrow(nxt)) return(NULL)
    tok <- nxt[1, ]
    quoted <- tok$token == "STR_CONST"
    if (!quoted && !pd$text[i] %in% DEPS_LOADER_SYMBOL_OK) return(NULL)
    data.frame(file = basename(path), line = pd$line1[i],
               package = gsub("^[\"']|[\"']$", "", tok$text),
               symbol = NA_character_, internal = FALSE,
               stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# Bare-symbol call sites of names brought in by importFrom(), attributed to a
# line so they can be mapped onto the enclosing function like any other call.
deps_import_symbol_calls <- function(path, symbol_to_pkg) {
  if (!length(symbol_to_pkg)) return(NULL)
  pd <- tryCatch(utils::getParseData(parse(path, keep.source = TRUE)),
                 error = function(e) NULL)
  if (is.null(pd) || !nrow(pd)) return(NULL)
  keep <- pd$terminal &
    pd$token %in% c("SYMBOL", "SYMBOL_FUNCTION_CALL", "SPECIAL") &
    pd$text %in% names(symbol_to_pkg)
  if (!any(keep)) return(NULL)
  hit <- pd[keep, , drop = FALSE]
  data.frame(file = basename(path), line = hit$line1,
             package = unname(symbol_to_pkg[hit$text]),
             symbol = hit$text, internal = FALSE,
             stringsAsFactors = FALSE)
}

# String constants naming a function of this package. The parallel workers are
# handed their environment through `.export = c("fn", ...)` lists, and the
# dispatchers self-qualify as `SEMseeker:::fn`. Both are real calls that no
# AST walk over the enclosing function body can see.
deps_string_function_refs <- function(path, own_funs) {
  if (!length(own_funs)) return(NULL)
  pd <- tryCatch(utils::getParseData(parse(path, keep.source = TRUE)),
                 error = function(e) NULL)
  if (is.null(pd) || !nrow(pd)) return(NULL)
  hit <- pd[pd$terminal & pd$token == "STR_CONST", , drop = FALSE]
  if (!nrow(hit)) return(NULL)
  val <- gsub('^["\']|["\']$', "", hit$text)
  keep <- val %in% own_funs
  if (!any(keep)) return(NULL)
  data.frame(file = basename(path), line = hit$line1[keep], callee = val[keep],
             stringsAsFactors = FALSE)
}

# Package names that appear only as string constants: the annotation layer
# dispatches on lookup tables of package names and then calls
# requireNamespace(pkg) on the variable, which no call-site scan can resolve.
# Restricted to names that are actually declared, so ordinary strings do not
# turn into phantom dependencies.
deps_string_references <- function(path, known_packages) {
  if (!length(known_packages)) return(NULL)
  pd <- tryCatch(utils::getParseData(parse(path, keep.source = TRUE)),
                 error = function(e) NULL)
  if (is.null(pd) || !nrow(pd)) return(NULL)
  hit <- pd[pd$terminal & pd$token == "STR_CONST", , drop = FALSE]
  if (!nrow(hit)) return(NULL)
  val <- gsub('^["\']|["\']$', "", hit$text)
  keep <- val %in% known_packages
  if (!any(keep)) return(NULL)
  data.frame(file = basename(path), line = hit$line1[keep],
             package = val[keep], symbol = NA_character_, internal = FALSE,
             stringsAsFactors = FALSE)
}
