# 03_join_minimal_core.R: where the two graphs meet.
#
# Every external dependency is kept alive by a set of internal functions. Cross
# that set with the reachability computed in 02 and each dependency lands in one
# of five buckets:
#
#   blocking   not installable from CRAN or Bioconductor: submission stopper
#   core       reached by all three endpoints, or by the SEM layer alone
#   layer      reached by one downstream layer only (assoc or enrich)
#   periphery  used only by code no exported function can reach
#   unused     no call site in R/ at all
#
# The counterfactual install closure ("what if only core stayed in Imports")
# is what makes "how much is improvable" a number rather than an opinion.

deps_classify <- function(external, internal) {
  ext <- external$deps
  fp  <- internal$fun_pkg
  funs <- internal$functions
  reach <- setNames(funs$reach_exported, funs$fun)
  ep_cols <- paste0("reach_", ENDPOINTS)

  layer_of <- function(pkg) {
    users <- fp$fun[fp$package == pkg]
    if (!length(users)) return(NA_character_)
    rows <- funs[funs$fun %in% users, , drop = FALSE]
    hit <- vapply(ep_cols, function(cc) any(rows[[cc]]), logical(1))
    if (!any(hit)) return("none")
    paste(sub("^reach_", "", ep_cols[hit]), collapse = "+")
  }

  ext$user_functions <- vapply(ext$package, function(p)
    length(unique(fp$fun[fp$package == p])), integer(1))
  ext$reachable_users <- vapply(ext$package, function(p) {
    u <- unique(fp$fun[fp$package == p])
    sum(reach[u], na.rm = TRUE)
  }, integer(1))
  ext$endpoints <- vapply(ext$package, layer_of, character(1))

  ext$bucket <- with(ext, ifelse(
    field %in% c("Depends", "Imports") &
      !origin %in% c("CRAN", "Bioconductor", "base"), "blocking",
    ifelse(!used_in_R & field %in% c("Depends", "Imports"), "unused",
      ifelse(reachable_users == 0 & used_in_R, "periphery",
        ifelse(!is.na(endpoints) & grepl("[+]", endpoints), "core",
          ifelse(!is.na(endpoints) & endpoints != "none", "layer", "periphery"))))))
  ext$bucket[ext$field == "Suggests" & !ext$used_anywhere] <- "unused"
  ext$bucket[ext$field == "Suggests" & ext$used_anywhere] <- "suggests"
  ext
}

# What the install closure would be if only the packages in `keep` stayed hard
# dependencies. This is the improvement headroom, in packages.
deps_counterfactual <- function(keep, paths) {
  db <- deps_available_db(paths)
  cl <- unique(unlist(lapply(keep, deps_closure, db = db)))
  length(union(keep, cl))
}

run_join_minimal_core <- function(pkg_root = ".", external, internal) {
  paths <- deps_paths(pkg_root)
  deps_msg("[03] joining the two graphs")

  ext <- deps_classify(external, internal)
  deps_write_csv(ext, file.path(paths$out, "dependency_buckets.csv"))

  hard <- ext$package[ext$field %in% c("Depends", "Imports")]
  keep_core  <- ext$package[ext$bucket == "core"]
  keep_layer <- ext$package[ext$bucket == "layer"]

  actual        <- deps_counterfactual(hard, paths)
  core_only     <- deps_counterfactual(keep_core, paths)
  core_plus_sem <- deps_counterfactual(union(keep_core, keep_layer), paths)

  funs <- internal$functions
  minimal_core_funs <- funs$fun[funs$reach_endpoints]

  # Per-endpoint cost: the code each layer reaches and the install closure of
  # the Imports that code actually touches. This is the table the minimal
  # end-to-end target is built on.
  fp <- internal$fun_pkg
  # Strong evidence only here: a package whose name merely appears in a string
  # inside a reachable function is not thereby a dependency of that function.
  fp_strong <- fp[fp$kind %in% c("ns", "import", "loader"), , drop = FALSE]
  imports_for <- function(fun_set) {
    intersect(unique(fp_strong$package[fp_strong$fun %in% fun_set]),
              ext$package[ext$field %in% c("Depends", "Imports")])
  }
  per_layer <- lapply(c(ENDPOINTS, "UNION"), function(ep) {
    sel <- if (ep == "UNION") funs$reach_endpoints else funs[[paste0("reach_", ep)]]
    imps <- imports_for(funs$fun[sel])
    list(endpoint = ep,
         functions = sum(sel),
         files     = length(unique(funs$file[sel])),
         loc       = sum(funs$n_lines[sel], na.rm = TRUE),
         imports   = length(imps),
         closure   = deps_counterfactual(imps, paths))
  })
  names(per_layer) <- c(ENDPOINTS, "UNION")
  reachable_imports <- imports_for(minimal_core_funs)

  summary <- list(
    buckets = as.list(table(ext$bucket)),
    blocking  = sort(ext$package[ext$bucket == "blocking"]),
    unused    = sort(ext$package[ext$bucket == "unused"]),
    periphery = sort(setdiff(ext$package[ext$bucket == "periphery"], DEPS_BASE_PKGS)),
    layer_only = setNames(
      as.list(ext$endpoints[ext$bucket == "layer"]),
      ext$package[ext$bucket == "layer"]
    ),
    per_layer                = unname(per_layer),
    imports_never_reached    = sort(setdiff(
      ext$package[ext$field %in% c("Depends", "Imports")], reachable_imports)),
    closure_now              = actual,
    closure_endpoint_reachable = deps_counterfactual(reachable_imports, paths),
    closure_core_only        = core_only,
    closure_core_plus_layer  = core_plus_sem,
    packages_removable       = actual - deps_counterfactual(reachable_imports, paths),
    minimal_core_functions   = length(minimal_core_funs),
    minimal_core_files       = length(unique(funs$file[funs$reach_endpoints])),
    total_functions          = nrow(funs)
  )
  jsonlite::write_json(summary, file.path(paths$out, "join_summary.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  deps_msg("  install closure now %d → %d if only endpoint-reachable Imports stayed hard (−%d)",
           actual, summary$closure_endpoint_reachable,
           actual - summary$closure_endpoint_reachable)
  invisible(list(deps = ext, summary = summary))
}

# ------------------------------------------------------------------- rendering

# Domain-level graph as DOT. Small enough to read; the function-level graph is
# 1500+ nodes and is left as CSV for whoever wants to load it into igraph.
deps_render_domain_dot <- function(domain_matrix, path) {
  lines <- c("digraph domains {", "  rankdir=LR;",
             "  node [shape=box, style=rounded, fontname=Helvetica];",
             "  edge [fontname=Helvetica, fontsize=9];")
  for (i in seq_len(nrow(domain_matrix))) {
    lines <- c(lines, sprintf('  "%s" -> "%s" [label="%d", penwidth=%.1f];',
                              domain_matrix$from[i], domain_matrix$to[i],
                              domain_matrix$weight[i],
                              max(0.5, min(6, log10(domain_matrix$weight[i] + 1) * 2))))
  }
  writeLines(c(lines, "}"), path)
  deps_msg("  wrote %s", basename(path))
}

# Mermaid version of the same graph, for pasting into the analysis document.
# graphviz is not installed everywhere and mermaid renders in Markdown.
deps_render_domain_mermaid <- function(domain_matrix, feedback_arcs, path,
                                       min_weight = 3) {
  fb <- vapply(feedback_arcs, function(a) paste(a$from, a$to), character(1))
  d <- domain_matrix[domain_matrix$from != domain_matrix$to &
                       domain_matrix$weight >= min_weight, , drop = FALSE]
  d <- d[order(-d$weight), , drop = FALSE]
  lines <- c("flowchart LR")
  for (i in seq_len(nrow(d))) {
    back <- paste(d$from[i], d$to[i]) %in% fb
    arrow <- if (back) "-.->" else "-->"
    lines <- c(lines, sprintf("  %s %s|%d| %s", d$from[i], arrow, d$weight[i], d$to[i]))
  }
  lines <- c(lines, sprintf("  %%%% dashed = feedback arc (breaks the layering); edges below %d calls omitted",
                            min_weight))
  writeLines(lines, path)
  deps_msg("  wrote %s", basename(path))
}
