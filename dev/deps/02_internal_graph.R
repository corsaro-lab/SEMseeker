# 02_internal_graph.R: the internal call graph.
#
# Builds a function-level directed graph of R/ (who calls whom), lifts it to the
# file and domain levels, and computes reachability from the three public
# endpoints. Reachability is what the minimal end-to-end surface is derived
# from: everything the endpoints cannot reach is, by construction, outside it.
#
# Static analysis only. Dispatch through variables, do.call() with a computed
# name, and get() are invisible to it: the reachable set is therefore a lower
# bound and the unreachable set must be read as "candidates", not as proof.

ENDPOINTS <- c("semseeker", "association_analysis", "enrichment_analysis")

# Called by R itself, never by package code: roots of the graph even though
# nothing points at them.
NAMESPACE_HOOKS <- c(".onLoad", ".onUnload", ".onAttach", ".onDetach", ".Last.lib")

deps_collect_functions <- function(r_dir) {
  files <- list.files(r_dir, pattern = "[.][Rr]$", full.names = TRUE)
  tab <- do.call(rbind, lapply(files, deps_file_functions))
  tab$domain <- deps_domain_of(tab$fun)
  tab$file_domain <- deps_domain_of(tab$file)
  tab
}

# Function-to-function edges, restricted to functions defined in the package.
deps_call_edges <- function(fun_env, own_funs) {
  rows <- list()
  for (nm in ls(fun_env, all.names = TRUE)) {
    f <- get(nm, envir = fun_env)
    if (!is.function(f)) next
    globals <- tryCatch(codetools::findGlobals(f, merge = FALSE)$functions,
                        error = function(e) character())
    callees <- intersect(globals, own_funs)
    callees <- setdiff(callees, nm)
    if (!length(callees)) next
    rows[[nm]] <- data.frame(from = nm, to = callees, stringsAsFactors = FALSE)
  }
  if (!length(rows)) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  unique(do.call(rbind, rows))
}

# Map each external call site (file, line) onto the function whose span contains it.
deps_attribute_call_sites <- function(call_sites, funs) {
  if (!nrow(call_sites)) {
    return(cbind(call_sites, fun = character()))
  }
  owner <- vapply(seq_len(nrow(call_sites)), function(i) {
    cand <- funs[funs$file == call_sites$file[i] &
                   funs$line_from <= call_sites$line[i] &
                   funs$line_to   >= call_sites$line[i], , drop = FALSE]
    if (!nrow(cand)) return(NA_character_)
    cand$fun[which.min(cand$n_lines)]  # innermost span wins
  }, character(1))
  call_sites$fun <- owner
  call_sites
}

run_internal_graph <- function(pkg_root = ".", external = NULL) {
  paths <- deps_paths(pkg_root)
  deps_msg("[02] internal call graph")

  funs <- deps_collect_functions(paths$r_dir)
  imports <- deps_namespace_imports(paths$namespace)
  fun_env <- deps_function_env(paths$r_dir)
  failed <- attr(fun_env, "failed")
  if (length(failed)) deps_msg("  WARNING: %d file(s) not fully materialised: %s",
                               length(failed), paste(failed, collapse = ", "))

  own <- unique(funs$fun)
  edges <- deps_call_edges(fun_env, own)
  edges$kind <- "ast"

  # Two call forms the AST walk cannot see, both of them real:
  #   SEMseeker:::fn(...)   self-qualified dispatch inside %dofuture% blocks
  #   .export = c("fn")     function names handed to workers as strings
  # Without them the reachable set is wrong, not merely conservative.
  extra <- list()
  for (path in list.files(paths$r_dir, pattern = "[.][Rr]$", full.names = TRUE)) {
    self <- deps_qualified_calls(path)
    if (!is.null(self)) {
      self <- self[self$package %in% c("SEMseeker", "semseeker") &
                     !is.na(self$symbol) & self$symbol %in% own, , drop = FALSE]
      if (nrow(self)) extra[[length(extra) + 1L]] <-
        data.frame(file = self$file, line = self$line, callee = self$symbol,
                   kind = "self_ns", stringsAsFactors = FALSE)
    }
    st <- deps_string_function_refs(path, own)
    if (!is.null(st)) extra[[length(extra) + 1L]] <-
      data.frame(file = st$file, line = st$line, callee = st$callee,
                 kind = "string", stringsAsFactors = FALSE)
  }
  if (length(extra)) {
    extra <- do.call(rbind, extra)
    caller <- vapply(seq_len(nrow(extra)), function(i) {
      cand <- funs[funs$file == extra$file[i] &
                     funs$line_from <= extra$line[i] &
                     funs$line_to   >= extra$line[i], , drop = FALSE]
      if (!nrow(cand)) return(NA_character_)
      cand$fun[which.min(cand$n_lines)]
    }, character(1))
    extra$from <- caller
    extra <- extra[!is.na(extra$from) & extra$from != extra$callee, , drop = FALSE]
    if (nrow(extra)) {
      add <- unique(data.frame(from = extra$from, to = extra$callee,
                               kind = extra$kind, stringsAsFactors = FALSE))
      add <- add[!paste(add$from, add$to) %in% paste(edges$from, edges$to), ,
                 drop = FALSE]
      edges <- rbind(edges, add)
    }
  }
  funs$exported <- funs$fun %in% imports$exported

  # a name defined in two files is one vertex; the duplication is itself a
  # finding, so it is reported rather than silently merged
  duplicated_names <- sort(unique(funs$fun[duplicated(funs$fun)]))
  g <- igraph::graph_from_data_frame(edges[c("from", "to")], directed = TRUE,
                                     vertices = own)
  funs$fan_out <- igraph::degree(g, mode = "out")[funs$fun]
  funs$fan_in  <- igraph::degree(g, mode = "in")[funs$fun]

  # ---- reachability from the endpoints and from the whole export surface
  reach_from <- function(roots) {
    roots <- intersect(roots, igraph::V(g)$name)
    if (!length(roots)) return(character())
    r <- igraph::subcomponent(g, roots[1], mode = "out")
    for (rt in roots[-1]) r <- igraph::union(r, igraph::subcomponent(g, rt, mode = "out"))
    unique(igraph::V(g)$name[as.integer(r)])
  }
  reach_endpoints <- reach_from(ENDPOINTS)
  reach_exported  <- reach_from(c(imports$exported, NAMESPACE_HOOKS))

  funs$reach_endpoints <- funs$fun %in% reach_endpoints
  funs$reach_exported  <- funs$fun %in% reach_exported

  # per-endpoint reachability: which layer keeps a function alive
  for (ep in ENDPOINTS) {
    funs[[paste0("reach_", ep)]] <- funs$fun %in% reach_from(ep)
  }

  # ---- domain-level graph and cycles
  dom_of <- setNames(funs$domain, funs$fun)
  dom_edges <- data.frame(
    from = unname(dom_of[edges$from]),
    to   = unname(dom_of[edges$to]),
    stringsAsFactors = FALSE
  )
  dom_edges <- dom_edges[!is.na(dom_edges$from) & !is.na(dom_edges$to), ]
  dom_w <- aggregate(list(weight = rep(1, nrow(dom_edges))),
                     dom_edges[c("from", "to")], sum)
  dom_g <- igraph::graph_from_data_frame(dom_w[dom_w$from != dom_w$to, ], directed = TRUE)
  dom_scc <- igraph::components(dom_g, mode = "strong")
  dom_cycles <- split(names(dom_scc$membership), dom_scc$membership)
  dom_cycles <- dom_cycles[lengths(dom_cycles) > 1]

  # How far is the domain graph from being a DAG? The minimum feedback arc set,
  # weighted by how many function calls ride on each domain edge, is the honest
  # answer: cut these edges and the domains layer cleanly.
  fas_idx <- igraph::feedback_arc_set(dom_g, weights = igraph::E(dom_g)$weight,
                                      algo = "exact_ip")
  fas <- igraph::as_data_frame(igraph::subgraph_from_edges(dom_g, fas_idx,
                                                           delete.vertices = FALSE))
  fas <- fas[order(-fas$weight), , drop = FALSE]

  # ---- file-level graph and its strongly connected components
  file_of <- setNames(funs$file, funs$fun)
  file_edges <- unique(data.frame(
    from = unname(file_of[edges$from]),
    to   = unname(file_of[edges$to]),
    stringsAsFactors = FALSE
  ))
  file_edges <- file_edges[file_edges$from != file_edges$to, ]
  file_g <- igraph::graph_from_data_frame(file_edges, directed = TRUE,
                                          vertices = unique(funs$file))
  file_scc <- igraph::components(file_g, mode = "strong")
  file_cycles <- split(names(file_scc$membership), file_scc$membership)
  file_cycles <- file_cycles[lengths(file_cycles) > 1]

  # ---- external packages per function (feeds 03)
  fun_pkg <- NULL
  if (!is.null(external)) {
    cs <- deps_attribute_call_sites(external$usage$qualified, funs)
    fun_pkg <- unique(cs[!is.na(cs$fun), c("fun", "package", "kind")])
    deps_write_csv(cs, file.path(paths$out, "call_sites_attributed.csv"))
  }

  deps_write_csv(funs, file.path(paths$out, "internal_functions.csv"))
  deps_write_csv(edges, file.path(paths$out, "internal_edges.csv"))
  deps_write_csv(dom_w, file.path(paths$out, "domain_matrix.csv"))

  summary <- list(
    n_files      = length(unique(funs$file)),
    n_functions  = nrow(funs),
    n_exported   = sum(funs$exported),
    n_edges      = nrow(edges),
    edges_by_kind = as.list(table(edges$kind)),
    median_fun_lines = stats::median(funs$n_lines, na.rm = TRUE),
    functions_over_200_lines = sort(funs$fun[funs$n_lines > 200]),
    reachable_from_endpoints = length(reach_endpoints),
    reachable_from_exports   = length(reach_exported),
    unreachable_from_exports = sort(funs$fun[!funs$reach_exported]),
    domain_cycles = unname(lapply(dom_cycles, sort)),
    domain_cross_edges = sum(dom_w$weight[dom_w$from != dom_w$to]),
    domain_feedback_arcs = if (nrow(fas))
      Map(function(f, t, w) list(from = f, to = t, calls = w),
          fas$from, fas$to, fas$weight) else list(),
    domain_feedback_calls = sum(fas$weight),
    n_file_cycles = length(file_cycles),
    largest_file_cycle = if (length(file_cycles))
      sort(file_cycles[[which.max(lengths(file_cycles))]]) else character(),
    top_fan_in = as.list(sort(setNames(funs$fan_in, funs$fun), decreasing = TRUE)[1:15]),
    duplicated_definitions = duplicated_names,
    parse_failures = failed
  )
  jsonlite::write_json(summary, file.path(paths$out, "internal_summary.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  deps_msg("  %d functions in %d files, %d edges; %d reachable from the 3 endpoints",
           nrow(funs), length(unique(funs$file)), nrow(edges), length(reach_endpoints))

  invisible(list(functions = funs, edges = edges, graph = g,
                 domain_matrix = dom_w, fun_pkg = fun_pkg, summary = summary))
}
