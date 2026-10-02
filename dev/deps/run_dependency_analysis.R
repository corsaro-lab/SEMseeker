#!/usr/bin/env Rscript
# run_dependency_analysis.R: orchestrator.
#
#   Rscript dev/deps/run_dependency_analysis.R [--refresh] [--root <pkg-root>]
#
# --refresh re-fetches CRAN/Bioc metadata instead of using dev/deps/cache/.
#
# Writes machine-readable artifacts to dev/deps/output/ (gitignored) and one
# small metrics snapshot to dev/deps/snapshots/ (committed), so the same numbers
# can be recomputed later and compared: the point is a trend, not a one-off
# reading.

args <- commandArgs(trailingOnly = TRUE)
refresh <- "--refresh" %in% args
root <- if ("--root" %in% args) args[which(args == "--root") + 1L] else "."

here <- file.path(root, "dev", "deps")
source(file.path(here, "lib_deps.R"))
source(file.path(here, "01_external_deps.R"))
source(file.path(here, "02_internal_graph.R"))
source(file.path(here, "03_join_minimal_core.R"))

for (p in c("igraph", "jsonlite", "codetools")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("missing package: ", p)
}

paths <- deps_paths(root)
dir.create(paths$out, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$cache, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$snapshots, showWarnings = FALSE, recursive = TRUE)

t0 <- Sys.time()
external <- run_external_deps(root, refresh = refresh)
internal <- run_internal_graph(root, external = external)
joined   <- run_join_minimal_core(root, external, internal)

deps_render_domain_dot(internal$domain_matrix, file.path(paths$out, "domains.dot"))
deps_render_domain_mermaid(internal$domain_matrix, internal$summary$domain_feedback_arcs,
                           file.path(paths$out, "domains.mmd"))

desc <- read.dcf(paths$desc)[1, ]
snapshot <- list(
  package  = unname(desc[["Package"]]),
  version  = unname(desc[["Version"]]),
  date     = format(Sys.Date()),
  external = external$summary,
  internal = internal$summary,
  join     = joined$summary
)
# Never overwrite a snapshot: two runs on the same day are exactly the
# before/after of a change, and the "before" is the half worth keeping.
snap_path <- file.path(paths$snapshots,
                       sprintf("%s_metrics.json", format(Sys.Date())))
if (file.exists(snap_path)) {
  n <- 2L
  repeat {
    snap_path <- file.path(paths$snapshots,
                           sprintf("%s_metrics_%d.json", format(Sys.Date()), n))
    if (!file.exists(snap_path)) break
    n <- n + 1L
  }
}
jsonlite::write_json(snapshot, snap_path, auto_unbox = TRUE, pretty = TRUE)

deps_msg("")
deps_msg("snapshot: %s", snap_path)
deps_msg("elapsed:  %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs")))
