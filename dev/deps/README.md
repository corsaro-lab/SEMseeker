# dev/deps: dependency analysis pipeline

Measures two graphs and joins them:

1. **External**: what SEMseeker declares, where each package comes from
   (CRAN / Bioconductor / elsewhere), how large an install each one drags in,
   and whether it is used at all.
2. **Internal**: which function calls which, lifted to file and domain level,
   plus reachability from the three public endpoints.

The join answers the only question that matters for a submission: *which
declared dependency is kept alive by code a user can actually reach.*

## Run

```sh
Rscript dev/deps/run_dependency_analysis.R            # uses the cached repo metadata
Rscript dev/deps/run_dependency_analysis.R --refresh  # re-fetch CRAN + Bioc + r-multiverse
```

Needs `igraph`, `jsonlite`, `codetools`. Nothing is loaded or attached from
SEMseeker itself: the package is read, never run. Takes ~2 minutes, most of it
in `tools::package_dependencies()` over the full CRAN index.

## Output

`dev/deps/output/` (gitignored):

| file | contents |
|---|---|
| `external_deps.csv` | one row per declared dependency: origin, closure size, call sites, evidence kind |
| `external_closure.csv` | the full recursive install closure of Depends+Imports |
| `external_call_sites.csv` | every call site into another package |
| `call_sites_attributed.csv` | the same, with the enclosing function resolved |
| `internal_functions.csv` | every function: file, domain, size, fan-in/out, per-endpoint reachability |
| `internal_edges.csv` | the call graph, with how each edge was found |
| `domain_matrix.csv` | domain → domain call counts |
| `domains.dot` / `domains.mmd` | the domain graph, graphviz and mermaid |
| `*_summary.json` | the headline numbers |

`dev/deps/snapshots/YYYY-MM-DD_metrics.json` is committed on purpose: the point
is a trend across releases, not a single reading.

## Evidence strength

A package can be reached in four ways, and they are not equally good evidence:

| kind | what it is | strength |
|---|---|---|
| `ns` | `pkg::fn()` | proof |
| `import` | a symbol brought in by `importFrom()` | proof |
| `loader` | `requireNamespace("pkg")` with a literal | proof of an optional path |
| `string` | the name only ever appears as a string constant | **weak**: may be dynamic dispatch, may be dead |

The `evidence` column carries this. Anything resting on `string` alone needs a
human to look before it is called used or unused.

## What the analysis cannot see

Static analysis, so: dispatch through a variable, `do.call()` on a computed
name, and `get()` are invisible. Two forms that *would* have been invisible are
recovered explicitly, because this package leans on both: `SEMseeker:::fn()`
self-qualified calls, and function names handed to parallel workers as strings
in `.export =` lists. Everything unreachable is therefore a **candidate**, not a
proof; every candidate in the analysis document was checked by hand before being
reported.
