# semseeker — CI/CD Setup

## Overview

GitHub Actions CI runs `R CMD check --as-cran` on every push to any branch
except `main`/`master`, and on every PR targeting `main`, `master`, or `develop`.
Direct pushes to `main` are blocked by GitHub branch protection; changes must
go through a PR and all CI checks must pass before merge.

Workflow file: `.github/workflows/R-CMD-check.yml`

---

## Platform matrix

| Platform | R version | Status |
|---|---|---|
| `macos-latest` | release | ✅ active |
| `windows-latest` | release | ✅ active |
| `ubuntu-latest` | release | ✅ active |

`fail-fast: false` — all platforms always run to completion even if one fails.

Windows and Ubuntu were temporarily suspended while the polars 0.22→0.23
migration and `R CMD check` compliance fixes were being resolved. They were
re-enabled once macOS CI was confirmed green (run `23841500938`, 18m42s).

---

## Branch protection

`main` is protected via GitHub API:

- Direct pushes blocked
- A PR is required
- All active CI checks must pass before merge (`strict: true`)
- Configured with:
  ```bash
  gh api --method PUT repos/corsaro-lab/SEMseeker/branches/main/protection --input - <<'EOF'
  {
    "required_status_checks": {
      "strict": true,
      "contexts": [
        "R-CMD-check / macos-latest (release)",
        "R-CMD-check / windows-latest (release)",
        "R-CMD-check / ubuntu-latest (release)"
      ]
    },
    "enforce_admins": false,
    "required_pull_request_reviews": { "required_approving_review_count": 0 },
    "restrictions": null
  }
  EOF
  ```

---

## Non-CRAN dependencies

### polars

`polars` is not on CRAN. It is distributed via the **R-multiverse** community repository:

```r
Sys.setenv(NOT_CRAN = "true")
install.packages("polars", repos = "https://community.r-multiverse.org")
```

**Why `NOT_CRAN = "true"` is required:**
The polars package performs a runtime check and refuses to install in CRAN-like
sandboxed environments unless `NOT_CRAN` is explicitly set to `"true"`.

**DESCRIPTION configuration:**

```
Additional_repositories: https://community.r-multiverse.org
```

**Why `setup-r-dependencies@v2` (pak) cannot be used for polars:**
`pak` resolves ALL `Imports` from DESCRIPTION before configuring extra
repositories. This means polars is looked up before `Additional_repositories`
is configured, producing *"Can't find package called polars"* on every run.
Neither pre-install steps nor `.Rprofile` overrides help — pak spawns R
subprocesses with `--no-save --no-restore`, bypassing `.Rprofile`.

**Solution — use `remotes::install_deps()` instead:**
`remotes` reads `Additional_repositories` *before* resolving, so polars is
found correctly:

```yaml
- name: Install R dependencies
  run: |
    Sys.setenv(NOT_CRAN = "true")
    install.packages(
      c("remotes", "rcmdcheck", "testthat", "devtools"),
      repos = "https://packagemanager.posit.co/cran/latest"
    )
    remotes::install_deps(
      dependencies = NA,   # Imports/Depends/LinkingTo only — skip Suggests
      repos = c(
        rpolars = "https://community.r-multiverse.org",
        CRAN    = "https://packagemanager.posit.co/cran/latest"
      )
    )
  shell: Rscript {0}
```

`NOT_CRAN: "true"` is set as a global env var in the workflow job so it is
active for all R processes (including `R CMD check`).

**Why `dependencies = NA` (not `TRUE`):**
`dependencies = TRUE` installs Suggests too. semseeker's `Suggests` includes
`pathfindR`, which hard-requires `ggkegg` — a Bioconductor package that is not
available on CRAN or community.r-multiverse.org. This causes
`remotes::install_deps()` to emit *"dependency 'ggkegg' is not available"* and
the subsequent test run to fail.

Using `dependencies = NA` skips Suggests entirely. All code paths that need
a Suggests package guard themselves with `requireNamespace()` and the
corresponding tests call `testthat::skip()` when the package is absent, so CI
is not affected. See also _Known Issue #5_ below.

### polars API versioning (breaking changes 0.22 → 0.23)

CI installs the **latest** polars from community.r-multiverse.org; local
development typically uses an older pinned version (0.22.4 in renv.lock).
This mismatch has caused recurring CI failures. Known breaking changes:

| Old API (≤ 0.22.4) | New API (≥ 0.23) | Affected files |
|---|---|---|
| `df$to_data_frame()` | `as.data.frame(df)` | 11 R source files, 5 test files |
| `lf$cast(list(col = dtype))` | `lf$with_columns(pl$col("c")$cast(dtype))` | `annotate_position_pivots.R` |
| `lf$rename(list(old = "new"))` | `lf$select(pl$col("old")$alias("new"))` | `create_position_pivots.R` |
| `$sort(descending = c(F, F, ...))` | `$sort(descending = FALSE)` scalar only | `signal_save.R`, `signal_range_values.R`, `deltaX_get.R`, `annotate_position_pivots.R`, `create_position_pivots.R` |

**Rule:** never chain `$to_data_frame()` on a polars object; always use
`as.data.frame()`. When casting multiple columns, use `$with_columns()` with
per-column `$cast()`. This is forward-compatible with all polars R versions.

**Long-term fix:** pin a minimum polars version in DESCRIPTION once one is
confirmed to include all fixes, or switch to `pak::lockfile` to pin the exact
version in CI.

### ctdR

`ctdR` is a private/internal package hosted on GitHub at `drake69/ctdR`.
It is installed automatically via the `Remotes` field in `DESCRIPTION`:

```
Remotes: drake69/ctdR
```

`pak` / `setup-r-dependencies@v2` handles GitHub remotes natively.

---

## Known R CMD check issues resolved

| File | Issue | Fix | Version |
|---|---|---|---|
| `compare_markers.R:16` | `inference_details` typo (unused argument warning) | Renamed to `inference_detail` | 0.11.0 |
| `marker_quantization_metric.R:56,72,82` | `next` inside `%dorng%` foreach block (not a loop) | Replaced with `return(NULL)` | 0.11.0 |
| `exploratory_analysis.R` | Non-ASCII characters in R source | Known — fix with `\uXXXX` escapes | open |
| Many R files | `$to_data_frame()` removed in polars 0.23 | Replaced with `as.data.frame()` | 0.11.x |
| `annotate_position_pivots.R` | `$cast(list(...))` API changed in polars 0.23 | Replaced with `$with_columns(...$cast(...))` | 0.11.x |
| `create_position_pivots.R` | `$rename(list(...))` API changed | Replaced with `$select(...$alias(...))` | 0.11.x |
| Multiple test files | Internal functions called without `semseeker:::` — works under `devtools::load_all()` but fails in `R CMD check` installed-package mode | Added `semseeker:::` prefix | 0.11.x |
| `test-0-init_env.R` | Path assertion `== paste0(folder,"/Data")` fails on macOS when `normalizePath()` resolves `/var` symlink to `/private/var` | Changed to `normalizePath(file.path(...), mustWork=FALSE)` | 0.11.x |
| `test-coverage.yml` | `dependencies = TRUE` triggers `pathfindR` install → requires `ggkegg` (Bioconductor) → not available → test run crashes | Changed to `dependencies = NA`; pathway tests skip via `requireNamespace()` guard | 0.11.x |

---

## Disabled / legacy workflows

Two legacy workflow files were disabled to avoid duplicate runs and parse errors:

| File | Reason |
|---|---|
| `all-actions.yaml` | Was entirely commented out — no `on:` key, GitHub rejected it instantly |
| `cmd_check.yaml` | Duplicate of `R-CMD-check.yml`, macOS-only, less strict |

Both now contain only a `workflow_dispatch` trigger with a no-op placeholder job.
They can be safely deleted once the history is no longer needed.

---

## testthat configuration

`Config/testthat/parallel: false` is set in `DESCRIPTION`.

**Why:** testthat parallel mode spawns worker processes that do **not** share the
global environment initialised by `tests/testthat/setup.R`. Variables set with
`<<-` (e.g. `signal_data`, `signal_thresholds`, `mySampleSheet`) are invisible
in parallel workers, causing *"object not found"* warnings that `R CMD check`
treats as test failures. Sequential execution avoids this entirely.

If parallelism is re-enabled in the future, `setup.R` must be refactored to use
`withr::local_options()` or testthat fixtures so each worker initialises its own
state independently.

---

## Local test execution

### Option A — devtools (fastest, no container)

Bypass renv to run tests locally against the system library:

```bash
RENV_CONFIG_SANDBOX_ENABLED=FALSE RENV_ACTIVATE_PROJECT=0 \
  Rscript --vanilla -e 'devtools::test()'
```

Or per-file:

```bash
RENV_CONFIG_SANDBOX_ENABLED=FALSE RENV_ACTIVATE_PROJECT=0 \
  Rscript --vanilla -e '
    devtools::load_all(".")
    source("tests/testthat/setup.R")
    testthat::test_file("tests/testthat/test-6-semseeker.R")
  '
```

**Important caveat:** `devtools::load_all()` puts ALL internal functions into
the global environment, so calls like `pivot_file_name_parquet()` work without
`semseeker:::`. This masks a class of bugs that only appear under `R CMD check`
(installed package mode). To catch these, always run Option B before pushing.

### Option B — local Docker container (matches CI exactly)

`Dockerfile.ci` + `ci-local.sh` provide a container that mirrors GitHub Actions:
same R version (4.5), same polars from community.r-multiverse.org, same
`R CMD check --as-cran` flags.

```bash
# First run: builds image + installs all R packages (~5–10 min, layers cached)
./ci-local.sh

# Subsequent runs: reinstalls package + runs check (~1–2 min)
./ci-local.sh
```

**Error collection:**

| Method | Command | What you get |
|---|---|---|
| Stdout only | `./ci-local.sh 2>&1 \| tee ci-local.log` | Full terminal output saved to `ci-local.log` |
| Full artefacts | `./ci-local.sh logs` | `./ci-check-output/semseeker.Rcheck/` directory |
| Interactive debug | `./ci-local.sh shell` | Bash shell inside container at `/pkg` |

The `.Rcheck` directory contains:

| File | Contents |
|---|---|
| `00check.log` | Complete `R CMD check` output |
| `00install.out` | Package installation transcript |
| `semseeker/tests/testthat.Rout` | All `test_that` output, including passed tests |
| `semseeker-Ex.Rout` | Output from running `\examples{}` in man pages |

**Difference vs GitHub CI:** container runs on Linux (Debian), not macOS. This
catches polars API issues, `semseeker:::` visibility bugs, and `R CMD check`
compliance. It will not reproduce macOS-specific issues such as symlink
resolution differences (`/var` → `/private/var`), which are handled separately
in the test assertions.

---

## Codecov coverage

A second workflow `.github/workflows/test-coverage.yml` runs `covr::package_coverage()`
on macOS and uploads results to Codecov via `codecov/codecov-action@v4`.

`CODECOV_TOKEN` must be set as a GitHub Actions secret (already configured).
`covr` installs the package fresh in an isolated environment, so the same
`semseeker:::` visibility rules apply as in `R CMD check`.

The coverage workflow also uses `dependencies = NA` (same rationale as
`R-CMD-check.yml` above — avoid the `pathfindR` → `ggkegg` Bioconductor chain).

---

## Adding a Suggests package with Bioconductor transitive deps

If a new optional feature requires a CRAN package that itself depends on
Bioconductor packages (e.g. `pathfindR` → `ggkegg`), follow this checklist:

1. Add the CRAN package to `Suggests:` in `DESCRIPTION` only — **never**
   `Imports:` or `Depends:`.
2. Guard every call site in R source:
   ```r
   if (!requireNamespace("mypkg", quietly = TRUE)) {
     message("Install mypkg to use this feature: install.packages('mypkg')")
     return(invisible(NULL))
   }
   mypkg::some_function(...)
   ```
3. Guard the test:
   ```r
   test_that("feature works when mypkg is available", {
     if (!requireNamespace("mypkg", quietly = TRUE))
       testthat::skip("mypkg not installed")
     # ... assertions
   })
   ```
4. Do **not** add Bioconductor repositories to the CI `repos` vector. CI
   intentionally skips Suggests; the test will be skipped, not failed.
5. Document the optional dependency in `README.md` under an *Optional features*
   section so users know how to install it manually.
