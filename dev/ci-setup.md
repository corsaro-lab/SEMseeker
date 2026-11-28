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

`main` and `develop` are protected with the same settings (read back from the
GitHub API on 2026-10-06):

- direct pushes are blocked: every change arrives through a pull request;
- no approving review is required (`required_approving_review_count: 0`), but
  the rules apply to administrators too (`enforce_admins: true`);
- `strict: false`: a branch does not have to be up to date with its base
  before merging;
- four required checks. The names carry no workflow prefix:
  - `macos-latest (release)`, `windows-latest (release)`,
    `ubuntu-latest (release)` from `R-CMD-check.yml`;
  - `internal-references` from `security.yml`.

The repository allows squash, merge commits and auto-merge. Feature pull
requests into `develop` are squashed; the release pull request from `develop`
into `main` is merged with a merge commit, so that the two branches keep a
common ancestor and the next release does not conflict.

To add a required check without touching the others:

```bash
gh api --method POST \
  repos/corsaro-lab/SEMseeker/branches/<branch>/protection/required_status_checks/contexts \
  -f 'contexts[]=<check name>'
```

---

## Internal references check

The repository is public, and references to the private project-management
tooling (backlog item identifiers, the paths of the files that hold them) are
refused everywhere in the tracked tree. One script holds the rule,
`.github/scripts/check-internal-references.sh`, and runs in two places:

- in CI, job `internal-references` of `security.yml`, on the whole tree;
- locally, as a pre-commit hook on the lines a commit adds. Enable it once per
  clone:

  ```bash
  git config core.hooksPath .githooks
  ```

`.gitignore` and `.Rbuildignore` are exempt: they name the agent working
folder precisely to keep it out.

---

## Dependency installation

Every `Imports` and `Depends` entry is on CRAN or Bioconductor: since the move
of the data engine to DuckDB no hard dependency comes from another repository,
and `DESCRIPTION` declares neither `Additional_repositories` nor `Remotes`.

The step *Install R dependencies* of `R-CMD-check.yml` runs in four calls:

1. `remotes::install_deps(dependencies = NA)` with the repositories configured
   by `setup-r@v2` (on Linux the distribution-specific P3M URL, which serves
   binaries). `dependencies = NA` installs `Imports`, `Depends` and
   `LinkingTo` only: the full `Suggests` set pulls in chains such as
   `pathfindR` -> `ggkegg` that are not needed for `R CMD check` and have
   repeatedly broken CI.
2. `BiocManager::install(c("limma", "KEGGREST"))` on its own. The Illumina
   annotation packages load `minfi` while they install, and `minfi` needs
   `limma`; `KEGGREST` is loaded by `bumphunter` and `AnnotationDbi` during the
   same load tests. BiocManager installs a vector in the order given, not in
   dependency order, so these two go first.
3. The Bioconductor packages the test suite uses: the three Illumina annotation
   packages, `AnnotationHub`, `TxDb.Hsapiens.UCSC.hg19.knownGene`,
   `org.Hs.eg.db`, `qvalue` and `ctdR` (on Bioconductor since 3.24). On Linux
   `ctdR` needs the system library `libmagick++-dev`: it imports `GSVA`, which
   imports `SpatialExperiment`, which imports `magick`, and the P3M binary of
   `magick` links `libMagick++`. Without it `ctdR` installs but cannot load,
   and the meta-test below reports it as missing.
4. The CRAN `Suggests` the test suite uses: `Rfast`, `mockery`, `ggpubr`,
   `circlize`, `nanoparquet`.

`options(Ncpus = 1L)` keeps the installs serial: on the Windows runner a
parallel install load-tested a package before its dependencies were built.

The two lists in steps 3 and 4 are enumerated on purpose and checked by
`tests/testthat/test-pkg_suggests_installed.R`, which scans every
`skip_if_not_installed("pkg")` call in the tests and fails when the named
package is missing, so a test cannot skip silently because CI forgot to install
it. `test-coverage.yml` keeps its own install list, also with
`dependencies = NA`: nothing ties the two lists together, they are aligned by
hand and have drifted before.

`pathfindR`, `WebGestaltR` and the other enrichment backends are `Suggests`
that CI does not install. The code guards each call with
`requireNamespace()`; the tests that need one of them skip.

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
| `test-core_init_env.R` | Path assertion `== paste0(folder,"/Data")` fails on macOS when `normalizePath()` resolves `/var` symlink to `/private/var` | Changed to `normalizePath(file.path(...), mustWork=FALSE)` | 0.11.x |
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
    testthat::test_file("tests/testthat/test-e2e_semseeker.R")
  '
```

**Important caveat:** `devtools::load_all()` puts ALL internal functions into
the global environment, so calls like `pivot_file_name_parquet()` work without
`semseeker:::`. This masks a class of bugs that only appear under `R CMD check`
(installed package mode). To catch these, always run Option B before pushing.

### Option B — local Docker container (matches CI exactly)

`Dockerfile.ci` + `ci-local.sh` provide a container that mirrors GitHub Actions:
same R version (4.6), the same Bioconductor annotation packages, the same
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
catches data-engine issues, `semseeker:::` visibility bugs, and `R CMD check`
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
