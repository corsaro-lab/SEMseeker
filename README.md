
<!-- README.md is generated from README.Rmd. Please edit that file -->

# semseeker

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5095416.svg)](https://doi.org/10.5281/zenodo.5095416)
[![](https://img.shields.io/badge/devel%20version-0.99.4-blue.svg)](https://github.com/corsaro-lab/SEMseeker)
[![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![](https://codecov.io/gh/corsaro-lab/SEMseeker/branch/main/graph/badge.svg)](https://codecov.io/gh/corsaro-lab/SEMseeker)
[![](https://img.shields.io/github/last-commit/corsaro-lab/SEMseeker.svg)](https://github.com/corsaro-lab/SEMseeker/commits/main)
[![R build status](https://github.com/corsaro-lab/SEMseeker/workflows/R-CMD-check/badge.svg)](https://github.com/corsaro-lab/SEMseeker/actions)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

<!-- badges: end -->

**semseeker** detects **Stochastic Epigenetic Mutations (SEMs)** — aberrant,
sample-specific methylation changes that deviate from the expected population
signal. Starting from a normalised beta-value (or M-value) matrix, semseeker
identifies hypomethylated and hypermethylated mutations at single-probe
resolution, aggregates them into lesions (genomic clusters of co-occurring
mutations), and computes per-region delta metrics — the signed deviation
beyond the reference interval (`DELTAS`), its ratio to the interval width
(`DELTAR`), and their equal-width and quantile ranked variants
(`DELTAP`/`DELTARP`, `DELTAQ`/`DELTARQ`) — that quantify stochastic epigenetic
variability per region.

Supported input technologies:

- **Illumina microarrays**: EPIC (850k), 450k, 27k — detected automatically
- **Bisulfite sequencing (WGBS / RRBS)** and **long-read sequencing (ONT / PacBio)**: planned (see [documents/requirements.md](documents/requirements.md))

---

## Installation

Install from GitHub using devtools:

```r
install.packages("devtools")
library(devtools)
install_github("corsaro-lab/SEMseeker")
```

### System requirements

#### macOS

R 4.6's base `tcltk` package is linked against XQuartz's
`/opt/X11/lib/libX11.6.dylib` at install time. Plain SEMseeker usage
(`library(SEMseeker)`, `semseeker(...)`, `association_analysis(...)`)
does NOT load tcltk and works without XQuartz. But some Bioconductor
optional features (e.g. annotation chains pulling `minfi` / `GEOquery`)
transitively load tcltk, and without XQuartz the load segfaults with
`invalid permissions`. To avoid the surprise:

```sh
brew install --cask xquartz
```

or download the installer from <https://www.xquartz.org/>. A logout/login
may be required for the path to be picked up.

#### Linux (build-from-source)

Building from source on Linux requires `libuv` ≥ 1, a build-time dependency
of the `fs` package (a transitive dep of `rmarkdown`/`pkgdown`/`bslib`).
Either install it system-wide:

```sh
sudo apt-get install libuv1-dev          # Debian / Ubuntu
sudo dnf install libuv-devel             # Fedora / RHEL
```

…or skip the system dep entirely by exporting `USE_BUNDLED_LIBUV=1` before
installing — `fs` will then build its own bundled libuv. This is what the
project's CI workflows do.

### Dependencies

semseeker relies on several CRAN packages. Install them before the first use:

```r
install.packages(c(
  "doFuture", "FactoMineR", "factoextra", "FSA", "fst",
  "future", "future.apply", "ggplot2", "gridExtra", "gtools",
  "Hmisc", "lqmm", "progressr", "reshape2", "zoo",
  "VennDiagram", "effsize", "pwr", "combinat", "philentropy",
  "parallelly", "knitr", "lmtest", "usethis", "wordcloud",
  "WebGestaltR", "sqldf", "caret", "openai", "glmnet"
))
```

---

## Quick Example

### Step 1 — Load and normalise IDAT files with ChAMP

```r
library(ChAMP)

idat_folder  <- "~/source_idat/"
result_folder <- "~/result/"

myLoad <- champ.load(
  directory      = idat_folder,
  method         = "minfi",
  methValue      = "B",
  autoimpute     = TRUE,
  filterDetP     = TRUE,
  ProbeCutoff    = 0,
  SampleCutoff   = 0.1,
  detPcut        = 0.01,
  filterBeads    = TRUE,
  beadCutoff     = 0.05,
  filterNoCG     = TRUE,
  filterSNPs     = TRUE,
  filterMultiHit = TRUE,
  filterXY       = TRUE,
  arraytype      = "EPIC"   # or "450K"
)

myNorm <- champ.norm(
  beta       = myLoad$beta,
  rgSet      = myLoad$rgSet,
  mset       = myLoad$mset,
  resultsDir = result_folder,
  method     = "SWAN",
  arraytype  = "EPIC",
  cores      = parallel::detectCores() - 1
)

saveRDS(myNorm, "~/normalizedData.rds")
```

### Step 2 — Run semseeker

```r
library(semseeker)

normalizedData <- readRDS("~/normalizedData.rds")
sample_sheet   <- read.csv2("~/sample_sheet.csv")

semseeker(
  sample_sheet  = sample_sheet,
  signal_data   = normalizedData,
  result_folder = "~/semseeker_result/"
)
```

Results (BED files, pivot tables, delta-metric summaries) are written to
`result_folder`.

---

## Complete Example

A fully worked example using public data from Gene Expression Omnibus (GEO) is
available in the repository's `example/` folder.

---

## Input Requirements

### Beta-value matrix

- Rows = probes (Illumina probe IDs, e.g. `cg...`)
- Columns = samples
- Values = beta values in [0, 1] or M-values (detected automatically)

### Sample sheet

A data frame (or CSV readable with `read.csv2`) with **at minimum** these columns:

| Column | Description |
|---|---|
| `Sample_ID` | Unique sample identifier (must match matrix column names) |
| `Sample_Name` | Human-readable label |
| `Sample_Group` | Must be one of `"Case"`, `"Control"`, `"Reference"` |

> **Note:** If no Reference population is available, duplicate the Control
> rows and assign `Sample_Group = "Reference"` to the duplicates.

> **Note:** `Sample_ID` values are uppercased internally by semseeker
> (`name_cleaning()`); output file names and pivot column headers reflect the
> uppercased form.

---

## Output

For each sample semseeker writes to `result_folder/Data/`:

- `*_MUTATIONS_HYPO.bed` / `*_MUTATIONS_HYPER.bed` — probe-level mutations
- `*_LESIONS_HYPO.bed` / `*_LESIONS_HYPER.bed` — genomic lesion clusters
- `*_DELTAS_*.bed` / `*_DELTAR_*.bed` — per-region signed deviation (`DELTAS`) and its ratio to the reference interval (`DELTAR`)
- Pivot tables aggregating all samples per marker

---

## Citation

If you use semseeker in published work, please cite the package:

```r
citation("semseeker")
```

The package is archived on Zenodo under the concept DOI
[10.5281/zenodo.5095416](https://doi.org/10.5281/zenodo.5095416), which always
resolves to the latest archived version; every release tag also gets its own
version DOI.
