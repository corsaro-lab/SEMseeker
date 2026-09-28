# Getting started — epimutations and delta metrics

## What is semseeker?

**semseeker** is an R package for detecting, quantifying and analysing
**Stochastic Epigenetic Mutations (SEMs)** from Illumina array-based DNA
methylation data (Corsaro 2025).

### Background: what is a SEM?

A **Stochastic Epigenetic Mutation** is a rare, sample-specific
aberration in DNA methylation at a single CpG site — a probe that
deviates from the population-level signal in an individual sample beyond
what can be explained by normal biological or technical variability
(Gentilini et al. 2015).

SEMs were first formally defined and quantified by Gentilini *et al.*
(2015), who showed that SEM burden increases exponentially with age
(Gentilini et al. 2015). Subsequent work (Corsaro et al. 2023; Spada et
al. 2020) extended the framework to disease association studies,
introducing the concept of *lesions* (genomic clusters of co-occurring
SEMs) and *delta metrics* for population-level variability summaries.

**semseeker** implements this methodology and wraps it into a fully
reproducible R workflow, from raw beta-value matrix to association
analysis and pathway enrichment (Corsaro 2025; Corsaro et al. 2023).

### What semseeker computes

Starting from a normalised beta-value (or M-value) matrix, semseeker:

1.  Derives **signal thresholds** (IQR-based) per probe from the
    reference population.
2.  Identifies **hypomethylated** and **hypermethylated** SEMs at
    single-probe resolution.
3.  Aggregates adjacent SEMs into **lesions** (genomic clusters of
    co-occurring mutations).
4.  Computes per-region **delta metrics**: the signed deviation beyond
    the reference interval (`DELTAS`), its ratio to the interval width
    (`DELTAR`), and their binned/ranked variants (`DELTAP`/`DELTARP`,
    `DELTAQ`/`DELTARQ`).
5.  Optionally analyses the **raw methylation signal** directly (mean
    per region, `SIGNAL_MEAN`) for differential methylation without
    mutation calling.

Supported input platforms: Illumina EPIC (850k), 450k, and 27k arrays
(detected automatically from probe IDs).

------------------------------------------------------------------------

## Installation

``` r

# Install remotes if needed
install.packages("remotes")

# Install SEMseeker from GitHub
remotes::install_github("drake69/SEMseeker")
```

> **Note:** semseeker requires the `polars` package from R-multiverse.
> It is installed automatically as a dependency. If you get an
> installation error, try:
>
> ``` r
>
> Sys.setenv(NOT_CRAN = "true")
> install.packages("polars", repos = "https://community.r-multiverse.org")
> ```

### System dependencies

**macOS — XQuartz.** R 4.6’s base `tcltk` package is linked against
`/opt/X11/lib/libX11.6.dylib` at install time. Plain SEMseeker usage
([`library(SEMseeker)`](https://github.com/drake69/semseeker),
`semseeker(...)`, `association_analysis(...)`) does not load tcltk.
However, some Bioconductor optional features that transitively pull
`minfi` or `GEOquery` (e.g. running the annotation-concordance checks)
load tcltk and segfault if XQuartz is not installed. Install XQuartz
once via:

``` sh
brew install --cask xquartz
```

or download from <https://www.xquartz.org/>. You may need to log out and
back in for the path to be picked up. Linux and Windows users do not
need XQuartz.

**Linux — libuv.** Building from source on Linux pulls `fs` which needs
`libuv` ≥ 1 at compile time. Install via
`sudo apt-get install libuv1-dev` (Debian/Ubuntu) or
`sudo dnf install libuv-devel` (Fedora/RHEL), or set
`USE_BUNDLED_LIBUV=1` to compile against a bundled copy.

------------------------------------------------------------------------

## Step 1 — Prepare your data

### Beta-value matrix

semseeker expects a **probe × sample** matrix of beta values (0–1) or
M-values. Probe IDs must be Illumina CpG identifiers
(e.g. `cg00000029`).

``` r

# Example: load a pre-normalised matrix (rows = probes, columns = samples).
# Here we use the small synthetic matrix prepared in the hidden setup
# chunk above; in practice you would readRDS("normalised_beta.rds").

dim(beta_matrix)       # e.g. 866895 × 50 on a real EPIC dataset
#> [1] 77 10
beta_matrix[1:3, 1:3]
#>               CTRL01    CTRL02    CTRL03
#> cg11753499 0.4156308 0.3945168 0.4139144
#> cg13210239 0.6379258 0.6371134 0.6563844
#> cg17985533 0.5566887 0.5739348 0.5200354
```

If you are starting from raw IDAT files, normalise first with
[ChAMP](https://bioconductor.org/packages/ChAMP/):

``` r

library(ChAMP)

myLoad <- champ.load(
  directory  = "path/to/idat/",
  arraytype  = "EPIC"
)

myNorm <- champ.norm(
  beta      = myLoad$beta,
  rgSet     = myLoad$rgSet,
  mset      = myLoad$mset,
  method    = "SWAN",
  arraytype = "EPIC"
)

beta_matrix <- myNorm
```

### Sample sheet

A data frame with at minimum:

| Column | Description |
|----|----|
| `Sample_ID` | Unique sample identifier — must match column names of the beta matrix |
| `Sample_Name` | Human-readable label |
| `Sample_Group` | One of `"Case"`, `"Control"`, `"Reference"` |

``` r

# Example: load a sample sheet (here we use the synthetic one prepared
# in the hidden setup chunk; in practice you would read.csv2("sample_sheet.csv")).

# Minimal required columns
head(sample_sheet[, c("Sample_ID", "Sample_Name", "Sample_Group")])
#>   Sample_ID Sample_Name Sample_Group
#> 1    CTRL01      CTRL01    Reference
#> 2    CTRL02      CTRL02    Reference
#> 3    CTRL03      CTRL03    Reference
#> 4    CTRL04      CTRL04    Reference
#> 5    CTRL05      CTRL05    Reference
#> 6    CTRL06      CTRL06    Reference
```

> **Tip:** If you have no independent reference population, duplicate
> the Control rows and set `Sample_Group = "Reference"` for the
> duplicates.

------------------------------------------------------------------------

## Step 2 — Detect epimutations

The main entry point is
[`semseeker()`](https://drake69.github.io/semseeker/reference/semseeker.md).
A single call runs the full pipeline: normalisation → mutation calling →
lesion aggregation → delta metric computation.

``` r

library(semseeker)

semseeker(
  input         = beta_matrix,
  sample_sheet  = sample_sheet,
  result_folder = "~/semseeker_results/"
)
```

### Key parameters

| Parameter | Default | Description |
|----|----|----|
| `parallel_strategy` | `"multicore"` | Parallelisation backend (`"sequential"`, `"multicore"`, `"multisession"`) |
| `maxResources` | `90` | % of CPU cores to use |
| `iqrTimes` | `3` | IQR multiplier for aberration threshold |
| `alpha` | `0.05` | Significance threshold |
| `areas` | `c("GENE","CHR","ISLAND","PROBE")` | Genomic aggregation levels |
| `marker` | `c("MUTATIONS","LESIONS","DELTARQ",...)` | Which markers to compute |
| `sex_chromosome_remove` | `TRUE` | Remove sex chromosomes from analysis |

Example — run only probe-level mutations, sequential mode, 2 cores:

``` r

semseeker(
  input             = beta_matrix,
  sample_sheet      = sample_sheet,
  result_folder     = "~/semseeker_results/",
  parallel_strategy = "sequential",
  areas             = c("PROBE"),
  marker            = c("MUTATIONS", "LESIONS")
)
```

### Genomic areas & subareas

After per-probe markers are computed, SEMseeker annotates each probe to
the genomic **areas** (and their **subareas**) it belongs to, so that
epimutations can be aggregated region-by-region. The available subareas
per area:

| Area | Subareas |
|----|----|
| `GENE` | `BODY`, `TSS1500`, `TSS200`, `1STEXON`, `5UTR`, `3UTR`, `EXONBND`, `WHOLE` |
| `ISLAND` | `ISLAND`, `N_SHORE`, `S_SHORE`, `N_SHELF`, `S_SHELF`, `OPENSEA`, `WHOLE` |
| `CHR` | `WHOLE`, `CYTOBAND` |
| `DMR` | `WHOLE`, `DMR` |
| `PROBE` | `WHOLE` |

#### CpG-island context mapping (`ISLAND`)

The `ISLAND` subareas mirror Illumina’s `Relation_to_Island`, which
classifies every probe into one of six mutually exclusive categories.
SEMseeker maps them as follows — and, exactly as `GENE_WHOLE` covers the
**whole gene**, `ISLAND_WHOLE` covers the **whole island
neighbourhood**:

| Illumina `Relation_to_Island` | SEMseeker subarea | Region label |
|----|----|----|
| `Island` | `ISLAND` | island coordinate (`chr:start-end`) |
| `N_Shore` / `S_Shore` | `N_SHORE` / `S_SHORE` | island coordinate (±2 kb flank) |
| `N_Shelf` / `S_Shelf` | `N_SHELF` / `S_SHELF` | island coordinate (2–4 kb flank) |
| `OpenSea` | `OPENSEA` | gap coordinate (`chr:start-end`) |
| *(union of the above except OpenSea)* | `WHOLE` | island coordinate |

Notes:

- **`WHOLE` = `ISLAND` + shores + shelves** (the whole CpG-island
  neighbourhood, island core ±4 kb). It is **not** the island core alone
  — use `ISLAND` for the core.
- **`OPENSEA`** are the CpGs outside every neighbourhood. Each open-sea
  probe is grouped by the genomic **gap** between neighbourhoods that
  contains it, labelled by that gap’s coordinate `chr:start-end`. A gap
  never crosses a chromosome boundary, so no open-sea region spans
  chromosomes.
- The same definitions apply to WGBS / long-read data, where the island
  neighbourhoods are reconstructed from UCSC CpG-island coordinates
  instead of the Illumina manifest, so array and sequencing analyses
  share identical region semantics.

------------------------------------------------------------------------

## Step 3 — Understand the output

Results are written under `result_folder/Data/`. Each sample produces:

    Data/
    ├── SAMPLEID_MUTATIONS_HYPO.bed       ← probe-level hypomethylated mutations
    ├── SAMPLEID_MUTATIONS_HYPER.bed      ← probe-level hypermethylated mutations
    ├── SAMPLEID_LESIONS_HYPO.bed         ← genomic lesion clusters (hypo)
    ├── SAMPLEID_LESIONS_HYPER.bed        ← genomic lesion clusters (hyper)
    ├── SAMPLEID_SIGNAL_MEAN.PROBE.bedgraph
    └── ...

Population-level pivot tables aggregate all samples per marker:

    Data/Pivots/
    ├── SIGNAL/SIGNAL_MEAN_PROBE_WHOLE_hg19.parquet       ← mean signal per probe, all samples
    ├── MUTATIONS/MUTATIONS_HYPO_PROBE_WHOLE_hg19.parquet
    ├── DELTARP/DELTARP_HYPER_GENE_WHOLE_hg19.parquet     ← DELTARP (hyper) per gene, all samples
    └── ...

Pivot file names follow
`<MARKER>_<FIGURE>_<AREA>_<SUBAREA>_<genome_build>.parquet` under
`Data/Pivots/<MARKER>/`. The marker name is invariant across aggregation
levels — only the `AREA` segment (`PROBE`, `GENE`, `ISLAND`, …) changes.

### Reading results

``` r

library(polars)

# Read the DELTARP pivot table (the marker name is identical at every
# aggregation level; only the AREA segment changes, e.g. GENE vs PROBE)
deltarp <- as.data.frame(
  pl$read_parquet("~/semseeker_results/Data/Pivots/DELTARP/DELTARP_HYPER_GENE_WHOLE_hg19.parquet")
)

head(deltarp[, 1:6])
```

------------------------------------------------------------------------

## Step 4 — Delta metrics explained

semseeker computes delta metrics for each genomic position relative to a
reference signal interval `[inferior, superior]`, derived per probe from
the control population as IQR-based lower and upper thresholds:

| Marker | Full name | Captures |
|----|----|----|
| `DELTAS` | Delta Signal | Signed deviation beyond the reference threshold (`value − superior` for hyper, `inferior − value` for hypo) |
| `DELTAR` | Delta Ratio | `DELTAS` divided by the reference interval width (`superior − inferior`) — the deviation as a proportion of the reference range |
| `DELTAP` / `DELTARP` | Delta (ratio) ranked, equal-width | `DELTAS` / `DELTAR` discretised into equal-width bins; each position gets an integer rank weight `1..B` (default `B = 4`) |
| `DELTAQ` / `DELTARQ` | Delta (ratio) ranked, quantile | `DELTAS` / `DELTAR` discretised into quantile bins; each position gets an integer rank weight `1..Q` (default `Q = 4`) |

`DELTARP` is the recommended metric for downstream analysis: by ranking
the relative deviation (`DELTAR`) into equal-width bins it is robust to
outliers and directly comparable across regions. The marker name is the
same at every aggregation level (locus, gene, island, …) — only the AREA
segment of the pivot file name changes.

``` r

# Load DELTARP at gene level (hyper figure)
deltarp_gene <- as.data.frame(
  pl$read_parquet("~/semseeker_results/Data/Pivots/DELTARP/DELTARP_HYPER_GENE_WHOLE_hg19.parquet")
)

# Top regions by mean DELTARP across cases
case_cols  <- grep("^CASE", names(deltarp_gene), value = TRUE)
deltarp_gene$mean_deltarp <- rowMeans(deltarp_gene[, case_cols], na.rm = TRUE)

head(deltarp_gene[order(-deltarp_gene$mean_deltarp),
                  c("CHR", "START", "END", "mean_deltarp")], 10)
```

------------------------------------------------------------------------

------------------------------------------------------------------------

## Step 5 — Differential signal analysis

Beyond mutation calling, semseeker also tracks the **raw methylation
signal** per region via the `SIGNAL_MEAN` marker:

| Marker | File pattern | Description |
|----|----|----|
| `SIGNAL_MEAN` | `Pivots/SIGNAL/SIGNAL_MEAN_*.parquet` | Per-region mean beta value per sample |

This pivot table has exactly the same structure as the mutation/delta
tables and can be passed directly to
[`association_analysis()`](https://drake69.github.io/semseeker/reference/association_analysis.md)
— allowing you to test whether mean methylation differs by phenotype,
without any mutation-calling step.

``` r

# inference_details specifies the PREDICTOR and the test — not a model formula.
# `independent_variable` is the phenotype column, `family_test` the statistical
# test. The dependent marker (SIGNAL_MEAN, MUTATIONS, …) is NOT named here: it
# comes from the markers you ran through semseeker(); association_analysis()
# iterates over them and reads the matching pivot tables internally.
inference_details <- data.frame(
  independent_variable = "Sample_Group",
  family_test          = "wilcoxon",   # two independent groups, non-parametric
  transformation_y     = "none",
  depth_analysis       = 3,            # per genomic area (region-level / EWAS-style)
  filter_p_value       = TRUE
)

association_analysis(
  inference_details = inference_details,
  result_folder     = "~/semseeker_results/",
  parallel_strategy = "multicore",
  maxResources      = 90
)
```

This is particularly useful when you want to detect **differential
methylation** at the region level (classic EWAS-style analysis)
alongside the epimutation burden analysis, using a single unified
framework.

------------------------------------------------------------------------

## Subset rationale & full-pipeline runnable

The runnable demo above uses a **subset** of probes (only those around
the BWS-related imprinting DMRs `KCNQ1OT1:TSS-DMR` and
`H19_IGF2:IG-DMR`) so the vignette builds in seconds during
CRAN/Bioconductor checks. This is **not** representative of a real-world
workflow: a production analysis would normalise IDAT files for an entire
cohort, correct for biological confounders (age, cell composition for
whole-blood samples), and then run semseeker on the full beta matrix.

The block below — kept as commented `eval = FALSE` reference — shows the
**complete end-to-end pipeline** as you would run it on a real GEO
dataset such as Beckwith-Wiedemann GSE95486. Uncomment and adapt to your
study.

``` r

# ─── Full pipeline from raw GEO data (~30-60 min on a laptop) ──────────────

## Step 0 — Download IDAT files from GEO (e.g. BWS GSE95486)
library(GEOquery)
gse <- getGEO("GSE95486", GSEMatrix = FALSE, destdir = "./raw")
GEOquery::getGEOSuppFiles("GSE95486", baseDir = "./raw", makeDirectory = TRUE)
# (then extract the IDATs into ./raw/GSE95486/IDATs/)

## Step 1 — Load & SWAN-normalise IDAT files via ChAMP
library(ChAMP)
myLoad <- champ.load(
  directory      = "./raw/GSE95486/IDATs",
  arraytype      = "EPIC",
  filterDetP     = TRUE,
  filterSNPs     = TRUE,
  filterMultiHit = TRUE,
  filterXY       = TRUE
)
myNorm <- champ.norm(
  beta       = myLoad$beta,
  rgSet      = myLoad$rgSet,
  mset       = myLoad$mset,
  method     = "SWAN",
  arraytype  = "EPIC",
  cores      = parallel::detectCores() - 1L
)

## Step 2 — Age correction
##   Confounding by chronological age is the strongest known driver of SEM
##   burden (Gentilini et al., 2015 — SEMs accumulate exponentially with
##   age). Always regress age out before disease-vs-control comparison.
library(limma)
age_residuals <- removeBatchEffect(
  myNorm,
  covariates = sample_sheet$Age   # numeric, years
)

## Step 3 — Cell composition deconvolution (whole-blood samples ONLY)
##   In peripheral whole blood, methylation is strongly confounded by the
##   relative proportions of B cells, CD4/CD8 T, monocytes, NK and granulocytes.
##   For tissues other than whole blood, skip this step.
library(FlowSorted.Blood.EPIC)
cells <- estimateCellCounts2(
  myLoad$rgSet,
  compositeCellType = "Blood",
  referencePlatform = "IlluminaHumanMethylationEPIC"
)
# Regress out the cell-fraction covariates from the (age-residualised) matrix
adjusted <- removeBatchEffect(
  age_residuals,
  covariates = cells$counts
)

## Step 4 — Run SEMseeker on the fully preprocessed matrix
SEMseeker::semseeker(
  input         = adjusted,
  sample_sheet  = sample_sheet,
  result_folder = "./bws_full_results/"
)
# ──────────────────────────────────────────────────────────────────────────
```

**Pre-processing checklist** before running
[`semseeker()`](https://drake69.github.io/semseeker/reference/semseeker.md):

Cohort normalisation (SWAN / BMIQ / Functional Normalisation)

Detection-p filter, multi-hit/SNP/cross-reactive probe removal

**Age** residualisation (mandatory if cohort age range \> a few years)

**Cell composition** deconvolution (whole-blood only —
`FlowSorted.Blood.EPIC`)

Batch / processing-date correction (if known)

Sample QC pass (detection rate, bead count, predicted-sex check)

semseeker assumes upstream preprocessing has already been performed.
Skipping these steps inflates the false-positive SEM rate — particularly
in older cohorts and whole-blood studies.

------------------------------------------------------------------------

## References

------------------------------------------------------------------------

## Next steps

- **Association analysis** (all statistical models): see the
  [Association analysis
  vignette](https://drake69.github.io/semseeker/articles/association-analysis.md).
- **Pathway and enrichment analysis**: see the [Pathway analysis
  vignette](https://drake69.github.io/semseeker/articles/pathway-analysis.md).
- **Full function reference**: see the
  [Reference](https://drake69.github.io/semseeker/reference/index.md)
  page.

------------------------------------------------------------------------

## Session info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: aarch64-apple-darwin23
#> Running under: macOS Tahoe 26.6.2
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: UTC
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.61         cachem_1.1.0      knitr_1.52        htmltools_0.5.9  
#>  [9] rmarkdown_2.32    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.1     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.6.1    tools_4.6.1       ragg_1.5.2        bslib_0.12.0     
#> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
#> [25] rlang_1.3.0       fs_2.1.0          htmlwidgets_1.6.4
```

Corsaro, Luigi. 2025. “SEMSeeker: Un Pacchetto R Per Condurre Studi Di
Associazione Epigenetici Basati Su Epimutazione Stocastica.” PhD thesis,
Università degli Studi di Pavia.
<https://hdl.handle.net/20.500.14242/189932>.

Corsaro, Luigi, Davide Gentilini, Luciano Calzari, and Vincenzo Gambino.
2023. “Notch, SUMOylation, and ESR-Mediated Signalling Are the Main
Molecular Pathways Showing Significantly Different Epimutation Scores
Between Expressing or Not Oestrogen Receptor Breast Cancer in Three
Public EWAS Datasets.” *Cancers* 15 (16): 4109.
<https://doi.org/10.3390/cancers15164109>.

Gentilini, Davide, Paolo Garagnani, Serena Pisoni, et al. 2015.
“Stochastic Epigenetic Mutations (DNA Methylation) Increase
Exponentially in Human Aging and Correlate with X Chromosome
Inactivation Skewing in Females.” *Aging* 7 (8): 568–78.
<https://doi.org/10.18632/aging.100792>.

Spada, Elena, Luciano Calzari, Luigi Corsaro, et al. 2020. “Epigenome
Wide Association and Stochastic Epigenetic Mutation Analysis on Cord
Blood of Preterm Birth.” *International Journal of Molecular Sciences*
21 (14): 5044. <https://doi.org/10.3390/ijms21145044>.
