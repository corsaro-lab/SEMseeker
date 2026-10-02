# Pathway and enrichment analysis

## Overview

After
[`association_analysis()`](https://corsaro-lab.github.io/SEMseeker/reference/association_analysis.md)
has identified genomic regions whose SEM metrics (ΔRP, mutation counts,
…) are statistically associated with a phenotype, the next step is to
ask: **which biological pathways and molecular networks are enriched
among those regions?**

semseeker provides six pathway/enrichment functions that all accept the
output of
[`association_analysis()`](https://corsaro-lab.github.io/SEMseeker/reference/association_analysis.md)
as input and run the enrichment for every combination of genomic area,
marker type and sample group that was analysed:

| Function | Backend | Type of analysis |
|----|----|----|
| `enrich_WebGestalt()` | WebGestaltR (Liao et al. 2019) | Gene Ontology ORA / GSEA |
| `enrich_STRINGdb()` | STRINGdb (Szklarczyk et al. 2023) | Protein-protein interaction network |
| `enrich_pathfindR()` | pathfindR (Ulgen et al. 2019) | Active subnetwork enrichment |
| `enrich_ctdR()` | ctdR / CTD (Davis et al. 2023) | Chemical-gene-disease associations |
| `enrich_Phenolyzer_WebGestalt()` | Phenolyzer + WebGestaltR | Disease-prioritised ORA / GSEA |
| `enrich_Phenolyzer_STRINGdb()` | Phenolyzer + STRINGdb | Disease-prioritised network |

All functions are **session-aware**: they read the active semseeker
session (initialised by
[`association_analysis()`](https://corsaro-lab.github.io/SEMseeker/reference/association_analysis.md))
and therefore do **not** require passing `result_folder` again — they
pick it up automatically via `core_get_session_info()`.

------------------------------------------------------------------------

## Workflow

    semseeker()               →  epimutation detection + delta metrics
          ↓
    association_analysis()    →  statistical associations (CSV in Inference/)
          ↓
    pathway_*()               →  enriched pathways / networks (CSV in Pathway/)
          ↓
    meta_enrichment_compare_studies()       →  cross-cohort Venn diagrams
    meta_enrichment_overlaps_subsamples()   →  stability across subsamples

The `inference_details` object passed to the pathway functions is the
**same** data frame used in
[`association_analysis()`](https://corsaro-lab.github.io/SEMseeker/reference/association_analysis.md).
It tells the pathway functions which statistical results to use as the
gene list for enrichment.

------------------------------------------------------------------------

## Installing optional backends

Each pathway function requires a separate optional package. Install only
what you need:

``` r

# WebGestaltR — ORA / GSEA
install.packages("WebGestaltR")

# STRINGdb — protein network (Bioconductor)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("STRINGdb")

# pathfindR — active subnetwork enrichment
install.packages("pathfindR")

# ctdR — chemical-disease associations (GitHub)
remotes::install_github("drake69/ctdR")
```

> All backends are in `Suggests` — semseeker installs and runs without
> them. Each function checks for its backend with
> [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) and logs
> an informative ERROR if it is missing.

------------------------------------------------------------------------

## 1. `enrich_WebGestalt()` — Gene Ontology ORA and GSEA

WebGestalt (Wang et al. 2017; Liao et al. 2019) runs Over-Representation
Analysis (ORA) (Boyle et al. 2004) or Gene Set Enrichment Analysis
(GSEA) (Subramanian et al. 2005) against Gene Ontology (Ashburner et al.
2000), KEGG (Kanehisa and Goto 2000), Reactome (Jassal et al. 2020), and
other databases.

``` r

library(semseeker)

# Run after association_analysis() — session is already active
inference_details <- data.frame(
  formula               = "DELTARP ~ age",
  family                = "spearman",
  transformation_y      = "none",
  depth_analysis        = 3,
  filter_p_value        = TRUE,
  stringsAsFactors      = FALSE
)

association_analysis(
  inference_details = inference_details,
  result_folder     = "~/semseeker_results/",
  parallel_strategy = "multicore"
)

# Gene Ontology ORA on Biological Process and Molecular Function
enrich_WebGestalt(
  study          = "my_study",
  types          = c("BP", "MF"),       # GO namespaces: BP, MF, CC
  enrich_methods = c("ORA"),            # "ORA" or "GSEA"
  pvalue_column  = "PVALUE_ADJ_ALL_BH", # which association p-value to use
  inference_detail = inference_details[1, ],
  significance   = TRUE                 # only use significant genes
)

# GSEA on KEGG pathways
enrich_WebGestalt(
  study          = "my_study",
  types          = c("KEGG"),
  enrich_methods = c("GSEA"),
  pvalue_column  = "PVALUE_ADJ_ALL_BH",
  inference_detail = inference_details[1, ],
  significance   = TRUE
)
```

### Parameters

| Parameter | Default | Description |
|----|----|----|
| `study` | — | Study label, used in output file names |
| `types` | `c("BP","MF")` | GO namespaces or databases: `"BP"`, `"MF"`, `"CC"`, `"KEGG"`, `"Reactome"` |
| `enrich_methods` | `c("ORA")` | `"ORA"` (over-representation) or `"GSEA"` (ranked) |
| `pvalue_column` | `"PVALUE_ADJ_ALL_BH"` | Column from association results used to select / rank genes |
| `inference_detail` | — | One row from `inference_details` |
| `significance` | — | `TRUE` = use only FDR-significant regions as gene list |
| `adjust_per_area` | `FALSE` | Re-adjust p-values within each genomic area |
| `adjust_globally` | `FALSE` | Re-adjust p-values globally across all areas |
| `adjustment_method` | `"BH"` | Multiple testing correction method |

### Output

Results land in `result_folder/Pathway/WebGestalt/` with one CSV per
analysis key. After `enrich_analysy_add_category()` normalisation, each
result has a unified `SS_CATEGORY` column:

| `SS_CATEGORY` | Source                           |
|---------------|----------------------------------|
| `GO-BP`       | Gene Ontology Biological Process |
| `GO-MF`       | Gene Ontology Molecular Function |
| `GO-CC`       | Gene Ontology Cellular Component |
| `KEGG-ORA`    | KEGG via ORA                     |

------------------------------------------------------------------------

## 2. `enrich_STRINGdb()` — Protein-protein interaction network

STRINGdb (Szklarczyk et al. 2019, 2023) maps the significant genes onto
the STRING protein-protein interaction network and identifies enriched
functional categories, GO terms, KEGG pathways and tissue expressions.

``` r

enrich_STRINGdb(
  study              = "my_study",
  inference_details  = inference_details,
  pvalue_column      = "PVALUE_ADJ_ALL_BH",
  significance       = TRUE,
  stringDBVersion    = "12.0",           # STRING database version
  statistic_parameter = ""               # "" = presence/absence; column name = use score
)
```

### Parameters

| Parameter | Default | Description |
|----|----|----|
| `study` | — | Study label |
| `inference_details` | — | Full `inference_details` data frame (multi-row supported) |
| `pvalue_column` | `"PVALUE_ADJ_ALL_BH"` | Association p-value column |
| `significance` | `TRUE` | Use only significant regions |
| `stringDBVersion` | `"12.0"` | STRING database version (`"11.5"`, `"12.0"`) |
| `statistic_parameter` | `""` | `""` = unweighted; column name = use that value as interaction weight |
| `adjust_per_area` / `adjust_globally` | `FALSE` | Additional p-value correction |

### Network score threshold

The default score threshold is **200** (medium confidence). Increase to
400 or 700 for higher confidence interactions:

``` r

# Higher confidence network (score >= 700)
# Modify via the STRINGdb$new() call inside the function — currently fixed at 200
# Use statistic_parameter to weight genes by their association score:
enrich_STRINGdb(
  study               = "my_study",
  inference_details   = inference_details,
  statistic_parameter = "SCORE",   # weight by Spearman rho or regression coefficient
  significance        = TRUE
)
```

### Output categories

| `SS_CATEGORY` | Description             |
|---------------|-------------------------|
| `GO-BP`       | GO Biological Process   |
| `GO-MF`       | GO Molecular Function   |
| `GO-CC`       | GO Cellular Component   |
| `KEGG-ORA`    | KEGG pathway enrichment |

------------------------------------------------------------------------

## 3. `enrich_pathfindR()` — Active subnetwork enrichment

pathfindR (Ulgen et al. 2019) identifies enriched pathways by searching
for **active subnetworks** in a protein interaction network (PIN) rather
than simple gene list overlap. This approach captures pathway activity
even when only a subset of pathway genes pass the significance
threshold.

``` r

enrich_pathfindR(
  study              = "my_study",
  path_dbs           = c("KEGG", "Reactome", "GO-BP", "GO-MF", "GO-CC"),
  iterations         = 20,              # number of active subnetwork iterations
  inference_details  = inference_details,
  pvalue_column      = "PVALUE_ADJ_ALL_BH",
  significance       = TRUE,
  statistic_parameter = ""
)
```

### Parameters

| Parameter | Default | Description |
|----|----|----|
| `study` | — | Study label |
| `path_dbs` | — | Databases: `"KEGG"`, `"Reactome"`, `"GO-BP"`, `"GO-MF"`, `"GO-CC"` |
| `iterations` | `20` | Number of active subnetwork search iterations — higher = more stable but slower |
| `inference_details` | — | Full `inference_details` data frame |
| `pvalue_column` | `"PVALUE_ADJ_ALL_BH"` | Association p-value column |
| `significance` | `TRUE` | Filter to significant regions |
| `statistic_parameter` | `""` | Optional: column to use as gene-level score |

### Output categories

| `SS_CATEGORY`  | Description                               |
|----------------|-------------------------------------------|
| `GO-BP`        | GO Biological Process (active subnetwork) |
| `GO-MF`        | GO Molecular Function (active subnetwork) |
| `KEGG-NTA`     | KEGG via network topology (NTA)           |
| `REACTOME-NTA` | Reactome via network topology             |

### Installation note

pathfindR requires `ggkegg` (Bioconductor). Install the full chain:

``` r

BiocManager::install("ggkegg")
install.packages("pathfindR")
```

------------------------------------------------------------------------

## 4. `enrich_ctdR()` — Chemical-gene-disease associations

`enrich_ctdR()` uses the **Comparative Toxicogenomics Database** (CTD)
(Davis et al. 2023) via the `ctdR` package to find chemical compounds
and diseases associated with the significant genes. This is particularly
useful for epigenetic studies involving environmental exposures.

``` r

enrich_ctdR(
  study              = "my_study",
  inference_details  = inference_details,
  pvalue_column      = "PVALUE_ADJ_ALL_BH",
  significance       = TRUE,
  statistic_parameter = ""
)
```

### Parameters

| Parameter | Default | Description |
|----|----|----|
| `study` | — | Study label |
| `inference_details` | — | Full `inference_details` data frame |
| `pvalue_column` | `"PVALUE_ADJ_ALL_BH"` | Association p-value column |
| `significance` | `TRUE` | Use only significant regions |
| `statistic_parameter` | `""` | Optional signal weight column |

### Output

Results are normalised to `SS_CATEGORY = "CHEMICAL"`, linking
significant epigenetic regions to chemical exposures and associated
diseases.

------------------------------------------------------------------------

## 5 & 6. Phenolyzer variants — disease-prioritised enrichment

The Phenolyzer variants (`enrich_Phenolyzer_WebGestalt()` and
`enrich_Phenolyzer_STRINGdb()`) first **prioritise genes** using the
Phenolyzer disease-gene scoring system and then run enrichment on the
prioritised list.

Use these when you have a specific disease hypothesis and want to focus
the enrichment on genes already known to be relevant to that disease:

``` r

# Phenolyzer + WebGestalt
enrich_Phenolyzer_WebGestalt(
  study            = "my_study",
  disease          = "breast cancer",   # disease name for Phenolyzer prioritisation
  types            = c("BP", "MF"),
  enrich_methods   = c("ORA"),
  pvalue_column    = "PVALUE_ADJ_ALL_BH",
  inference_detail = inference_details[1, ],
  significance     = TRUE
)

# Phenolyzer + STRINGdb
enrich_Phenolyzer_STRINGdb(
  study               = "my_study",
  disease             = "breast cancer",
  inference_detail    = inference_details,
  pvalue_column       = "PVALUE_ADJ_ALL_BH",
  significance        = TRUE,
  statistic_parameter = ""
)
```

### Additional parameter: `disease`

| Parameter | Description |
|----|----|
| `disease` | Free-text disease name passed to Phenolyzer for gene prioritisation (e.g. `"lung cancer"`, `"Alzheimer"`, `"TCDD exposure"`) |

------------------------------------------------------------------------

## Running all backends in one script

``` r

library(semseeker)

result_folder <- "~/semseeker_results/"

# 1. SEM detection
semseeker(input = beta_matrix, sample_sheet = sample_sheet, result_folder = result_folder)

# 2. Association analysis
inference_details <- data.frame(
  formula          = "DELTARP ~ case_control",
  family           = "wilcoxon",
  transformation_y = "none",
  depth_analysis   = 3,
  filter_p_value   = TRUE,
  stringsAsFactors = FALSE
)
association_analysis(inference_details, result_folder)

# 3. Pathway enrichment — run all backends
enrich_WebGestalt(
  study="study1", types=c("BP","MF","KEGG"), enrich_methods=c("ORA"),
  pvalue_column="PVALUE_ADJ_ALL_BH",
  inference_detail=inference_details[1,], significance=TRUE)

enrich_STRINGdb(
  study="study1", inference_details=inference_details,
  pvalue_column="PVALUE_ADJ_ALL_BH", significance=TRUE,
  stringDBVersion="12.0")

enrich_pathfindR(
  study="study1", path_dbs=c("KEGG","Reactome","GO-BP"),
  iterations=20, inference_details=inference_details,
  pvalue_column="PVALUE_ADJ_ALL_BH", significance=TRUE)

enrich_ctdR(
  study="study1", inference_details=inference_details,
  pvalue_column="PVALUE_ADJ_ALL_BH", significance=TRUE)
```

------------------------------------------------------------------------

## Output structure

All results are written under `result_folder/Pathway/`:

    Pathway/
    ├── WebGestalt/
    │   └── <samples_filter>/<inference_filter>/
    │       └── <area>_<marker>_<figure>_<test>_<variable>.csv
    ├── STRINGdb/
    │   └── ...
    ├── pathfindR/
    │   └── ...
    ├── ctdR/
    │   └── ...
    ├── Phenolyzer_WebGestalt/
    │   └── ...
    └── Phenolyzer_STRINGdb/
        └── ...

Each CSV is processed by `enrich_analysy_add_category()` which: -
Normalises category labels to `SS_CATEGORY` (consistent across
backends) - Adds a `PHENOTYPE` column (`TRUE` if the term matches the
study name) - Applies multiple testing correction (BH or q-value) -
Sorts by adjusted p-value

### Common output columns

| Column | Description |
|----|----|
| `SS_CATEGORY` | Normalised category (`GO-BP`, `KEGG-ORA`, `KEGG-NTA`, `REACTOME-NTA`, `CHEMICAL`, …) |
| `ID` | Pathway / term identifier |
| `Description` | Human-readable term name |
| `SCORE` | Enrichment score (ORA: enrichment ratio; GSEA: NES; STRINGdb: FDR) |
| `PVALUE_ADJ_ALL_BH` | BH-adjusted p-value across all areas |
| `PHENOTYPE` | `TRUE` if term description matches study label |
| `SS_RANK` | Rank within `SS_CATEGORY` |

------------------------------------------------------------------------

## Cross-study and cross-subsample comparison

### `meta_enrichment_compare_studies()`

Compares enriched pathways across two or more independent cohorts and
produces Venn diagrams (Chen and Boutros 2011) of overlapping terms:

``` r

meta_enrichment_compare_studies(
  result_folder = "~/meta_results/",
  # studies A and B must already have Pathway/ results
)
```

### `meta_enrichment_overlaps_subsamples()`

Tests the **stability** of pathway enrichment under random subsampling
of the cohort. Run this to assess how robust the enriched pathways are
to sample variability:

``` r

meta_enrichment_overlaps_subsamples(
  inference_details  = inference_details,
  result_folder      = "~/semseeker_results/",
  pathway_package    = "WebGestalt",            # "" = compare all backends
  significance       = TRUE,
  pathways_sql_selection = "",                  # SQL to filter specific pathways
  old_label_samples_sql_condition    = "",
  new_label_samples_sql_condition    = "",
  association_pvalue_column          = "PVALUE_ADJ_ALL_BH"
)
```

------------------------------------------------------------------------

## Choosing the right backend

| Scenario | Recommended backend |
|----|----|
| First exploratory run, GO + KEGG | `enrich_WebGestalt()` with ORA |
| Small gene lists (\< 20 genes) | `enrich_WebGestalt()` with ORA |
| Ranked gene lists, pathway activity | `enrich_WebGestalt()` with GSEA |
| Protein interaction context | `enrich_STRINGdb()` |
| Partial pathway activation (sparse hits) | `enrich_pathfindR()` |
| Environmental/chemical exposure study | `enrich_ctdR()` |
| Known disease hypothesis | `enrich_Phenolyzer_WebGestalt()` or `enrich_Phenolyzer_STRINGdb()` |
| Replication across cohorts | [`meta_enrichment_compare_studies()`](https://corsaro-lab.github.io/SEMseeker/reference/meta_enrichment_compare_studies.md) |
| Bootstrap stability check | [`meta_enrichment_overlaps_subsamples()`](https://corsaro-lab.github.io/SEMseeker/reference/meta_enrichment_overlaps_subsamples.md) |

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

------------------------------------------------------------------------

## References

Ashburner, Michael et al. 2000. “Gene Ontology: Tool for the Unification
of Biology.” *Nature Genetics* 25: 25–29.
<https://doi.org/10.1038/75556>.

Boyle, Matthew R. et al. 2004. “GO:: TermFinder–Open Source Software for
Accessing Gene Ontology Information and Finding Significantly Enriched
Gene Ontology Terms Associated with a List of Genes.” *Bioinformatics*
20 (18): 3710–15. <https://doi.org/10.1093/bioinformatics/bth456>.

Chen, Hanbo, and Paul C. Boutros. 2011. “VennDiagram: A Package for the
Generation of Highly-Customizable Venn and Euler Diagrams in R.” *BMC
Bioinformatics* 12: 35. <https://doi.org/10.1186/1471-2105-12-35>.

Davis, Allan Peter et al. 2023. “Comparative Toxicogenomics Database
(CTD): Update 2023.” *Nucleic Acids Research* 51 (D1): D1257–62.
<https://doi.org/10.1093/nar/gkac833>.

Jassal, Bijay et al. 2020. “The Reactome Pathway Knowledgebase.”
*Nucleic Acids Research* 48 (D1): D498–503.
<https://doi.org/10.1093/nar/gkz1031>.

Kanehisa, Minoru, and Susumu Goto. 2000. “KEGG: Kyoto Encyclopedia of
Genes and Genomes.” *Nucleic Acids Research* 28 (1): 27–30.
<https://doi.org/10.1093/nar/28.1.27>.

Liao, Yuxing, Jing Wang, Eric R. Jaehnig, Zhiao Shi, and Bing Zhang.
2019. “WebGestalt 2019: Gene Set Analysis Toolkit with Revamped UIs and
APIs.” *Nucleic Acids Research* 47 (W1): W199–205.
<https://doi.org/10.1093/nar/gkz401>.

Subramanian, Aravind et al. 2005. “Gene Set Enrichment Analysis: A
Knowledge-Based Approach for Interpreting Genome-Wide Expression
Profiles.” *Proceedings of the National Academy of Sciences* 102 (43):
15545–50. <https://doi.org/10.1073/pnas.0506580102>.

Szklarczyk, Damian et al. 2019. “STRING V11: Protein–Protein Association
Networks with Increased Coverage.” *Nucleic Acids Research* 47 (D1):
D607–13. <https://doi.org/10.1093/nar/gky1131>.

Szklarczyk, Damian et al. 2023. “STRING V12.0: A Protein-Protein
Association Network with Increased Coverage.” *Nucleic Acids Research*
51 (D1): D638–46. <https://doi.org/10.1093/nar/gkac1000>.

Ulgen, Ege, Ozan Ozisik, and Ugur Sezerman. 2019. “pathfindR: An R
Package for Comprehensive Identification of Enriched Pathways in Omics
Data Through Active Subnetworks.” *Frontiers in Genetics* 10.
<https://doi.org/10.3389/fgene.2019.00858>.

Wang, Jing, Duncan Vasquez, Zhiao Shi, and Bing Zhang. 2017. “WEB-based
GEne SeT AnaLysis Toolkit (WebGestalt): Update 2017.” *Nucleic Acids
Research* 45 (W1): W130–37. <https://doi.org/10.1093/nar/gkx356>.
