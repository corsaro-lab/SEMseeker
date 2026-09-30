#' dmr_annotation
#'
#' Differentially methylated region (DMR) annotations for CpG probes on
#' Illumina methylation arrays.  Contains only probes that overlap at least one
#' known imprinted or disease-associated DMR.  This lightweight table
#' (~1,600 rows) supplements the Bioconductor array annotation packages used by
#' \code{\link{anno_probe_annotation_build}}, which do not carry DMR-level
#' annotations.
#'
#' @format A data frame with three columns:
#' \describe{
#'   \item{PROBE}{Illumina probe identifier (e.g. \code{"cg00000924"}).}
#'   \item{DMR_WHOLE}{Genomic region label of the enclosing DMR
#'     (e.g. \code{"KCNQ1OT1:TSS-DMR"}).}
#'   \item{DMR_DMR}{Fine-grained DMR label, identical to \code{DMR_WHOLE} for
#'     most probes.}
#' }
"dmr_annotation"


#' cytoband_hg19
#'
#' Cytogenetic band coordinates for the hg19 human genome assembly.
#' Used by \code{\link{anno_probe_annotation_build}} to assign a \code{CHR_CYTOBAND}
#' label (e.g. \code{"q12.2"}) to each CpG probe based on its chromosomal
#' position.  Band boundaries are derived from probe positions in the full
#' Illumina EPIC annotation and cover all autosomes and sex chromosomes.
#'
#' @format A data frame with 829 rows and four columns:
#' \describe{
#'   \item{CHR}{Chromosome identifier without \code{"chr"} prefix
#'     (e.g. \code{"1"}, \code{"X"}).}
#'   \item{START}{Approximate start position of the cytogenetic band (bp).}
#'   \item{END}{Approximate end position of the cytogenetic band (bp).}
#'   \item{CYTOBAND}{ISCN band label (e.g. \code{"q12.2"}, \code{"p36.33"}).}
#' }
"cytoband_hg19"


#' metrics_properties
#'
#' Metadata table describing the statistical properties of each SEM metric.
#' Used internally to determine ranking direction and scaling behaviour during
#' association analysis and pathway enrichment.
#'
#' @format A data frame with one row per metric and columns including
#'   \code{Metric} (metric name), \code{Higher_the_Better} (logical: whether
#'   higher values indicate stronger signal), and \code{Affected_by_Scaling}
#'   (logical: whether the metric is affected by data transformation).
"metrics_properties"


#' ssEnv
#'
#' Internal session environment object persisted between SEMseeker analysis
#' steps. Stores runtime parameters (result folder paths, technology flag,
#' alpha threshold, etc.) set by \code{\link{core_init_env}} and retrieved by
#' \code{core_get_session_info()}.
#'
#' @format An \code{environment} containing named slots for session-level
#'   analysis parameters. Users should not modify this object directly;
#'   use \code{\link{core_init_env}} and \code{core_set_env_variable()} instead.
"ssEnv"


#' test_master_features
#'
#' Master probe-features fixture used by the test suite and by the
#' \code{getting-started} vignette as a runnable synthetic example.
#' Contains 20,000 real Illumina EPIC (K850) probe IDs, sampled
#' deterministically around 58 curated human imprinting DMRs (KCNQ1OT1,
#' H19/IGF2, MEG3/DLK1, GNAS, PEG3, SNURF, PLAGL1, NNAT, MEST, ...).
#'
#' Construction strategy (see \code{data-raw/build_test_master_features.R}):
#' \enumerate{
#'   \item All unique probe IDs in \code{\link{dmr_annotation}} are kept as
#'     the biological \emph{signal} layer (~820 probes labelled with their
#'     parent DMR in \code{DMR_LABEL}).
#'   \item Each imprinting DMR is expanded into a genomic window of
#'     \eqn{\pm 500\,\textrm{kb}} to capture the surrounding imprinted gene
#'     cluster; all EPIC probes within these windows are added.
#'   \item If the resulting pool is below 20,000, the remainder is filled
#'     with EPIC probes sampled deterministically from outside any imprinting
#'     window (background layer, \code{DMR_LABEL = NA}).
#' }
#'
#' This makes the fixture biology-aware (BWS / SRS / PWS-AS / GNAS-related
#' regions are guaranteed to be represented) while preserving statistical
#' validity (20k probes is the minimum size used by the test suite for IQR
#' robustness).  Because the fixture ships with real EPIC probe IDs, the
#' Bioconductor annotation package is not required at test-setup time --
#' eliminating the macOS \code{tcltk}/XQuartz segfault path triggered by
#' \code{requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")}.
#'
#' @format A data frame with 20,000 rows and 6 columns:
#' \describe{
#'   \item{PROBE}{Illumina EPIC probe identifier (e.g. \code{"cg00000924"}).}
#'   \item{CHR}{Chromosome label without \code{"chr"} prefix (e.g. \code{"11"},
#'     \code{"X"}).}
#'   \item{START}{1-based start position on hg19 (bp).}
#'   \item{END}{End position; equal to \code{START} for array probes.}
#'   \item{ABSOLUTE}{Concatenated chromosome/position key
#'     (\code{paste(CHR, START, sep = "_")}).}
#'   \item{DMR_LABEL}{Imprinting DMR identifier (e.g.
#'     \code{"KCNQ1OT1:TSS-DMR"}) for the ~820 \emph{signal} probes;
#'     \code{NA} for flanking and background probes.}
#' }
#'
#' @source Built from \code{\link{dmr_annotation}} and the
#'   \code{IlluminaHumanMethylationEPICanno.ilm10b4.hg19} Bioconductor
#'   annotation package with \code{set.seed(20210713)} (date of the v.0.1.9
#'   Zenodo software-archive release, DOI
#'   \href{https://doi.org/10.5281/zenodo.5095416}{10.5281/zenodo.5095416}).
#'   Re-generate with \code{Rscript data-raw/build_test_master_features.R}.
"test_master_features"


#' test_signal_gse133774
#'
#' Real EPIC methylation beta-value matrix from GEO series GSE133774
#' (Infinium MethylationEPIC 850k, GPL21145), filtered to the probe IDs
#' present in \code{\link{test_master_features}}.  Contains 10 samples from
#' a Beckwith-Wiedemann Syndrome (BWS) / Multi-Locus Imprinting Disturbance
#' (MLID) family study: 6 unrelated controls and 4 family members (L1 = BWS
#' proband, L2-L4 = siblings/parent with NLRP5 compound heterozygous variants).
#'
#' This fixture replaces all \code{rbeta()}-based synthetic signal generation
#' in the test suite and vignette.  Because the data contain real BWS
#' epimutations at imprinting DMRs (KCNQ1OT1, H19/IGF2, MEG3, ...), running
#' the full SEMseeker pipeline on this matrix detects biologically expected
#' hypo-epimutation events in the Case samples without any artificial injection.
#'
#' @format A numeric matrix with 18,089 rows (EPIC probe IDs) and 10 columns
#'   (samples: CTRL01-CTRL06, L1-L4). Values are beta coefficients in
#'   \eqn{[0, 1]}; a small number of probes may have \code{NA} (QC-filtered
#'   positions in the original GEO submission).
#'
#' @source GEO accession GSE133774, series matrix file parsed by
#'   \code{data-raw/build_test_signal_fixture.R}.  Original study: Docherty
#'   \emph{et al.} (2020), NLRP5 variants associated with MLID and BWS.
"test_signal_gse133774"


#' test_samplesheet_gse133774
#'
#' Sample sheet for \code{\link{test_signal_gse133774}} in the canonical
#' SEMseeker three-class design (Reference / Control / Case).
#'
#' Control samples (CTRL01-CTRL06) appear twice: once as \code{Reference}
#' (population baseline for IQR threshold estimation) and once as
#' \code{Control} (comparison group).  Family samples (L1-L4) are \code{Case};
#' L1 is the BWS proband.  This is the \emph{Reference-reuse pattern}
#' documented in the getting-started vignette.
#'
#' @format A data frame with 16 rows and 2 columns:
#' \describe{
#'   \item{Sample_ID}{Sample identifier matching column names of
#'     \code{\link{test_signal_gse133774}} (e.g. \code{"CTRL01"}, \code{"L1"}).}
#'   \item{Sample_Group}{One of \code{"Reference"}, \code{"Control"},
#'     or \code{"Case"}.}
#' }
#'
#' @source Derived from GEO accession GSE133774 metadata.
"test_samplesheet_gse133774"
