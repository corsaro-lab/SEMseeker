#' Run SEMseeker on methylation data from any supported source
#'
#' Public entry point. Accepts a wide range of inputs — bedmethyl files
#' (modkit / nanopolish), coordinate-based data frames (WGBS / long-read),
#' Illumina probe-indexed matrices, or already-loaded data frames — normalises
#' them to the internal format, validates the tech / genome_build combination,
#' and delegates to the core pipeline \code{\link{sem_core}}. Input
#' values are passed through unchanged: if you need to convert M-values to
#' beta, call \code{\link{sem_mvalue_to_beta}} explicitly before \code{semseeker()}.
#'
#' Supported \code{input} forms:
#' \itemize{
#'   \item Character vector of bedmethyl file paths (\code{.bed}/\code{.tsv}/
#'     \code{.bedmethyl}) — parsed via \code{\link{io_bedmethyl_read}}.
#'   \item Data frame with \code{CHR}/\code{START}[\code{/END}] columns
#'     (WGBS / long-read coordinate format) — normalised via
#'     \code{\link{io_normalize_signal_input}}.
#'   \item Matrix or data frame with probe-ID rownames (Illumina array) —
#'     passed through unchanged.
#'   \item List of any of the above, one element per batch.
#' }
#'
#' @param input Input signal data. See details for supported forms.
#' @param sample_sheet Data frame (or list of data frames) with a
#'   \code{Sample_ID} column identifying samples.
#' @param result_folder Output directory.
#' @param input_type One of \code{"auto"} (default), \code{"bedmethyl"},
#'   \code{"coord_df"}, \code{"matrix"}. When \code{"auto"}, the type is
#'   inferred from the object class / file extension.
#' @param tech Optional technology label (\code{"K850"}, \code{"K450"},
#'   \code{"K27"}, \code{"WGBS"}, \code{"LONGREAD"}). If \code{NULL} (default)
#'   it is auto-detected downstream by \code{core_get_meth_tech()}.
#' @param genome_build Reference genome build. One of \code{"hg19"} (default),
#'   \code{"hg38"}, \code{"mm10"}, \code{"legacy"}.
#' @param strict_build_check If \code{TRUE} (default), impossible tech /
#'   genome_build combinations (e.g. \code{tech="LONGREAD"} with
#'   \code{genome_build="hg19"}) raise an error. If \code{FALSE}, a warning is
#'   emitted and the pipeline proceeds.
#' @param ... Additional arguments forwarded to \code{\link{sem_core}}
#'   and \code{core_init_env()} (e.g. \code{parallel_strategy}, \code{alpha},
#'   \code{LESIONS_BP}, \code{marker}, \code{areas}). \code{LESIONS_BP}
#'   (default 2000) is the maximum bp distance between two probes for them to
#'   be in the same LESIONS enrichment window — replaces the legacy
#'   \code{sliding_window_size} probe-count parameter (removed in 0.99.2).
#'   \code{coverage_minimum} (default 80) is the minimum percentage of input
#'   positions that must be present in the reference annotation: the coverage
#'   charts are written on every run, and the run stops below this threshold
#'   Lower it explicitly for deliberate cross-technology runs.
#'   The per-sample statistics are no longer selected here: they are
#'   built on demand by \code{sem_study_summary_get(regions = ...)} and by
#'   \code{association_analysis(scope = "SAMPLE")} over the region classes that
#'   call declares, so a region class that was not foreseen at run time costs one
#'   scan instead of a whole rerun. Each position is counted once, even when it is annotated to
#'   several genes.
#'
#' @return Invisibly \code{NULL}; writes output files to \code{result_folder}.
#'
#' @examples
#' # Stub: see vignette('imprinting-disorders', package = 'SEMseeker') for a
#' # runnable Beckwith-Wiedemann workflow on the GSE133774 subset.
#' invisible(NULL)
#' \dontrun{
#' # Bedmethyl (Nanopore / modkit):
#' semseeker(
#'   input = list.files("bedmethyl/", pattern = "\\.bed$", full.names = TRUE),
#'   sample_sheet = ss,
#'   result_folder = tempdir(),
#'   tech = "LONGREAD",
#'   genome_build = "hg38"
#' )
#'
#' # Illumina beta matrix:
#' semseeker(
#'   input = beta_matrix,
#'   sample_sheet = ss,
#'   result_folder = tempdir()
#' )
#'
#' # WGBS coordinate data frame:
#' semseeker(
#'   input = wgbs_df,              # columns: CHR, START, END, sample1, ...
#'   sample_sheet = ss,
#'   result_folder = tempdir(),
#'   tech = "WGBS"
#' )
#' }
#'
#' @export
semseeker <- function(input,
                      sample_sheet,
                      result_folder,
                      input_type = c("auto", "bedmethyl", "coord_df", "matrix"),
                      tech = NULL,
                      genome_build = "hg19",
                      strict_build_check = TRUE,
                      ...) {

  input_type <- match.arg(input_type)

  # ---- Step 1: core_init_env FIRST — cleans folder (start_fresh), sets up
  #      parallel plan, configures session. Must happen before any
  #      data processing or file I/O. -----------------------------------
  core_init_env(
    result_folder = result_folder,
    tech          = if (is.null(tech)) "" else tech,
    genome_build  = genome_build,
    ...
  )

  # ---- Step 2: validate tech × genome_build ----------------------------
  .sem_validate_tech_build(tech, genome_build, strict_build_check)

  # ---- Step 3: normalise input to SEMseeker signal_data ----------------
  # Values are passed through unchanged. The beta-vs-M-value flag is
  # detected and stored downstream by core_get_meth_tech() (ssEnv$beta).
  signal_data <- .sem_dispatch_one(input, input_type)

  # ---- Step 4: delegate to core pipeline -------------------------------
  sem_core(
    sample_sheet  = sample_sheet,
    signal_data   = signal_data,
    result_folder = result_folder
  )
}

# --- Internal helpers --------------------------------------------------------

#' @keywords internal
.sem_dispatch_one <- function(input, input_type) {

  # Auto-detect when requested
  if (identical(input_type, "auto")) {
    input_type <- .sem_detect_input_type(input)
  }

  signal_data <- switch(
    input_type,
    bedmethyl = io_bedmethyl_read(input),
    coord_df  = io_normalize_signal_input(input),
    matrix    = input,
    stop(".sem_dispatch_one(): unknown input_type '", input_type, "'.")
  )

  # coord_df path: io_normalize_signal_input() returns probe-ID-indexed df already
  # bedmethyl path: returns coord df → also run through io_normalize_signal_input
  if (identical(input_type, "bedmethyl")) {
    signal_data <- io_normalize_signal_input(signal_data)
  }

  signal_data
}

#' @keywords internal
.sem_detect_input_type <- function(input) {
  # List of inputs → apply per element; dispatcher delegates one-by-one
  # (sem_core itself iterates lists, so we just need each element
  # normalised before being wrapped again.)
  if (is.character(input)) {
    exts <- tolower(tools::file_ext(input))
    if (all(exts %in% c("bed", "tsv", "bedmethyl", "gz")))
      return("bedmethyl")
    stop(".sem_detect_input_type(): character input with unsupported extensions: ",
         paste(unique(exts), collapse = ", "))
  }
  if (is.data.frame(input) && io_is_coord_format(input))
    return("coord_df")
  if (is.matrix(input) || is.data.frame(input))
    return("matrix")
  stop(".sem_detect_input_type(): cannot infer input_type for object of class '",
       class(input)[1L], "'. Pass input_type explicitly.")
}

#' @keywords internal
.sem_validate_tech_build <- function(tech, genome_build, strict) {
  if (is.null(tech) || !nzchar(tech)) return(invisible(NULL))

  bad <- (identical(tech, "LONGREAD") && identical(genome_build, "hg19"))

  if (bad) {
    msg <- paste0(
      "tech = 'LONGREAD' with genome_build = 'hg19' is almost certainly wrong: ",
      "long-read methylation pipelines (Nanopore / PacBio) align to GRCh38. ",
      "Set genome_build = 'hg38' in core_init_env() / semseeker()."
    )
    if (isTRUE(strict)) stop(msg) else warning(msg)
  }

  invisible(NULL)
}
