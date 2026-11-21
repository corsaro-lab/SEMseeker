#' Run SEMseeker on methylation data from any supported source
#'
#' Public entry point. Accepts a wide range of inputs - bedmethyl files
#' (modkit / nanopolish), coordinate-based data frames (WGBS / long-read),
#' Illumina probe-indexed matrices, or already-loaded data frames - normalises
#' them to the internal format, fixes the assembly of the run (see
#' \code{genome_build}),
#' and delegates to the core pipeline \code{\link{sem_core}}. Input
#' values are passed through unchanged: if you need to convert M-values to
#' beta, call \code{\link{sem_mvalue_to_beta}} explicitly before \code{semseeker()}.
#'
#' Supported \code{input} forms:
#' \itemize{
#'   \item Character vector of bedmethyl file paths (\code{.bed}/\code{.tsv}/
#'     \code{.bedmethyl}) - parsed via \code{\link{io_bedmethyl_read}}.
#'   \item Data frame with \code{CHR}/\code{START}[\code{/END}] columns
#'     (WGBS / long-read coordinate format) - normalised via
#'     \code{\link{io_normalize_signal_input}}.
#'   \item Matrix or data frame with probe-ID rownames (Illumina array) -
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
#' @param genome_build Reference assembly of the coordinates: \code{"hg19"},
#'   \code{"hg38"}, \code{"mm10"} or \code{"legacy"}. A run works on one
#'   assembly, which names every pivot and is recorded in the session. For
#'   coordinate inputs (bedmethyl files, WGBS or long-read data frames) only
#'   whoever aligned the reads knows it, so it has no default and must be given.
#'   For an Illumina matrix it is the assembly of the array manifest (hg19 for
#'   the 27k, 450k and EPIC manifests) and can be omitted; a different value is
#'   refused, because the probe coordinates would keep the manifest's assembly
#'   while the areas were built on the declared one. A result folder keeps the
#'   assembly of its first run for good: a later run on the same folder that
#'   declares another one is refused, \code{start_fresh = TRUE} included, and
#'   one that declares none uses the recorded one. SEMseeker never converts
#'   coordinates from one assembly to another. Every refusal comes before
#'   anything is written.
#' @param ... Additional arguments forwarded to \code{\link{sem_core}}
#'   and \code{core_init_env()} (e.g. \code{parallel_strategy}, \code{alpha},
#'   \code{LESIONS_BP}, \code{marker}, \code{areas}). \code{LESIONS_BP}
#'   (default 2000) is the maximum bp distance between two probes for them to
#'   be in the same LESIONS enrichment window - replaces the legacy
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
#'   tech = "WGBS",
#'   genome_build = "hg38"
#' )
#' }
#'
#' @export
semseeker <- function(input,
                      sample_sheet,
                      result_folder,
                      input_type = c("auto", "bedmethyl", "coord_df", "matrix"),
                      tech = NULL,
                      genome_build = NULL,
                      ...) {

  input_type <- match.arg(input_type)

  # ---- Step 0: the assembly of the run, refused before any folder exists --
  genome_build <- .sem_genome_build_resolve(input, input_type, tech, genome_build, result_folder)
  .sem_regions_check(list(...)[["areas"]], list(...)[["subareas"]])

  # ---- Step 1: core_init_env FIRST - cleans folder (start_fresh), sets up
  #      parallel plan, configures session. Must happen before any
  #      data processing or file I/O. -----------------------------------
  core_init_env(
    result_folder = result_folder,
    tech          = if (is.null(tech)) "" else tech,
    genome_build  = genome_build,
    ...
  )

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

#' Fix the one assembly a run works on
#'
#' Every position SEMseeker writes is a coordinate, so the assembly is part of
#' its identity. An Illumina manifest carries its own assembly; a coordinate
#' input carries none that SEMseeker could read, and a silent default would
#' label the coordinates of another assembly.
#'
#' @param input the \code{input} of \code{semseeker()}: one source or a list of
#'   sources, one per batch.
#' @param input_type the \code{input_type} of \code{semseeker()}, already matched;
#'   \code{"auto"} is resolved here per source.
#' @param tech the declared technology, or \code{NULL} when it is left to
#'   auto-detection.
#' @param genome_build the assembly the user declared, or \code{NULL}.
#' @param result_folder the result folder of the run; when it already holds a
#'   run, its recorded assembly stands for an undeclared one, and
#'   core_init_env() refuses a different one.
#' @return The assembly of the run, a single string.
#' @keywords internal
.sem_genome_build_resolve <- function(input, input_type, tech, genome_build, result_folder = NULL) {
  if (is.null(genome_build))
    genome_build <- core_result_folder_genome_build(result_folder)
  allowed_builds <- .SS_DEFAULTS$genome_build$choices
  if (!is.null(genome_build) && !(length(genome_build) == 1L && genome_build %in% allowed_builds))
    stop("genome_build must be one of ", paste0("'", allowed_builds, "'", collapse = ", "), "; got '",
         paste(genome_build, collapse = "', '"), "'.", call. = FALSE)

  sources <- if (is.list(input) && !is.data.frame(input)) input else list(input)
  source_kinds <- vapply(sources, .sem_source_kind, character(1), input_type = input_type, tech = tech)

  if (any(source_kinds == "coordinates") && is.null(genome_build))
    stop("genome_build is required for coordinate input (bedmethyl, WGBS, long reads): declare the assembly ",
         "the reads were aligned to, e.g. genome_build = 'hg38'.", call. = FALSE)

  if (any(source_kinds == "illumina")) {
    manifests <- if (!is.null(tech) && tech %in% names(.ANNO_PKGS)) .ANNO_PKGS[[tech]] else unname(.ANNO_PKGS)
    manifest_build <- unique(.anno_pkg_genome_build(manifests))
    if (length(manifest_build) != 1L)
      stop("the Illumina manifests disagree on the assembly; declare tech so that one manifest is used.", call. = FALSE)
    if (!is.null(genome_build) && !identical(genome_build, manifest_build))
      stop("genome_build = '", genome_build, "' does not match the Illumina manifest (",
           paste(manifests, collapse = ", "), "), whose coordinates are ", manifest_build,
           ". Omit genome_build, or lift the data over to ", genome_build,
           " and pass it as coordinates.", call. = FALSE)
    genome_build <- manifest_build
  }
  genome_build
}

#' Refuse region classes SEMseeker does not build
#'
#' Checked before any folder exists: an unknown area or subarea used to reach
#' the end of core_init_env() as an unrecognised argument, after the folder and
#' its session had been created.
#'
#' @param areas the \code{areas} passed to \code{semseeker()}, or \code{NULL}
#'   for the default set; \code{"ALL"} asks for every area.
#' @param subareas the \code{subareas} passed to \code{semseeker()}, or
#'   \code{NULL}; \code{"ALL"} asks for every subarea.
#' @return Invisibly \code{NULL}; stops naming the unknown values and the
#'   values allowed.
#' @keywords internal
.sem_regions_check <- function(areas, subareas) {
  allowed <- util_areas_subareas_allowed()
  allowed_areas <- setdiff(unique(allowed$AREA), "POSITION")
  unknown_areas <- setdiff(areas, c(allowed_areas, "POSITION", "ALL"))
  if (length(unknown_areas) > 0L)
    stop("areas: SEMseeker builds no area ", paste0("'", unknown_areas, "'", collapse = ", "), ". Allowed: ",
         paste(allowed_areas, collapse = ", "), ", or 'ALL'.", call. = FALSE)
  allowed_subareas <- unique(allowed$SUBAREA)
  unknown_subareas <- setdiff(subareas, c(allowed_subareas, "ALL"))
  if (length(unknown_subareas) > 0L)
    stop("subareas: SEMseeker builds no subarea ", paste0("'", unknown_subareas, "'", collapse = ", "), ". Allowed: ",
         paste(allowed_subareas, collapse = ", "), ", or 'ALL'.", call. = FALSE)
  invisible(NULL)
}

#' Tell an Illumina source from a coordinate source
#'
#' @param one_input one source of \code{semseeker()}: file paths, a data frame
#'   or a matrix.
#' @param input_type as in \code{.sem_genome_build_resolve()}.
#' @param tech as in \code{.sem_genome_build_resolve()}.
#' @return \code{"illumina"} or \code{"coordinates"}.
#' @keywords internal
.sem_source_kind <- function(one_input, input_type, tech) {
  if (!is.null(tech) && nzchar(tech))
    return(if (tech %in% names(.ANNO_PKGS)) "illumina" else "coordinates")
  source_type <- if (identical(input_type, "auto")) .sem_detect_input_type(one_input) else input_type
  if (!identical(source_type, "matrix"))
    return("coordinates")
  # a matrix keyed by CHR_START is coordinate data, as core_get_meth_tech() reads it
  row_ids <- rownames(one_input)
  if (length(row_ids) > 0L && all(grepl("_", row_ids[seq_len(min(1000L, length(row_ids)))])))
    "coordinates" else "illumina"
}
