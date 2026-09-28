#' Compose the identity key of a computed artefact (internal)
#'
#' AI-255. **The** single compositor of the six-coordinate key. Everything that
#' names an artefact — a pivot file, a row of an inference result, an entry in a
#' resume or overlap comparison — goes through here, so two artefacts are the
#' same thing exactly when their keys are equal.
#'
#' \preformatted{
#' <MARKER>_<FIGURE>_<SCOPE>_<AREA>_<SUBAREA>_<AGGREGATION>_<GENOME_BUILD>
#'
#' MUTATIONS_HYPER_SAMPLE_PROBE_WHOLE_SUM_HG19      burden of the whole sample
#' MUTATIONS_HYPER_SAMPLE_GENE_WHOLE_SUM_HG19       same, over gene probes only
#' DELTAS_HYPO_INSTANCE_GENE_TSS200_MEDIAN_HG19     median per gene, TSS200 window
#' SIGNAL_BETA_INSTANCE_PROBE_WHOLE_VALUE_HG19      the beta value per probe
#' }
#'
#' Before AI-255 identity was the *combination* of columns, checked separately in
#' deduplication, in the resume match and in the cross-study overlaps. Every new
#' coordinate had to be remembered in each of those places — the work AI-248 did
#' by hand across three files — and forgetting one is exactly how a silent defect
#' is born. One composed string means the seventh coordinate will change one
#' function instead of five call sites.
#'
#' **The key is never parsed by position.** `core_name_cleaning()` maps every
#' non-alphanumeric character to `_`, so no separator hierarchy is available, and
#' two vocabularies carry an underscore inside their values: `SUBAREA` has
#' `N_SHORE` / `S_SHELF` (aligned to Illumina's `Relation_to_Island`, so not ours
#' to rename). `ISLAND_N_SHORE` is three tokens for two coordinates. Coordinates
#' travel as columns; if one ever has to be recovered from a key, it is matched
#' against the closed vocabularies, longest first — never by index.
#'
#' The genome build is part of the identity, not decoration: the same gene on
#' hg19 and on hg38 is not the same thing.
#'
#' @param marker,figure the quantity and its side (or, for `SIGNAL`, its scale).
#' @param scope `"SAMPLE"` (one number per sample) or `"INSTANCE"` (one number
#'   per instance of the region class).
#' @param area,subarea the region class.
#' @param aggregation how the positions are reduced; `"VALUE"` means the block is
#'   already a single position and nothing is reduced.
#' @param genome_build defaults to the session's build.
#' @return the key, uppercase via [core_name_cleaning()].
#' @keywords internal
#' @noRd
io_artefact_key <- function(marker, figure, scope, area, subarea, aggregation,
                            genome_build = NULL) {

  scope <- io_scope_validate(scope)

  if (is.null(genome_build) || !nzchar(as.character(genome_build))) {
    ssEnv <- core_get_session_info()
    genome_build <- if (!is.null(ssEnv$genome_build) && nzchar(ssEnv$genome_build))
      ssEnv$genome_build else "hg19"
  }

  parts <- c(marker, figure, scope, area, subarea, aggregation, genome_build)
  if (any(vapply(parts, function(p) is.null(p) || !nzchar(as.character(p)[1]),
                 logical(1))))
    stop("io_artefact_key(): every coordinate is required — marker, figure, ",
         "scope, area, subarea, aggregation. An artefact whose name omits one ",
         "cannot be told apart from an artefact that omits a different one.",
         call. = FALSE)

  core_name_cleaning(paste(as.character(parts), collapse = "_"))
}

#' The two legal values of the SCOPE coordinate (internal)
#'
#' AI-255. `SAMPLE` is the partition with one block — the whole sample reduced to
#' one number — and `INSTANCE` is the partition induced by the region class, one
#' number per gene, island, cytoband or probe. There is no third value: the
#' historical 1/2/3 depth scale tried to order a lattice on a line and lost the
#' pairs that are not comparable (a TSS200 window and an open-sea stretch refine
#' neither each other nor anything in between).
#'
#' @keywords internal
#' @noRd
io_scope_vocabulary <- function() c("SAMPLE", "INSTANCE")

#' @rdname io_scope_vocabulary
#' @keywords internal
#' @noRd
io_scope_validate <- function(scope) {
  if (is.null(scope) || length(scope) != 1L || !nzchar(as.character(scope)))
    stop("scope is required: \"SAMPLE\" (one number per sample) or ",
         "\"INSTANCE\" (one number per instance of the region class).",
         call. = FALSE)
  scope <- toupper(as.character(scope))
  if (!(scope %in% io_scope_vocabulary()))
    stop("scope = '", scope, "' is not a scope. Legal values: ",
         paste(io_scope_vocabulary(), collapse = ", "), ".", call. = FALSE)
  scope
}

#' Is this region class a single position? (internal)
#'
#' AI-255. `PROBE` and `POSITION` are the bottom of the refinement lattice: one
#' block per position. It is the only place where `VALUE` — the identity, "I do
#' not aggregate" — is a meaningful aggregation, and the only place where an
#' omitted aggregation can be resolved without guessing.
#'
#' @keywords internal
#' @noRd
io_area_is_single_position <- function(area) {
  toupper(as.character(area)) %in% c("PROBE", "POSITION")
}

#' Which single-position area this technology speaks (internal)
#'
#' AI-308. `PROBE` and `POSITION` are two names for the same bottom of the
#' lattice, and each technology has exactly one that means anything: Illumina
#' reports probe ids (`cg00000029`), long reads have no probe concept and are
#' keyed by coordinate. `assoc_run_marker()` already skips the other one, the
#' symmetric guard added by AI-098, so the two never both reach a result.
#'
#' At `SCOPE = SAMPLE` that guard is not enough on its own. The registry always
#' carries `POSITION` (`util_keys_create()` forces it in), so an Illumina run
#' that did not declare `PROBE` among its areas would have its whole-sample
#' burden built on `POSITION_WHOLE` and then dropped by the guard: a row
#' missing from the result, which is the failure mode that looks like an answer.
#' Naming the canonical area lets the collapsed branch normalise to it instead.
#'
#' @param tech the session technology; read from the session when omitted.
#' @return `"POSITION"` for WGBS/LONGREAD, `"PROBE"` otherwise.
#' @keywords internal
#' @noRd
io_single_position_area <- function(tech = NULL) {
  if (is.null(tech)) {
    ssEnv <- core_get_session_info()
    tech <- ssEnv$tech
  }
  if (!is.null(tech) && toupper(as.character(tech)[1]) %in% c("WGBS", "LONGREAD"))
    "POSITION" else "PROBE"
}
