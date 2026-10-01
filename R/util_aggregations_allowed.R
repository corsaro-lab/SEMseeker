#' Aggregations admissible for an artefact (internal)
#'
#' AI-248, extended by AI-255. Single source of truth of the AGGREGATION axis:
#' which ways of reducing a set of positions to one number are admissible for a
#' given artefact. The producer iterates over this to know what to compute, the
#' consumer validates against it, and [io_artefact_key()] names the result - so
#' the computed keys and the declared vocabulary cannot drift apart.
#'
#' The axis is **permissive by design**: every aggregation applies to every
#' marker. Whether the median of hypermethylated mutations *means* something is a
#' question for the biologist reading the table; making it expressible and
#' nameable in a uniform way is the package's job. This generalises the old
#' binary `DISCRETE` flag, which could only say "sum" or "mean".
#'
#' Three restrictions, none of them arbitrary:
#'
#' 1. **`VALUE` is the only operator on a single position.** At `PROBE` or
#'    `POSITION` granularity with `SCOPE = INSTANCE` the block holds one value:
#'    its sum, mean and median are all that same value, so naming them would
#'    invite the reader to believe a reduction happened.
#' 2. **`MODELOW` / `MODEHIGH` need `SCOPE = SAMPLE`.** They estimate the two
#'    peaks of a bimodal density, which needs the whole distribution *and* enough
#'    of it. Per instance there are ~19 probes in a gene and often 2-5 in a
#'    TSS200 window: the estimate would be dominated by the bandwidth and would
#'    return a plausible-looking number instead of refusing - see
#'    [util_signal_descriptors()], which now also guards on numerosity.
#' 3. **The two modes exist only on the bounded scale**, `SIGNAL` with figure
#'    `BETA`. On M-values the distribution is unimodal and the split around 0.5
#'    carries no meaning.
#'
#' @param marker,figure the pair. `figure` matters only for `SIGNAL`, where it
#'   carries the scale (see [io_signal_figure()]).
#' @param discrete the `DISCRETE` flag of the key. It no longer decides the
#'   operator - that is the point of this axis - but it still says which
#'   aggregation is produced **by default** for a marker that is not `SIGNAL`:
#'   a count is summed, a continuous deviation is averaged. Preserving that
#'   distinction is what keeps the migration from silently changing the numbers
#'   of the DELTA family.
#' @param default when `TRUE` (default) return only what is produced without
#'   being asked. When `FALSE` return everything the taxonomy admits.
#' @param scope `"SAMPLE"` or `"INSTANCE"` (see [io_scope_vocabulary()]).
#' @param area the region class; only needed to recognise the single-position
#'   case of restriction 1.
#' @return character vector of aggregation names.
#' @keywords internal
#' @noRd
util_aggregations_allowed <- function(marker, figure = NULL, discrete = TRUE,
                                      default = TRUE, scope = "SAMPLE",
                                      area = NULL) {

  marker <- toupper(as.character(marker))
  figure <- if (is.null(figure)) "" else toupper(as.character(figure))
  scope  <- io_scope_validate(scope)

  # Restriction 1: one position in the block, one operator that means anything.
  if (identical(scope, "INSTANCE") && !is.null(area) &&
      io_area_is_single_position(area))
    return("VALUE")

  # Two classes, and the class is what decides. SIGNAL, DELTAS and DELTAR are
  # all continuous - a value per position - and take the same operators; the
  # instability markers are counts, 0 or 1 per position.
  # SIGNAL is continuous by definition, whatever the caller passes: the raw
  # signal is never a count.
  if (identical(marker, "SIGNAL"))
    discrete <- FALSE

  # Restrictions 2 and 3 together: bounded scale AND one big group.
  bounded <- identical(marker, "SIGNAL") && identical(figure, "BETA") &&
             identical(scope, "SAMPLE")

  if (isTRUE(discrete)) {
    # SUM is the burden. MEAN of a 0/1 vector is the same burden divided by the
    # positions of the scope, i.e. the DENSITY of epimutations - the only form
    # comparable across regions of different size, which a raw burden is not.
    # MEDIAN, VARIANCE and IQR of a 0/1 vector are degenerate, so they are not
    # admissible: naming them would only add noise to the file.
    admissible <- c("SUM", "MEAN")
    return(if (isTRUE(default)) "SUM" else admissible)
  }

  admissible <- c("SUM", "MEAN", "MEDIAN", "VARIANCE", "IQR")
  if (bounded)
    admissible <- c(admissible, "MODELOW", "MODEHIGH")

  if (!isTRUE(default))
    return(admissible)

  # Every continuous marker is produced with the whole descriptor set: no
  # operator is privileged for a continuous value. This ADDS columns without
  # changing any existing number - the historical mean keeps its value and
  # merely gains the _MEAN suffix.
  setdiff(admissible, "SUM")
}
