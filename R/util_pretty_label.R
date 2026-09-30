#' Convert a snake_case / SCREAMING_SNAKE token into a display-friendly label.
#'
#' SEMseeker uses `core_name_cleaning()` on disk artifacts (file names, CSV
#' column names, AREA / MARKER / FIGURE identifiers): UPPERCASE with
#' `_` as the only legal separator. That is the right convention for
#' stable IDs but it reads poorly inside plot titles, axes, and legend
#' labels, where a human is meant to skim the chart at a glance.
#'
#' `util_pretty_label()` reverses just the cosmetic side of that contract:
#'
#'   1. Replace every `_` with a single space.
#'   2. Collapse runs of two or more whitespace characters into one.
#'   3. Trim leading / trailing whitespace.
#'
#' It does NOT lowercase, title-case, or replace any other character -
#' that keeps the function fully reversible visually (a reader can
#' mentally re-insert the underscores to find the matching CSV column).
#'
#' Use this in any plotting function (volcano, box, fitted-model, manhattan)
#' that consumes SEMseeker's canonical UPPER_WITH_UNDERSCORE identifiers.
#' Three current inline call sites should migrate to this helper to keep
#' a single source of truth for the display contract:
#'
#'   - `R/association_cross_studies_overlaps.R:221`
#'   - `R/association_cross_subsamples_overlaps.R:227, :287`
#'
#' @param x A character scalar or vector.
#' @return Character vector of the same length with the three substitutions
#'   applied. `NA` and empty strings pass through unchanged.
#'
#' @examples
#' \dontrun{
#' util_pretty_label("TUMOUR_STAGE_N")
#' #> "TUMOUR STAGE N"
#' util_pretty_label("PVALUE_ADJ_ALL_FDR")
#' #> "PVALUE ADJ ALL FDR"
#' util_pretty_label(c("DELTARP_GENE_TSS200", NA, ""))
#' #> "DELTARP GENE TSS200" NA                  ""
#' }
#'
#' @keywords internal
#' @noRd
util_pretty_label <- function(x) {
  if (length(x) == 0L) return(x)
  out <- ifelse(is.na(x) | !nzchar(x), x,
                trimws(gsub("\\s{2,}", " ",
                             gsub("_", " ", x, fixed = TRUE))))
  out
}
