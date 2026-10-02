# Internal utilities for accessing Illumina methylation annotation packages.
#
# The IlluminaHumanMethylation*anno packages expose their data as *separate*
# lazy data objects ("Locations", "Islands.UCSC", "Other") that must be loaded
# individually with data().  The top-level S4 wrapper object (which has the
# same name as the package) only contains lazy-loading descriptors in its @data
# slot - accessing slot(obj, "data")$Locations directly returns a list of the
# form list(what="Locations", envir="package:...") rather than the actual table.
#
# Correct access pattern (no minfi required):
#   env <- new.env(parent = emptyenv())
#   utils::data("Locations", package = pkg, envir = env)
#   locs <- as.data.frame(get("Locations", envir = env), stringsAsFactors = FALSE)

# Map from SEMseeker technology key to Bioconductor annotation package name
.ANNO_PKGS <- c(
  K850 = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  K450 = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  K27  = "IlluminaHumanMethylation27kanno.ilmn12.hg19"
)

#' Load a named sub-table from an Illumina annotation package
#'
#' Each data object (e.g. "Locations", "Islands.UCSC", "Other") is stored as an
#' independent lazy dataset in the package and must be loaded via \code{data()}.
#' Accessing the top-level S4 wrapper and then reading \code{@data$Locations}
#' only yields a lazy descriptor, not the actual table.
#'
#' @param pkg   Character scalar: annotation package name.
#' @param table Character scalar: name of the dataset to load (e.g. \code{"Locations"}).
#' @return A \code{data.frame} with probe IDs as rownames.
#' @keywords internal
.anno_pkg_load_table <- function(pkg, table) {
  env <- new.env(parent = emptyenv())
  utils::data(list = table, package = pkg, envir = env)
  as.data.frame(get(table, envir = env), stringsAsFactors = FALSE)
}

#' Get probe IDs from an Illumina annotation package
#'
#' Extracts the complete vector of probe identifiers (e.g. \code{cg00000029})
#' from an annotation package.
#'
#' @param pkg Character scalar: annotation package name.
#' @return Character vector of probe IDs (rownames of the Locations table).
#' @keywords internal
.anno_pkg_probe_ids <- function(pkg) {
  locs <- .anno_pkg_load_table(pkg, "Locations")
  rownames(locs)
}

#' Build a data.frame from an Illumina annotation package
#'
#' Loads the Locations, Islands.UCSC, and Other sub-tables from the annotation
#' package and combines them into a single data.frame with one row per probe.
#' Does not require minfi.
#'
#' @param pkg Character scalar: annotation package name.
#' @return A \code{data.frame} with probe IDs as rownames and columns:
#'   \code{chr}, \code{pos}, \code{strand}, \code{Islands_Name},
#'   \code{Relation_to_Island}, \code{UCSC_RefGene_Name},
#'   \code{UCSC_RefGene_Group}, and further columns from the Other table.
#' @keywords internal
.anno_pkg_to_df <- function(pkg) {

  # Core genomic positions - always present
  locs <- .anno_pkg_load_table(pkg, "Locations")

  # CpG island context - present in all three arrays
  islands <- tryCatch(
    .anno_pkg_load_table(pkg, "Islands.UCSC"),
    error = function(e) {
      data.frame(
        Islands_Name       = NA_character_,
        Relation_to_Island = NA_character_,
        row.names          = rownames(locs),
        stringsAsFactors   = FALSE
      )
    }
  )

  # Gene body and other annotations
  other <- tryCatch(
    .anno_pkg_load_table(pkg, "Other"),
    error = function(e) data.frame(row.names = rownames(locs))
  )

  # Combine - all tables share the same rownames (probe IDs)
  cbind(locs, islands, other)
}

#' Detect Illumina array technology by probe-ID overlap
#'
#' Queries each installed annotation package and counts how many probe IDs
#' from \code{probe_ids} are present in each array's probe list.  Returns the
#' technology key (\code{"K27"}, \code{"K450"}, or \code{"K850"}) with the
#' highest overlap count.
#'
#' Returns \code{""} if none of the annotation packages are installed.
#'
#' @param probe_ids Character vector of probe identifiers from the signal matrix.
#' @return Named integer vector of overlap counts, or \code{""} if no packages
#'   are available.
#' @keywords internal
.anno_detect_tech_from_anno <- function(probe_ids) {

  counts <- integer(0)

  for (tech in names(.ANNO_PKGS)) {
    pkg <- .ANNO_PKGS[[tech]]
    if (requireNamespace(pkg, quietly = TRUE)) {
      n <- sum(probe_ids %in% .anno_pkg_probe_ids(pkg))
      counts[[tech]] <- n
    }
  }

  if (length(counts) == 0 || max(counts) == 0)
    return("")

  names(which.max(counts))
}
