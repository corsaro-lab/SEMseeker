#' Smart split of a multi-gene AREA name with prefix recovery
#'
#' Bioconductor and external annotations sometimes
#' encode multiple genes in a compact slash-separated form where ALL
#' suffix tokens inherit a prefix from the first token. The naive
#' `strsplit(s, "/")` loses the prefix and produces unintelligible
#' single-letter or digit-only tokens:
#'
#'   `"HLA-A/B/C"`  →  `c("HLA-A", "B", "C")`        (BAD)
#'   `"KRT8/18"`    →  `c("KRT8", "18")`             (BAD)
#'
#' This function reconstructs the missing prefix using one of two
#' heuristics, applied in order:
#'
#'   - Strategy 1: if the first token contains `-`, take everything up
#'     to and including the LAST `-` as the prefix:
#'        `"HLA-A"` → prefix `"HLA-"`
#'   - Strategy 2: otherwise, take the leading run of alphabetic
#'     characters as the prefix:
#'        `"KRT8"` → prefix `"KRT"`
#'
#' Each non-first token is then prefixed only if it does not already
#' start with the recovered prefix (avoids double-prepending on inputs
#' like `"HLA-A/HLA-B"` or `"HBA1/HBA2"`).
#'
#' @param s A single character string (one AREA name).
#'
#' @return A character vector with the expansion. If `s` does not
#'   contain `"/"`, returns `s` unchanged in a length-1 vector.
#'
#' @examples
#' # Internal helper - not exported. Reach it via ::: so R CMD check --as-cran
#' # can run the examples block under CheckExEnv (which only sees exports).
#' SEMseeker:::.anno_smart_split_area_name("HLA-A/B/C")    # c("HLA-A","HLA-B","HLA-C")
#' SEMseeker:::.anno_smart_split_area_name("HBA1/HBA2")    # c("HBA1","HBA2")
#' SEMseeker:::.anno_smart_split_area_name("KRT8/18")      # c("KRT8","KRT18")
#' SEMseeker:::.anno_smart_split_area_name("TP53")         # "TP53"
#' SEMseeker:::.anno_smart_split_area_name("GENE-A-B/C/D") # c("GENE-A-B","GENE-A-C","GENE-A-D")
#'
#' @keywords internal
.anno_smart_split_area_name <- function(s) {
  if (length(s) != 1L || is.na(s) || !nzchar(s)) return(s)
  if (!grepl("/", s, fixed = TRUE)) return(s)

  parts <- strsplit(s, "/", fixed = TRUE)[[1]]
  first <- parts[1]

  prefix <- if (grepl("-", first, fixed = TRUE)) {
    # Strategy 1: everything up to and including the last "-"
    sub("-[^-]+$", "-", first)
  } else {
    # Strategy 2: leading run of letters
    sub("[^A-Za-z].*$", "", first)
  }

  if (!nzchar(prefix)) return(parts)

  out <- first
  for (suf in parts[-1]) {
    if (!startsWith(suf, prefix)) suf <- paste0(prefix, suf)
    out <- c(out, suf)
  }
  out
}
