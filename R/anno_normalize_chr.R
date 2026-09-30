#' Normalize chromosome names between internal and external formats
#'
#' SEMseeker uses bare chromosome names internally (\code{"1"}, \code{"X"})
#' for fast joins and compact storage.
#' External formats (BED, bedgraph) require the UCSC \code{"chr"} prefix.
#'
#' This function is idempotent: stripping an already-bare name or prefixing
#' an already-prefixed name is a no-op.
#'
#' @param x Character vector of chromosome names.
#' @param direction One of \code{"internal"} (strip \code{"chr"} prefix)
#'   or \code{"output"} (add \code{"chr"} prefix).
#'
#' @return Character vector of normalized chromosome names.
#'
#' @examples
#' anno_normalize_chr(c("chr1", "chrX", "chr22"), "internal")
#' # => c("1", "X", "22")
#'
#' anno_normalize_chr(c("1", "X", "22"), "output")
#' # => c("chr1", "chrX", "chr22")
#'
#' # Idempotent:
#' anno_normalize_chr(c("1", "X"), "internal")
#' # => c("1", "X")  - already bare, unchanged
#'
#' anno_normalize_chr(c("chr1", "chrX"), "output")
#' # => c("chr1", "chrX")  - already prefixed, unchanged
#'
#' @keywords internal
#' @noRd
anno_normalize_chr <- function(x, direction = c("internal", "output")) {
  direction <- match.arg(direction)
  x <- as.character(x)
  if (length(x) == 0L) return(character(0))
  if (direction == "internal") {
    sub("^chr", "", x, ignore.case = TRUE)
  } else {
    needs_prefix <- !grepl("^chr", x, ignore.case = TRUE)
    x[needs_prefix] <- paste0("chr", x[needs_prefix])
    x
  }
}
