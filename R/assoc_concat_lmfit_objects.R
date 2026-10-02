# AI-061+ (2026-06-09): concatenate a list of `limma::MArrayLM` objects
# produced by chunked-per-chr lmFit calls into a single MArrayLM ready
# for a global `limma::eBayes()` call.
#
# Cross-gene eBayes shrinkage requires ALL sigma / coefficients /
# stdev.unscaled values together. If each chunk runs its own eBayes,
# every chunk uses a different prior variance and the moderation is
# inconsistent across the dataset. The right architectural answer is
# therefore: chunked lmFit + GLOBAL eBayes on the concatenated fit.
#
# This helper makes the concatenation safe by:
#   1. Concatenating per-row fields (coefficients, sigma, Amean,
#      stdev.unscaled, df.residual, weights if present) along the row
#      axis - same column / list semantics as a monolithic lmFit.
#   2. Carrying forward design-side fields (cov.coefficients, design,
#      pivot, qr, rank, method, call) from the first chunk, asserting
#      they are equal across chunks. Same design matrix produces the
#      same QR decomposition every time.
#   3. Setting rownames on the combined object so downstream consumers
#      (eBayes, topTable) see a coherent probe index.

#' Concatenate per-chr limma fits into a single MArrayLM
#'
#' @param fit_list A non-empty list of `limma::MArrayLM` objects, all
#'   produced from the SAME design matrix on disjoint row subsets.
#'
#' @return A `limma::MArrayLM` whose row dimension is the sum of input
#'   row dimensions and whose design-side fields are inherited from the
#'   first chunk. The returned object is safe to feed to
#'   `limma::eBayes()` for global variance moderation.
#'
#' @keywords internal
#' @noRd
assoc_concat_lmfit_objects <- function(fit_list) {
  if (!is.list(fit_list) || length(fit_list) == 0L) {
    stop("assoc_concat_lmfit_objects: fit_list must be a non-empty list.",
         call. = FALSE)
  }
  # Drop NULL entries (chunks that produced no rows after filtering).
  fit_list <- fit_list[!vapply(fit_list, is.null, logical(1))]
  if (length(fit_list) == 0L) return(NULL)
  if (length(fit_list) == 1L) return(fit_list[[1]])

  fit_first <- fit_list[[1]]
  if (!inherits(fit_first, "MArrayLM")) {
    stop("assoc_concat_lmfit_objects: every element must inherit from MArrayLM.",
         call. = FALSE)
  }

  # ---- design-side fields: assert equal, take from first chunk -----
  # These are functions of the design matrix only - must match across
  # chunks. Comparing $design is sufficient because the rest derive
  # from it.
  for (i in seq.int(2L, length(fit_list))) {
    f <- fit_list[[i]]
    if (!isTRUE(all.equal(f$design, fit_first$design,
                            check.attributes = FALSE,
                            tolerance        = 1e-10))) {
      stop("assoc_concat_lmfit_objects: chunk ", i,
           " was fitted with a different design matrix - refusing to ",
           "concatenate (would corrupt cov.coefficients and stdev.unscaled).",
           call. = FALSE)
    }
  }

  # ---- per-gene fields: rbind / c() along the row axis --------------
  out <- fit_first

  # Helper: pick `f[[field]]` across chunks and rbind / c() them.
  .rbind_field <- function(field) {
    pieces <- lapply(fit_list, function(f) f[[field]])
    if (is.null(pieces[[1]])) return(NULL)
    if (is.matrix(pieces[[1]])) do.call(rbind, pieces) else unlist(pieces, use.names = TRUE)
  }

  per_gene_matrix <- c("coefficients", "stdev.unscaled", "t", "p.value",
                       "lods", "weights")
  per_gene_vector <- c("sigma", "Amean", "df.residual", "s2.post", "F",
                       "F.p.value")

  for (f in per_gene_matrix) {
    if (!is.null(fit_first[[f]])) out[[f]] <- .rbind_field(f)
  }
  for (f in per_gene_vector) {
    if (!is.null(fit_first[[f]])) out[[f]] <- .rbind_field(f)
  }

  # Names: rownames of the row-axis fields.
  all_names <- unlist(lapply(fit_list, function(f) rownames(f$coefficients)),
                      use.names = FALSE)
  if (!is.null(all_names)) {
    if (is.matrix(out$coefficients)) rownames(out$coefficients) <- all_names
    if (is.matrix(out$stdev.unscaled)) rownames(out$stdev.unscaled) <- all_names
    if (!is.null(out$sigma))         names(out$sigma) <- all_names
    if (!is.null(out$Amean))         names(out$Amean) <- all_names
  }

  out
}
