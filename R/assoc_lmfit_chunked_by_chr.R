# 2026-06-09: chunked-per-chromosome lmFit path for the
# batch-lazy dispatcher. Activated when the bulk_model memory gate
# decides the monolithic y_mat would exceed the RAM budget.
#
# Chromosome is the natural unit for chunking because:
#   1. it bounds the chunk size to a known fraction (~5%-20% of the
#      genome for human autosomes),
#   2. it composes naturally with the canonical CHR/START sort gate
#      written by io_signal_save() - each chr-block is already contiguous
#      and pre-sorted on disk, so the lazy filter is a single linear
#      scan,
#   3. it's tech-agnostic: Illumina probes have a chromosome
#      annotation via the Bioconductor manifest; long-read positions
#      carry chromosome in the coord-encoded AREA string ("chr1_<pos>").
#
# Filtering strategy (tech-aware, no pivot schema change required):
#   Illumina    : `pl$col("AREA")$is_in(probes_in_chr)`  (lookup via
#                 probe_features data.frame, ~12 MB R)
#   WGBS/LONGREAD: `pl$col("AREA")$str$starts_with("chr<i>_")` - chromosome
#                 parsed lazily from AREA string prefix, no manifest
#                 lookup needed.
#
# The cross-gene eBayes shrinkage is applied by the dispatcher AFTER
# `assoc_concat_lmfit_objects()` glues the per-chr fits together. Doing eBayes
# globally on the concatenated fit is numerically identical to running
# it on a monolithic lmFit because:
#   - cov.coefficients = (X' X)^{-1} depends only on design (identical
#     across chunks);
#   - eBayes consumes sigma / coefficients / stdev.unscaled / df.residual
#     all of which are concatenated verbatim by `assoc_concat_lmfit_objects()`.

#' Fit limma::lmFit one chromosome at a time on a lazy pivot
#'
#' @param pivot_lazy `polars_lazy_frame` already filtered for
#'   `area_to_remove` and transformed by `io_data_preparation_lazy()`.
#'   First column MUST be `AREA` (probe / position identifier); the
#'   remaining columns are samples.
#' @param sample_cols_kept Character vector - sample column names that
#'   survived the IV / covariates complete-cases filter. Used to subset
#'   the pivot lazily before materialisation, same as the monolithic
#'   path.
#' @param design Numeric matrix - passed to `limma::lmFit()` for every
#'   chunk. Must be the same design across chunks (the helper
#'   `assoc_concat_lmfit_objects()` asserts this).
#' @param engine Character - `"limma"` or `"voom"`. NOTE: voom chunked
#'   requires a global mean-variance trend, not implemented yet.
#'   Routing voom to this function emits a WARNING and falls back to
#'   monolithic (caller must ensure they don't hit this path).
#' @param key List with `MARKER`/`FIGURE`/`AREA`/`SUBAREA` - for log lines.
#' @param family_test Character - log identifier.
#' @param probe_features `data.frame` with columns `PROBE` and `CHR` -
#'   required for Illumina pivots (manifests the probe → chr lookup).
#'   Pass `NULL` for WGBS/LONGREAD where the AREA string is
#'   coord-encoded.
#' @param tech_is_longread Logical - TRUE for WGBS / LONGREAD, FALSE
#'   for Illumina. Decides the filter strategy (`is_in` vs
#'   `starts_with`).
#'
#' @return The concatenated pre-eBayes `MArrayLM` fit object across
#'   all chromosomes, or NULL if every chunk failed.
#'
#' @keywords internal
#' @noRd
assoc_lmfit_chunked_by_chr <- function(pivot_lazy,
                                  sample_cols_kept,
                                  design,
                                  engine,
                                  key,
                                  family_test,
                                  probe_features = NULL,
                                  tech_is_longread = FALSE) {

  if (engine == "voom") {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_lmfit_chunked_by_chr: voom chunking requires a global ",
              "mean-variance trend pass (not implemented yet). ",
              "Falling back to monolithic - caller must ensure the y_mat ",
              "fits within budget.")
    return(assoc_lmfit_monolithic_lazy(pivot_lazy, sample_cols_kept, design,
                                  engine, key, family_test))
  }

  # ---- enumerate chromosomes -----------------------------------------
  if (tech_is_longread) {
    # AREA is coord-encoded like "chr1_12345_12346". Extract CHR via
    # lazy string split (Rust-native, no R-side allocation).
    chrs <- as.character(as.data.frame(
      pivot_lazy$select(
        polars::pl$col("AREA")$str$split("_")$list$get(0L)$alias("CHR")
      )$unique()$collect()
    )$CHR)
    # Sort canonical: numeric chrs first, then X/Y/M.
    chr_key <- function(ch) {
      cl <- sub("^chr", "", ch)
      n  <- suppressWarnings(as.integer(cl))
      ifelse(!is.na(n), n,
             ifelse(cl == "X", 23L,
                    ifelse(cl == "Y", 24L,
                           ifelse(cl %in% c("M", "MT"), 25L, 99L))))
    }
    chrs <- chrs[order(chr_key(chrs))]
  } else {
    if (is.null(probe_features) ||
        !all(c("PROBE", "CHR") %in% colnames(probe_features))) {
      stop("assoc_lmfit_chunked_by_chr: Illumina path needs `probe_features` ",
           "with columns PROBE and CHR.",
           call. = FALSE)
    }
    probes_by_chr <- split(as.character(probe_features$PROBE),
                            as.character(probe_features$CHR))
    chrs <- names(probes_by_chr)
  }

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " assoc_lmfit_chunked_by_chr [", key$MARKER, "/", key$FIGURE, "/",
            key$AREA, "/", key$SUBAREA, "] engine=", engine,
            " chunks=", length(chrs), " (tech_is_longread=",
            tech_is_longread, ").")

  # ---- per-chr loop --------------------------------------------------
  fit_list <- vector("list", length(chrs))
  for (i in seq_along(chrs)) {
    ch <- chrs[i]
    t0 <- Sys.time()
    pivot_chunk <- if (tech_is_longread) {
      pivot_lazy$filter(
        polars::pl$col("AREA")$str$starts_with(paste0(ch, "_"))
      )
    } else {
      probes_chr <- probes_by_chr[[ch]]
      if (length(probes_chr) == 0L) next
      pivot_lazy$filter(
        polars::pl$col("AREA")$is_in(
          polars::pl$lit(probes_chr)$implode()
        )
      )
    }
    fit_list[[i]] <- assoc_lmfit_monolithic_lazy(pivot_chunk, sample_cols_kept,
                                            design, engine, key, family_test)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    n_rows  <- if (!is.null(fit_list[[i]]))
      nrow(fit_list[[i]]$coefficients) else 0L
    core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_lmfit_chunked_by_chr chunk ", ch, " (",
              i, "/", length(chrs), ") rows=", n_rows,
              " elapsed=", round(elapsed, 1), "s")
  }

  # ---- concat ---------------------------------------------------------
  fit <- assoc_concat_lmfit_objects(fit_list)
  if (is.null(fit)) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " assoc_lmfit_chunked_by_chr: every chunk produced an empty fit ",
              "for [", key$MARKER, "/", key$FIGURE, "/", key$AREA, "/",
              key$SUBAREA, "].")
    return(NULL)
  }
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " assoc_lmfit_chunked_by_chr [", key$MARKER, "/", key$FIGURE, "/",
            key$AREA, "/", key$SUBAREA, "] concat done, total rows=",
            nrow(fit$coefficients))
  fit
}
