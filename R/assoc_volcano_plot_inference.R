#' Volcano plot of association results for one inference_detail row.
#'
#' Reads the inference CSV that corresponds to the given `inference_detail`
#' row (already written by `association_analysis()` into
#' `<result_folder>/Inference/`), splits it by `(AREA, SUBAREA)`, and
#' produces one PNG volcano per combination under
#' `<result_folder>/Chart/VOLCANO/`. Naming mirrors the inference CSV
#' convention so the visual artifact is one-to-one traceable to the
#' analytic output:
#'
#'   `<MARKER>_<IV>_<transformation_y>_<family>_<covariates>_<areas_sql_condition>_<AREA>_<SUBAREA>.png`
#'
#' (passed through `core_name_cleaning()` which uppercases + replaces
#' comparison operators with `_GT_` / `_LT_` / `_EQ_` etc.)
#'
#' @param inference_detail A single row of `inference_details` (data.frame
#'   or list) - the same shape consumed by `association_analysis()`. Must
#'   carry `independent_variable`, `family_test`, `covariates`,
#'   `covariates_dummy`, `transformation_y`,
#'   `areas_sql_condition`, `samples_sql_condition`.
#' @param result_folder Project results folder (e.g.
#'   `~/.../results/GSE225845`). The function reads CSVs from
#'   `<result_folder>/Inference/` and writes PNGs under
#'   `<result_folder>/Chart/VOLCANO/`.
#' @param markers Character vector of marker names to plot for this scheda
#'   (e.g. `c("SIGNAL","DELTARP","DELTARQ")` for limma_2 PROBE). If NULL,
#'   inferred from the CSVs that exist in `Inference/` matching this
#'   scheda's metadata.
#' @param alpha Significance threshold drawn as a horizontal dashed line
#'   at `-log10(alpha)`. Points with `PVALUE_ADJ <= alpha` are coloured
#'   red (significant), the rest grey (non-significant). Default 0.05.
#' @param top_n_label How many of the top-significant points (by lowest
#'   `PVALUE_ADJ`) get an `AREA_OF_TEST` text label. Default 20. Set to 0
#'   to disable labels.
#' @param width,height,units Passed to `ggplot2::ggsave()`. Default
#'   9 × 9 inches.
#' @param dpi Plot resolution. Defaults to `ssEnv$plot_resolution_ppi`
#'   (typically 600).
#' @param overwrite If FALSE (default) skip PNGs that already exist.
#' @param pvalue_column Name of the p-value column to plot on the y axis. NULL
#'   lets the function pick the adjusted column matching the inference request
#'
#' @return Invisibly, a character vector of the PNG paths written
#'   (or that would have been written, when `overwrite = FALSE` and
#'   the file already exists).
#'
#' @export
#' @examples
#' # Stub: see vignette('imprinting-disorders', package = 'SEMseeker') for a
#' # runnable Beckwith-Wiedemann workflow on the GSE133774 subset.
#' invisible(NULL)
assoc_volcano_plot_inference <- function(inference_detail,
                                    result_folder,
                                    markers       = NULL,
                                    pvalue_column = NULL,
                                    alpha         = 0.05,
                                    top_n_label   = 20L,
                                    width         = 9,
                                    height        = 9,
                                    units         = "in",
                                    dpi           = NULL,
                                    overwrite     = FALSE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("assoc_volcano_plot_inference requires the 'ggplot2' package. Install it ",
         "with install.packages('ggplot2').")
  }
  use_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

  ssEnv <- tryCatch(core_get_session_info(), error = function(e) NULL)
  if (is.null(dpi)) {
    dpi <- if (!is.null(ssEnv$plot_resolution_ppi))
      as.numeric(ssEnv$plot_resolution_ppi) else 600
  }
  # AI-044 (2026-06-09): use SEMseeker's pastel palette (`color_palette`)
  # consistently with the rest of the package (plot_box_plot, fitted model
  # charts). `color_palette_darker` is a base-R name vector used for
  # accent strokes (e.g. trend lines); raw "red" / "blue" look harsh
  # against pastel-filled charts.
  palette_main <- if (!is.null(ssEnv$color_palette))
    ssEnv$color_palette else
    c("#b9e192", "#b3c7f7", "#f8b8d0", "#f194b8", "#ffefb6", "#cfebb6")
  palette_dark <- if (!is.null(ssEnv$color_palette_darker))
    ssEnv$color_palette_darker else c("blue", "red", "purple", "green")
  # Significant = vivid pink pastel (#f194b8, 4th); non-significant =
  # soft green pastel (#cfebb6, 6th). Threshold dashed line = soft blue
  # pastel (#b3c7f7, 2nd) so it reads as guide rather than focus.
  col_signif      <- palette_main[4]
  col_nonsignif   <- palette_main[6]
  col_threshold   <- palette_main[2]

  if (is.data.frame(inference_detail)) {
    inference_detail <- as.list(inference_detail[1L, , drop = FALSE])
  }

  # AI-044 (2026-06-09): pvalue_column resolution order:
  #   1. explicit `pvalue_column` argument (caller override)
  #   2. `inference_detail$pvalue_column` (forward-compat - not yet in
  #      assoc_validate_inference_schema but accepted if user supplies)
  #   3. default `PVALUE_ADJ_ALL_FDR` - the one column SEMseeker writes
  #      for every association run regardless of family/engine.
  if (is.null(pvalue_column) || !nzchar(pvalue_column)) {
    pvalue_column <- if (!is.null(inference_detail$pvalue_column) &&
                         nzchar(as.character(inference_detail$pvalue_column)))
      as.character(inference_detail$pvalue_column) else "PVALUE_ADJ_ALL_FDR"
  }

  iv               <- as.character(inference_detail$independent_variable)
  family_test      <- as.character(inference_detail$family_test)
  covariates       <- as.character(inference_detail$covariates)
  covariates_dummy <- as.character(inference_detail$covariates_dummy)
  transformation_y <- as.character(inference_detail$transformation_y)
  if (!nzchar(transformation_y)) transformation_y <- "none"
  areas_sql        <- as.character(inference_detail$areas_sql_condition)
  if (length(areas_sql) == 0L || !nzchar(areas_sql)) areas_sql <- ""

  inference_folder <- file.path(result_folder, "Inference")
  chart_folder     <- io_dir_check_and_create(file.path(result_folder, "Chart"),
                                            c("VOLCANO"))

  # If markers not supplied, scan Inference/ for CSVs that match the
  # current inference_detail shape (same IV and family) and extract
  # the leading MARKER token from the file basename.
  if (is.null(markers) || length(markers) == 0L) {
    candidate_csvs <- list.files(inference_folder, pattern = "\\.csv$",
                                  full.names = FALSE)
    iv_token       <- core_name_cleaning(iv)
    fam_token      <- core_name_cleaning(family_test)
    matches        <- candidate_csvs[
      grepl(iv_token, candidate_csvs, fixed = TRUE) &
      grepl(fam_token, candidate_csvs, fixed = TRUE)
    ]
    if (length(matches) == 0L) {
      core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                " assoc_volcano_plot_inference: no inference CSV matches IV=", iv,
                " family=", family_test, " in ", inference_folder)
      return(invisible(character(0)))
    }
    markers <- unique(vapply(matches, function(fn) {
      strsplit(fn, "_DEPTH_", fixed = TRUE)[[1]][1]
    }, character(1)))
  }

  written <- character(0)

  for (marker in markers) {
    # Use canonical SEMseeker inference filename builder so the lookup
    # mirrors what `association_analysis()` wrote (handles covariates_dummy
    # split-by-"+", dummy/pca tokens, samples_sql_condition subfolder).
    csv_path <- tryCatch(
      io_inference_file_name(inference_detail, marker, inference_folder,
                          file_extension = "csv"),
      error = function(e) {
        core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                  " assoc_volcano_plot_inference: io_inference_file_name failed for ",
                  marker, ": ", conditionMessage(e),
                  " - falling back to manual list.files() scan.")
        NULL
      }
    )
    # csv_basename derived from the actual path (works for both lookup
    # success above and the manual-scan path further down).
    csv_basename <- if (!is.null(csv_path))
      tools::file_path_sans_ext(basename(csv_path)) else
      paste(c(marker,
              core_name_cleaning(iv),
              core_name_cleaning(transformation_y),
              core_name_cleaning(family_test)), collapse = "_")

    if (is.null(csv_path) || !file.exists(csv_path)) {
      # Permissive fallback: scan Inference/ for a file starting with
      # `<MARKER>_<IV>_` and matching family. Useful when
      # the inference_detail row has changed slightly since the CSV was
      # written (e.g. additional dummy covariates added downstream).
      prefix <- core_name_cleaning(paste(c(marker, iv), collapse = "_"))
      fam_tok <- core_name_cleaning(family_test)
      candidates <- list.files(inference_folder, pattern = "\\.csv$",
                                full.names = TRUE)
      candidates <- candidates[startsWith(basename(candidates), prefix) &
                                grepl(fam_tok, basename(candidates), fixed = TRUE)]
      if (length(candidates) == 1L) {
        csv_path     <- candidates[1]
        csv_basename <- tools::file_path_sans_ext(basename(csv_path))
        core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                  " assoc_volcano_plot_inference: matched ", basename(csv_path),
                  " for marker '", marker, "' via prefix scan.")
      } else if (length(candidates) > 1L) {
        core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                  " assoc_volcano_plot_inference: ambiguous match for marker '",
                  marker, "' - ", length(candidates),
                  " CSVs share the prefix. Skipping.")
        next
      } else {
        core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                  " assoc_volcano_plot_inference: no CSV for marker '", marker,
                  "' in ", inference_folder, " - skipping.")
        next
      }
    }

    df <- tryCatch(
      utils::read.csv2(csv_path, header = TRUE, stringsAsFactors = FALSE,
                        check.names = FALSE),
      error = function(e) {
        core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                  " assoc_volcano_plot_inference: failed to read ", csv_path, ": ",
                  conditionMessage(e))
        NULL
      }
    )
    if (is.null(df) || nrow(df) == 0L) next

    needed_cols <- c("AREA", "SUBAREA", "AREA_OF_TEST", pvalue_column)
    missed <- setdiff(needed_cols, colnames(df))
    if (length(missed) > 0L) {
      core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                " assoc_volcano_plot_inference: ", csv_path,
                " missing columns: ", paste(missed, collapse = ", "),
                " (pvalue_column='", pvalue_column, "') - skipping.")
      next
    }

    # The X axis is the first non-intercept coefficient ESTIMATE. PVALUE
    # (top-level) was set by the engine to the matching coefficient's
    # raw p-value, and PVALUE_ADJ is BH-adjusted on PVALUE. So the
    # *_ESTIMATE column with the same coefficient name as the one driving
    # PVALUE is the right X.
    estimate_col <- .assoc_volcano_pick_estimate_col(df)
    if (is.null(estimate_col)) {
      core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                " assoc_volcano_plot_inference: ", csv_path,
                " - could not identify primary ESTIMATE column. Skipping.")
      next
    }

    df$._x      <- suppressWarnings(as.numeric(df[[estimate_col]]))
    df$._pval   <- suppressWarnings(as.numeric(df[[pvalue_column]]))
    df$._y      <- -log10(pmax(df$._pval, .Machine$double.xmin))
    df$._signif <- !is.na(df$._pval) & df$._pval <= alpha

    threshold_y <- -log10(alpha)

    for (area in unique(df$AREA)) {
      for (subarea in unique(df$SUBAREA[df$AREA == area])) {
        sub_df <- df[df$AREA == area & df$SUBAREA == subarea, , drop = FALSE]
        sub_df <- sub_df[is.finite(sub_df$._x) & is.finite(sub_df$._y), ,
                          drop = FALSE]
        if (nrow(sub_df) == 0L) next

        png_basename <- core_name_cleaning(paste(
          c(csv_basename, area, subarea), collapse = "_"
        ))
        png_path <- file.path(chart_folder, paste0(png_basename, ".png"))
        if (!overwrite && file.exists(png_path)) {
          written <- c(written, png_path)
          next
        }

        p <- ggplot2::ggplot(sub_df,
                              ggplot2::aes(x = .data[["._x"]],
                                            y = .data[["._y"]],
                                            colour = .data[["._signif"]])) +
          ggplot2::geom_point(alpha = 0.6, size = 1.4) +
          ggplot2::geom_hline(yintercept = threshold_y,
                               linetype = "dashed",
                               colour = col_threshold) +
          ggplot2::scale_colour_manual(
            values = c("FALSE" = col_nonsignif, "TRUE" = col_signif),
            labels = c("FALSE" = sprintf("%s > %.3g", pvalue_column, alpha),
                       "TRUE"  = sprintf("%s <= %.3g", pvalue_column, alpha)),
            name = NULL) +
          ggplot2::labs(
            title = util_pretty_label(sprintf(
              "%s - %s / %s - %s vs %s",
              marker, area, subarea, iv, family_test)),
            subtitle = util_pretty_label(sprintf(
              "areas_sql=%s | covariates=%s",
              if (nzchar(areas_sql)) areas_sql else "none",
              paste(c(covariates, covariates_dummy)[
                nzchar(c(covariates, covariates_dummy))], collapse = " + "))),
            x = util_pretty_label(sprintf("ESTIMATE (%s)", estimate_col)),
            y = util_pretty_label(sprintf("-log10(%s)", pvalue_column))
          ) +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position = "bottom")

        if (top_n_label > 0L) {
          n_pick <- min(top_n_label, sum(sub_df$._signif, na.rm = TRUE))
          if (n_pick > 0L) {
            top_idx <- order(sub_df$._pval, na.last = NA)[seq_len(n_pick)]
            top_df  <- sub_df[top_idx, , drop = FALSE]
            if (use_ggrepel) {
              p <- p + ggrepel::geom_text_repel(
                data = top_df,
                ggplot2::aes(label = .data[["AREA_OF_TEST"]]),
                size = 2.7, max.overlaps = top_n_label, show.legend = FALSE)
            } else {
              p <- p + ggplot2::geom_text(
                data = top_df,
                ggplot2::aes(label = .data[["AREA_OF_TEST"]]),
                size = 2.7, vjust = -0.5, show.legend = FALSE)
            }
          }
        }

        tryCatch(
          ggplot2::ggsave(filename = png_path, plot = p,
                          width = width, height = height,
                          units = units, dpi = dpi),
          error = function(e) {
            core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                      " assoc_volcano_plot_inference: ggsave failed for ", png_path,
                      ": ", conditionMessage(e))
          }
        )
        written <- c(written, png_path)
      }
    }
  }

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " assoc_volcano_plot_inference: wrote ", length(written),
            " PNG(s) under ", chart_folder)
  invisible(written)
}

# Pick the column to use as X for the volcano. Strategy:
#   1. Look for the coefficient that drove the top-level PVALUE column
#      (engines set PVALUE = first non-intercept coef p-value, the same
#      coefficient owns *_ESTIMATE which is the natural X).
#   2. If not directly inferable, take the first *_ESTIMATE column whose
#      *_PVALUE column equals PVALUE row-by-row (modulo rounding).
#   3. Fall back to the first non-INTERCEPT *_ESTIMATE column found.
.assoc_volcano_pick_estimate_col <- function(df) {
  est_cols <- grep("_ESTIMATE$", colnames(df), value = TRUE)
  if (length(est_cols) == 0L) return(NULL)
  est_cols <- est_cols[!grepl("^INTERCEPT|_INTERCEPT_", est_cols)]
  if (length(est_cols) == 0L) return(NULL)
  if (!"PVALUE" %in% colnames(df)) return(est_cols[1])

  for (ec in est_cols) {
    pc <- sub("_ESTIMATE$", "_PVALUE", ec)
    if (pc %in% colnames(df)) {
      a <- suppressWarnings(as.numeric(df[[pc]]))
      b <- suppressWarnings(as.numeric(df[["PVALUE"]]))
      ok <- !is.na(a) & !is.na(b)
      if (any(ok) && all(abs(a[ok] - b[ok]) < 1e-8, na.rm = TRUE)) {
        return(ec)
      }
    }
  }
  est_cols[1]
}
