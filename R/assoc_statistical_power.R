#' Statistical power of a per-probe regression on an exposure
#'
#' Computes the statistical power of associating a continuous signal matrix
#' (e.g. methylation beta values, probes x samples) with an **exposure**
#' variable taken from the sample sheet, under the per-probe linear model
#' \code{lm(beta ~ exposure + covariates)}. The exposure may be **continuous**
#' (a single 1-df term) or **categorical** (a block of \code{L - 1} dummy
#' terms, where \code{L} is the number of levels; the two-group case is the
#' special case \code{u = 1}).
#'
#' The effect size is Cohen's partial \eqn{f^2} of the exposure block, adjusted
#' for the covariates (Frisch-Waugh-Lovell residualisation). Power is obtained
#' from the non-central F distribution, identical to
#' \code{pwr::pwr.f2.test(u, v = n - k - 1, f2, sig.level = alpha)} but fully
#' vectorised across probes (no per-probe model loop).
#'
#' Two layers are returned:
#' \itemize{
#'   \item \strong{detail} (\code{per_probe}): the observed partial \eqn{f^2}
#'     and post-hoc power for every probe.
#'   \item \strong{summary indicator}: a single prospective \emph{design power}
#'     for a target effect size (\code{target_f2}) given the realised sample
#'     size and design, collapsed to a binary verdict
#'     \code{"correctly_powered"} / \code{"underpowered"} against
#'     \code{target_power}. This is a design (a-priori) quantity and does not
#'     suffer from the well-known circularity of observed (post-hoc) power.
#' }
#'
#' @param signal_data numeric matrix or data.frame, rows = probes
#'   (\code{rownames} = probe ids), columns = samples (\code{colnames} =
#'   \code{Sample_ID}). Continuous signal (e.g. beta values).
#' @param sample_sheet data.frame with at least a \code{Sample_ID} column plus
#'   the \code{exposure} column and any \code{covariates} columns.
#' @param exposure character. Name of the sample-sheet column to use as the
#'   independent variable (the exposure). NOT \code{Sample_Group} unless that
#'   is genuinely the exposure of interest.
#' @param covariates character vector of sample-sheet column names to adjust
#'   for. Default none.
#' @param alpha numeric significance level. If \code{NULL} (default) it is read
#'   from \code{ssEnv$alpha} when a session exists, else falls back to 0.05.
#' @param target_power numeric. Power threshold for the binary verdict
#'   (default 0.8).
#' @param target_f2 numeric. Target Cohen's \eqn{f^2} for the prospective
#'   design-power indicator (default 0.15, Cohen's "medium"; 0.02 small,
#'   0.35 large).
#' @param exposure_type one of \code{"auto"}, \code{"continuous"},
#'   \code{"categorical"}. \code{"auto"} treats numeric exposures as continuous
#'   and factor/character exposures as categorical.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{per_probe}}{data.frame with columns \code{PROBE}, \code{n},
#'       \code{u}, \code{f2_partial}, \code{power}, \code{powered}.}
#'     \item{\code{summary}}{list with the single \code{verdict}
#'       (\code{"correctly_powered"}/\code{"underpowered"}), the prospective
#'       \code{design_power}, the realised design (\code{n}, \code{u}, \code{k},
#'       \code{v}), \code{alpha}, \code{target_f2}, \code{target_power}, and
#'       descriptive observed-power statistics
#'       (\code{median_observed_power}, \code{frac_probes_powered}).}
#'   }
#'
#' @keywords internal
#' @noRd
assoc_statistical_power <- function(signal_data, sample_sheet, exposure,
                                     covariates = character(0), alpha = NULL,
                                     target_power = 0.8, target_f2 = 0.15,
                                     exposure_type = c("auto", "continuous",
                                                       "categorical")) {

  exposure_type <- match.arg(exposure_type)

  # --- resolve alpha (session default when available) -----------------------
  if (is.null(alpha)) {
    ssEnv <- tryCatch(core_get_session_info(), error = function(e) NULL)
    alpha <- if (!is.null(ssEnv) && !is.null(ssEnv$alpha))
      as.numeric(ssEnv$alpha) else 0.05
  }

  # --- input validation -----------------------------------------------------
  if (!"Sample_ID" %in% colnames(sample_sheet))
    stop("assoc_statistical_power: sample_sheet must have a 'Sample_ID' column.")
  if (length(exposure) != 1L || !exposure %in% colnames(sample_sheet))
    stop("assoc_statistical_power: 'exposure' must name a single sample_sheet column.")
  missing_cov <- setdiff(covariates, colnames(sample_sheet))
  if (length(missing_cov))
    stop("assoc_statistical_power: covariates not found in sample_sheet: ",
         paste(missing_cov, collapse = ", "))
  if (target_f2 <= 0)
    stop("assoc_statistical_power: 'target_f2' must be > 0.")

  signal_data <- as.matrix(signal_data)

  # --- align samples: sample_sheet rows <-> signal_data columns -------------
  samples <- intersect(as.character(sample_sheet$Sample_ID), colnames(signal_data))
  if (length(samples) < 3L)
    stop("assoc_statistical_power: fewer than 3 samples shared between ",
         "sample_sheet$Sample_ID and signal_data columns.")
  meta <- sample_sheet[match(samples, as.character(sample_sheet$Sample_ID)), ,
                       drop = FALSE]

  # drop samples with missing exposure or covariate values (listwise)
  needed <- c(exposure, covariates)
  complete <- stats::complete.cases(meta[, needed, drop = FALSE])
  samples  <- samples[complete]
  meta     <- meta[complete, , drop = FALSE]
  n <- length(samples)
  if (n < 4L)
    stop("assoc_statistical_power: fewer than 4 complete-case samples after ",
         "dropping missing exposure/covariate values.")

  # --- coerce exposure type -------------------------------------------------
  ex <- meta[[exposure]]
  is_categorical <- switch(exposure_type,
    auto        = is.factor(ex) || is.character(ex),
    continuous  = FALSE,
    categorical = TRUE)
  if (is_categorical) {
    ex <- droplevels(as.factor(ex))
    if (nlevels(ex) < 2L)
      stop("assoc_statistical_power: categorical exposure has < 2 levels after ",
           "complete-case filtering.")
  } else {
    ex <- as.numeric(ex)
  }
  meta[[exposure]] <- ex

  # --- build reduced (covariates only) and full (covariates + exposure) -----
  # Exposure terms are appended LAST so the added block is identifiable as the
  # trailing columns; u = ncol(full) - ncol(reduced).
  rhs_red  <- if (length(covariates)) paste(covariates, collapse = " + ") else "1"
  rhs_full <- paste(c(rhs_red, exposure), collapse = " + ")
  Xr <- stats::model.matrix(stats::as.formula(paste("~", rhs_red)),  data = meta)
  Xf <- stats::model.matrix(stats::as.formula(paste("~", rhs_full)), data = meta)

  u <- ncol(Xf) - ncol(Xr)              # df of the exposure block
  k <- ncol(Xf) - 1L                    # total predictors (excl. intercept)
  v <- n - ncol(Xf)                     # residual df = n - k - 1
  if (u < 1L)
    stop("assoc_statistical_power: exposure contributes no degrees of freedom ",
         "(constant within complete cases?).")
  if (v < 1L)
    stop("assoc_statistical_power: not enough residual degrees of freedom ",
         "(n - k - 1 = ", v, "); reduce covariates or add samples.")

  # --- vectorised per-probe partial f^2 via FWL residualisation -------------
  # Yt: samples x probes. qr.resid() residualises every probe column at once.
  Yt <- t(signal_data[, samples, drop = FALSE])
  qr_f <- qr(Xf); qr_r <- qr(Xr)
  rss_full <- colSums(qr.resid(qr_f, Yt)^2)
  rss_red  <- colSums(qr.resid(qr_r, Yt)^2)

  # Cohen's partial f^2 of the exposure block: (R2_full - R2_red)/(1 - R2_full)
  #   = (RSS_red - RSS_full) / RSS_full   (TSS cancels).
  f2 <- (rss_red - rss_full) / rss_full
  f2[!is.finite(f2)] <- NA_real_
  f2[f2 < 0] <- 0                       # clamp tiny negatives from rounding

  # --- observed (post-hoc) power per probe: non-central F (== pwr.f2.test) ---
  power <- .f2_power(f2, u = u, v = v, alpha = alpha)

  per_probe <- data.frame(
    PROBE      = rownames(signal_data),
    n          = n,
    u          = u,
    f2_partial = f2,
    power      = power,
    powered    = power >= target_power,
    stringsAsFactors = FALSE
  )

  # --- single prospective design-power indicator ----------------------------
  # Same non-central F machinery as the per-probe path, evaluated at the target
  # effect size; identical to pwr::pwr.f2.test(u, v, f2 = target_f2).
  design_power <- .f2_power(target_f2, u = u, v = v, alpha = alpha)
  verdict <- if (design_power >= target_power) "correctly_powered" else "underpowered"

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " Statistical power on exposure '", exposure, "' (",
            if (is_categorical) "categorical" else "continuous",
            ", u=", u, ", n=", n, ", alpha=", alpha,
            "): design power for f2=", target_f2, " is ",
            round(design_power, 3), " -> ", verdict)

  summary <- list(
    verdict               = verdict,
    design_power          = design_power,
    exposure              = exposure,
    exposure_type         = if (is_categorical) "categorical" else "continuous",
    n                     = n,
    u                     = u,
    k                     = k,
    v                     = v,
    alpha                 = alpha,
    target_f2             = target_f2,
    target_power          = target_power,
    n_probes              = nrow(per_probe),
    n_probes_na           = sum(is.na(power)),
    median_observed_power = stats::median(power, na.rm = TRUE),
    frac_probes_powered   = mean(power >= target_power, na.rm = TRUE)
  )

  list(per_probe = per_probe, summary = summary)
}

#' Power of an F-test from Cohen's f^2 (vectorised, non-central F)
#'
#' Closed-form power for the general linear-model F-test, identical to
#' \code{pwr::pwr.f2.test(u, v, f2, sig.level)$power} but vectorised over
#' \code{f2} and free of partial-argument-matching warnings. Used both for the
#' per-probe observed power and the scalar prospective design power so the two
#' layers are guaranteed mutually consistent.
#'
#' @param f2 numeric (vector) Cohen's \eqn{f^2} effect size(s).
#' @param u numerator degrees of freedom (size of the tested block).
#' @param v denominator degrees of freedom (\eqn{n - k - 1}).
#' @param alpha significance level.
#' @return numeric vector of power values, same length as \code{f2}.
#' @keywords internal
#' @noRd
.f2_power <- function(f2, u, v, alpha) {
  crit   <- stats::qf(alpha, df1 = u, df2 = v, lower.tail = FALSE)
  lambda <- f2 * (u + v + 1)
  stats::pf(crit, df1 = u, df2 = v, ncp = lambda, lower.tail = FALSE)
}
