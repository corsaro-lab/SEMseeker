# AI-185. The function is internal (@keywords internal @noRd), so it is NOT on
# the search path of the installed package: every call needs the ::: form.
# devtools::load_all() puts internals in scope and hides this, which is why
# these tests passed locally and every one of them errored in CI with
# "could not find function". The original commit had the same defect, so this
# file had never once run against an installed package.
# AI-185 (2026-06-26): unit tests for SEMseeker:::assoc_statistical_power().
#
# This exercises the STATISTICAL MACHINERY (per-probe partial f^2 + non-central
# F power + prospective design-power verdict), so it deliberately uses a small
# SYNTHETIC, fully-controlled signal matrix with a known association structure.
# (The real-GEO-fixture rule applies to biological claims, not to testing the
# math of a power formula.) The key correctness anchor is exact agreement with
# pwr::pwr.f2.test() on the same observed f^2.

make_signal <- function(n, beta_effect, seed = 1) {
  set.seed(seed)
  exposure <- scale(seq_len(n))[, 1]                 # continuous, mean 0 sd 1
  # two probe classes: associated (effect) and null (noise only)
  assoc <- t(sapply(1:20, function(i) beta_effect * exposure + rnorm(n, sd = 1)))
  null  <- t(sapply(1:20, function(i)               rnorm(n, sd = 1)))
  signal <- rbind(assoc, null)
  rownames(signal) <- c(paste0("assoc_", 1:20), paste0("null_", 1:20))
  colnames(signal) <- paste0("S", seq_len(n))
  ss <- data.frame(Sample_ID = colnames(signal), exposure = exposure,
                   stringsAsFactors = FALSE)
  list(signal = signal, ss = ss)
}

test_that("associated probes get higher power than null probes (continuous exposure)", {
  d <- make_signal(n = 60, beta_effect = 0.8)
  res <- SEMseeker:::assoc_statistical_power(d$signal, d$ss, exposure = "exposure",
                                  alpha = 0.05)

  pp <- res$per_probe
  expect_equal(nrow(pp), 40L)
  expect_true(all(pp$u == 1L))                       # continuous -> 1 df
  expect_true(all(pp$n == 60L))

  mean_assoc <- mean(pp$power[grepl("^assoc", pp$PROBE)])
  mean_null  <- mean(pp$power[grepl("^null",  pp$PROBE)])
  expect_gt(mean_assoc, mean_null)
  expect_gt(mean_assoc, 0.8)                          # strong effect, n=60
})

test_that("per-probe power matches pwr::pwr.f2.test on the observed f2", {
  d <- make_signal(n = 50, beta_effect = 0.5)
  res <- SEMseeker:::assoc_statistical_power(d$signal, d$ss, exposure = "exposure",
                                  alpha = 0.05)
  pp <- res$per_probe
  v  <- res$summary$v

  # independent recomputation via pwr for a handful of probes
  for (i in c(1, 5, 21, 40)) {
    # suppressWarnings: pwr.f2.test internally calls pf(..., lower=) which trips
    # R's partial-argument-match warning — a cosmetic pwr-package quirk, not ours.
    expected <- suppressWarnings(
      pwr::pwr.f2.test(u = 1, v = v, f2 = pp$f2_partial[i],
                       sig.level = 0.05, power = NULL)$power)
    expect_equal(pp$power[i], expected, tolerance = 1e-6,
                 info = paste("probe", pp$PROBE[i]))
  }
})

test_that("design-power verdict flips with sample size", {
  small <- make_signal(n = 12, beta_effect = 0.5)
  large <- make_signal(n = 400, beta_effect = 0.5)

  # target_f2 = 0.15 (medium): underpowered at n=12, correctly powered at n=400
  v_small <- SEMseeker:::assoc_statistical_power(small$signal, small$ss,
                                      exposure = "exposure", alpha = 0.05)$summary
  v_large <- SEMseeker:::assoc_statistical_power(large$signal, large$ss,
                                      exposure = "exposure", alpha = 0.05)$summary

  expect_equal(v_small$verdict, "underpowered")
  expect_equal(v_large$verdict, "correctly_powered")
  # the single indicator is the prospective design power, not the observed one
  expect_equal(v_large$design_power,
               suppressWarnings(pwr::pwr.f2.test(u = 1, v = v_large$v, f2 = 0.15,
                                sig.level = 0.05, power = NULL)$power),
               tolerance = 1e-6)
})

test_that("categorical exposure uses an (L-1)-df block", {
  d <- make_signal(n = 90, beta_effect = 0.4)
  # 3-level categorical exposure
  d$ss$grp <- factor(rep(c("A", "B", "C"), length.out = nrow(d$ss)))
  res <- SEMseeker:::assoc_statistical_power(d$signal, d$ss, exposure = "grp",
                                  alpha = 0.05)
  expect_equal(res$summary$exposure_type, "categorical")
  expect_equal(res$summary$u, 2L)                    # L - 1 = 2
  expect_true(all(res$per_probe$u == 2L))
})

test_that("covariate adjustment removes covariate-driven signal (partial f2)", {
  set.seed(7)
  n <- 80
  conf <- scale(rnorm(n))[, 1]                        # confounder
  exposure <- scale(rnorm(n))[, 1]
  # signal driven ONLY by the confounder, NOT by the exposure
  sig <- t(sapply(1:30, function(i) 1.5 * conf + rnorm(n, sd = 1)))
  rownames(sig) <- paste0("p_", 1:30)
  colnames(sig) <- paste0("S", seq_len(n))
  ss <- data.frame(Sample_ID = colnames(sig), exposure = exposure, conf = conf,
                   stringsAsFactors = FALSE)

  unadj <- SEMseeker:::assoc_statistical_power(sig, ss, exposure = "exposure", alpha = 0.05)
  adj   <- SEMseeker:::assoc_statistical_power(sig, ss, exposure = "exposure",
                                    covariates = "conf", alpha = 0.05)

  # exposure is null; adjusting for the confounder must not inflate its f2.
  # Median partial f2 of the (null) exposure stays near zero either way.
  expect_lt(stats::median(adj$per_probe$f2_partial), 0.1)
  expect_equal(adj$summary$k, 2L)                     # exposure + conf
  expect_equal(unadj$summary$k, 1L)
})
