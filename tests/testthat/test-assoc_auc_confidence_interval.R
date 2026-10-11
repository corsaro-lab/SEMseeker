# The AUC, its interval, and the power: one number that was reported under three
# names, two of which were wrong about it, and a power computed from the wrong
# effect size in two branches.

.auc_ci <- function(...) SEMseeker:::assoc_auc_confidence_interval(...)
.pow    <- function(...) SEMseeker:::assoc_power_from_d(...)

# ---- the interval -----------------------------------------------------------

test_that("assoc_auc_confidence_interval: the interval brackets the estimate and stays in [0,1]", {
  for (a in c(0.05, 0.3, 0.5, 0.76, 0.95)) {
    ci <- .auc_ci(a, 30, 40)
    expect_lt(ci$lower, a); expect_gt(ci$upper, a)
    expect_gte(ci$lower, 0); expect_lte(ci$upper, 1)
  }
})

test_that("assoc_auc_confidence_interval: it cannot leave the bounds even near them", {
  # The reason the interval is built on the logit rather than symmetrically on
  # the AUC: a Wald interval on a bounded quantity walks out of the bounds.
  ci <- .auc_ci(0.99, 6, 6)
  expect_lte(ci$upper, 1)
  ci <- .auc_ci(0.01, 6, 6)
  expect_gte(ci$lower, 0)
})

test_that("assoc_auc_confidence_interval: complete separation answers the point, not an interval", {
  for (a in c(0, 1)) {
    ci <- .auc_ci(a, 20, 20)
    expect_equal(ci$lower, a); expect_equal(ci$upper, a); expect_equal(ci$se, 0)
  }
})

test_that("assoc_auc_confidence_interval: a group with nobody in it answers NA", {
  for (args in list(list(0.7, 0, 10), list(0.7, 10, 0), list(NA, 10, 10)))
    expect_true(is.na(do.call(.auc_ci, args)$lower))
})

test_that("assoc_auc_confidence_interval: it narrows as the groups grow", {
  w <- vapply(c(10, 30, 100, 400), function(n) { ci <- .auc_ci(0.7, n, n); ci$upper - ci$lower }, 0)
  expect_true(all(diff(w) < 0))
})

test_that("assoc_auc_confidence_interval: it covers its nominal level", {
  # The property the formula is for, asserted rather than assumed. The logit form
  # was chosen over the symmetric one because the symmetric one covers 90% at ten
  # per group; this pins that the chosen one does not regress to that.
  skip_on_cran()
  set.seed(101)
  n1 <- 30; n2 <- 40; delta <- 1
  true_auc <- stats::pnorm(delta / sqrt(2))
  hit <- replicate(600, {
    g1 <- stats::rnorm(n1, delta); g2 <- stats::rnorm(n2, 0)
    a  <- mean(outer(g1, g2, ">")) + 0.5 * mean(outer(g1, g2, "=="))
    ci <- .auc_ci(a, n1, n2)
    !is.na(ci$lower) && true_auc >= ci$lower && true_auc <= ci$upper
  })
  expect_gt(mean(hit), 0.90)
})

# ---- the power --------------------------------------------------------------

test_that("assoc_power_from_d: no effect answers the significance level", {
  # The defect this replaces: the wilcoxon branch passed A, so the exact null
  # read as d = 0.5 and the reported power was 0.53.
  expect_equal(.pow(0, 30, 40, 0.05), 0.05, tolerance = 1e-6)
  expect_gt(.pow(0.5, 30, 40, 0.05), 0.5)   # what A = 0.5 used to produce
})

test_that("assoc_power_from_d: an infinite effect answers one, and an undefined one NA", {
  expect_equal(.pow(Inf, 30, 40, 0.05), 1, tolerance = 1e-6)
  expect_equal(.pow(-Inf, 30, 40, 0.05), 1, tolerance = 1e-6)
  expect_true(is.na(.pow(NA, 30, 40, 0.05)))
  # a group of one has no standardised difference
  expect_true(is.na(.pow(1, 1, 40, 0.05)))
})

test_that("assoc_power_from_d: the sign of the effect does not change the power", {
  expect_equal(.pow(1.2, 20, 20, 0.05), .pow(-1.2, 20, 20, 0.05))
})

# ---- A, 2A-1 and the AUC are three different claims about one number --------

test_that("the rank-biserial correlation is 2A-1 and not A", {
  # Measured: W/(n1*n2) equals Vargha-Delaney A exactly, so the column named for
  # the rank-biserial correlation held the AUC. The two differ in value and, here,
  # in sign.
  set.seed(42)
  g1 <- stats::rnorm(30, 1); g2 <- stats::rnorm(40, 2)
  w  <- stats::wilcox.test(g1, g2)
  a  <- as.numeric(effsize::VD.A(g1, g2)$estimate)

  expect_equal(as.numeric(w$statistic) / (30 * 40), a, tolerance = 1e-9)
  expect_false(isTRUE(all.equal(a, 2 * a - 1)))
  expect_lt(2 * a - 1, 0)
  expect_gt(a, 0)
})

# ---- through assoc_test_model, which is what a consumer reads ----------------

.tg <- function(seed = 3L, n = 20L, shift = 1.5) {
  set.seed(seed)
  data.frame(
    BURDEN = c(stats::rnorm(n, 1), stats::rnorm(n, 1 + shift)),
    GROUP  = factor(c(rep("ctrl", n), rep("case", n)), levels = c("ctrl", "case")))
}
.k <- function() list(AREA = "GENE", SUBAREA = "TSS200",
                      MARKER = "MUTATIONS", FIGURE = "K850")

test_that("assoc_test_model wilcoxon: the AUC is reported with its interval", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df  <- .tg()
  res <- SEMseeker:::assoc_test_model("wilcoxon", df, stats::as.formula("BURDEN ~ GROUP"),
                                      "BURDEN", "GROUP", "", FALSE, "", .k())

  for (col in c("C_STATISTIC_AUC", "C_STATISTIC_AUC_CI_LOWER", "C_STATISTIC_AUC_CI_UPPER"))
    expect_true(col %in% colnames(res), info = col)

  expect_lt(res$C_STATISTIC_AUC_CI_LOWER, res$C_STATISTIC_AUC)
  expect_gt(res$C_STATISTIC_AUC_CI_UPPER, res$C_STATISTIC_AUC)
  expect_gte(res$C_STATISTIC_AUC_CI_LOWER, 0)
  expect_lte(res$C_STATISTIC_AUC_CI_UPPER, 1)
})

test_that("assoc_test_model wilcoxon: RANK_BISERIAL_CORRELATION is 2A-1, not A", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df  <- .tg()
  res <- SEMseeker:::assoc_test_model("wilcoxon", df, stats::as.formula("BURDEN ~ GROUP"),
                                      "BURDEN", "GROUP", "", FALSE, "", .k())

  expect_equal(as.numeric(res$RANK_BISERIAL_CORRELATION),
               2 * as.numeric(res$C_STATISTIC_AUC) - 1, tolerance = 1e-9)
  # and it is no longer the same number as the AUC, which is what it used to hold
  expect_false(isTRUE(all.equal(as.numeric(res$RANK_BISERIAL_CORRELATION),
                                as.numeric(res$C_STATISTIC_AUC))))
  # the controls are the first level and the cases are shifted up, so A is below
  # a half from the controls' point of view and the rank-biserial is negative:
  # the sign is decided by the order of the levels, which the request declares
  expect_lt(res$RANK_BISERIAL_CORRELATION, 0)
})

test_that("assoc_test_model wilcoxon: two identical groups report the significance level, not 0.53", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # A = 0.5 exactly: the two groups are the same values, so there is no effect
  # and no power beyond the significance level. This used to report 0.53, the
  # power of a medium effect, because A was handed to pwr as Cohen's d.
  df <- data.frame(BURDEN = c(1:10, 1:10),
                   GROUP  = factor(rep(c("a", "b"), each = 10), levels = c("a", "b")))
  res <- SEMseeker:::assoc_test_model("wilcoxon", df, stats::as.formula("BURDEN ~ GROUP"),
                                      "BURDEN", "GROUP", "", FALSE, "", .k())

  expect_equal(as.numeric(res$C_STATISTIC_AUC), 0.5, tolerance = 1e-9)
  expect_equal(as.numeric(res$RANK_BISERIAL_CORRELATION), 0, tolerance = 1e-9)
  expect_lt(as.numeric(res$power), 0.1)
})

test_that("assoc_test_model t.test: the power no longer moves with the scale of the data", {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # The same comparison read on two scales. Cohen's d is the same either way, so
  # the power has to be. It was 0.053 on one scale and 0.171 on the other.
  set.seed(53); n <- 30
  b <- c(stats::rnorm(n, 0.50, 0.03), stats::rnorm(n, 0.55, 0.03))
  g <- factor(rep(c("a", "b"), each = n), levels = c("a", "b"))

  on_scale <- function(y) {
    res <- SEMseeker:::assoc_test_model("t.test", data.frame(BURDEN = y, GROUP = g),
                                        stats::as.formula("BURDEN ~ GROUP"),
                                        "BURDEN", "GROUP", "", FALSE, "", .k())
    as.numeric(res$power)
  }
  p_beta <- on_scale(b)
  p_mval <- on_scale(log2(b / (1 - b)))

  expect_equal(p_beta, p_mval, tolerance = 0.02)
  # and a difference this clear has power near one, where the raw difference of
  # the means on the beta scale used to report the significance level
  expect_gt(p_beta, 0.9)
})
