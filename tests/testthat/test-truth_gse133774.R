# =============================================================================
# The biological truth: it runs FIRST and stops everything when it fails.
#
# dev/run-suite.R runs this file before the rest of the suite, and the CI runs
# it as a step before R CMD check: if the method no longer finds the known
# epimutations, or finds them where they are not, nothing the other files say
# about the code matters. It takes seconds: thresholds on the Reference and a
# count of the probes below them, no pipeline run.
#
# The acceptance shares (90%, 80%, 10%) are a choice of this test, not of the
# paper; their rationale is still to be settled with care.
# =============================================================================

.truth_session <- function(env = parent.frame()) {
  tf <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE,
                            parallel_strategy = "sequential")
  ssEnv <- SEMseeker:::core_get_session_info()
  ssEnv$beta <- TRUE
  SEMseeker:::core_update_session_info(ssEnv)
  withr::defer({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, envir = env)
}

# -----------------------------------------------------------------------------
# The biological check: GSE133774 (Sparago et al., Clin Epigenetics 2019,
# PMID 31829238). L1 is the BWS proband (III-1), L3 his healthy brother with
# MLID (III-2), L2 the mother (II-2, carrier), L4 the father (II-1, healthy,
# not a carrier). CTRL01-06 are the Reference.
#
# What the paper established, and what must therefore hold:
#   - L1 and L3: loss of methylation at KCNQ1OT1:TSS-DMR (COBRA, pyrosequencing,
#     array); MEST in both; PLAGL1 in L1 and NOT in L3 (pyrosequencing);
#     H19/IGF2 in both, more in L3.
#   - L4, the father: no imprinting defect.
# An epimutation is a beta below the inferior threshold, counted here by hand.
# -----------------------------------------------------------------------------

test_that("GSE133774: the affected siblings are found and the father is not", {
  .truth_session()
  utils::data("test_master_features", package = "SEMseeker", envir = environment())
  utils::data("test_signal_gse133774", package = "SEMseeker", envir = environment())
  features <- as.data.frame(test_master_features, stringsAsFactors = FALSE)
  beta <- as.matrix(test_signal_gse133774)
  beta <- beta[intersect(rownames(beta), features$PROBE), ]

  reference <- beta[, paste0("CTRL0", 1:6)]
  th <- SEMseeker:::sem_signal_range_values(reference, "gse133774_truth",
                                            features[features$PROBE %in% rownames(beta), ])
  low <- stats::setNames(th$signal_inferior_thresholds, th$PROBE)

  region <- function(label) features$PROBE[which(features$DMR_LABEL == label)]
  lom <- function(sample, label) {
    p <- intersect(region(label), names(low))
    mean(beta[p, sample] < low[p])
  }

  # found: the share of the region's probes below the threshold
  expect_gte(lom("L1", "KCNQ1OT1:TSS-DMR"), 0.9)
  expect_gte(lom("L3", "KCNQ1OT1:TSS-DMR"), 0.9)
  expect_gte(lom("L1", "MEST:alt-TSS-DMR"), 0.8)
  expect_gte(lom("L3", "MEST:alt-TSS-DMR"), 0.8)
  expect_gte(lom("L1", "PLAGL1:alt-TSS-DMR"), 0.8)
  expect_gt(lom("L3", "H19_IGF2:IG-DMR"), lom("L1", "H19_IGF2:IG-DMR"))

  # not found
  expect_lte(lom("L3", "PLAGL1:alt-TSS-DMR"), 0.1)
  for (label in c("KCNQ1OT1:TSS-DMR", "H19_IGF2:IG-DMR", "PLAGL1:alt-TSS-DMR",
                  "MEST:alt-TSS-DMR", "MEG3_DLK1:IG-DMR", "GNAS A_B:TSS-DMR"))
    expect_lte(lom("L4", label), 0.1, label = paste("father at", label))
})
