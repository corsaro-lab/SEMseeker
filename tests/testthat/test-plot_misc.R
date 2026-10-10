## Coverage for pure helpers behind the plotting layer.
## The plot RENDERERS themselves read pipeline pivots / inference files, so they
## need full-pipeline fixtures and are not exercised here. The marker-value and
## threshold charts are the exception: their drawings take their data as
## arguments and are covered in test-plot_lollipop_build.R,
## test-sem_marker_value_per_probe_plot.R, test-sem_marker_value_per_sample_plot.R and
## test-sem_signal_threshold_per_sample_plot.R. These two internal helpers are pure and
## deterministic.
##
## Covered:
##   .assoc_volcano_pick_estimate_col()  choose the estimate column to plot
##   .plot_comparison_pvalue_label()     format a group-comparison p-value label

# ---------------------------------------------------------------------------
# .assoc_volcano_pick_estimate_col
# ---------------------------------------------------------------------------

test_that(".assoc_volcano_pick_estimate_col returns NULL when there is no usable estimate", {
  expect_null(SEMseeker:::.assoc_volcano_pick_estimate_col(data.frame(X = 1, Y = 2)))
  # intercept estimates are excluded
  expect_null(SEMseeker:::.assoc_volcano_pick_estimate_col(
    data.frame(INTERCEPT_ESTIMATE = 1)))
})

test_that(".assoc_volcano_pick_estimate_col falls back to the first estimate w/o PVALUE", {
  expect_equal(
    SEMseeker:::.assoc_volcano_pick_estimate_col(data.frame(AGE_ESTIMATE = 0.5)),
    "AGE_ESTIMATE")
})

test_that(".assoc_volcano_pick_estimate_col matches the estimate whose PVALUE equals PVALUE", {
  df <- data.frame(
    AGE_ESTIMATE = 0.5, AGE_PVALUE = 0.01,
    BMI_ESTIMATE = 0.2, BMI_PVALUE = 0.40,
    PVALUE       = 0.01)
  expect_equal(SEMseeker:::.assoc_volcano_pick_estimate_col(df), "AGE_ESTIMATE")
})

# ---------------------------------------------------------------------------
# .plot_comparison_pvalue_label
# ---------------------------------------------------------------------------

test_that(".plot_comparison_pvalue_label formats a p-value for supported tests", {
  set.seed(1)
  df <- data.frame(
    VALUE = c(stats::rnorm(10, 0), stats::rnorm(10, 4)),
    GRP   = rep(c("a", "b"), each = 10))

  expect_match(
    SEMseeker:::.plot_comparison_pvalue_label(df, "GRP", "VALUE", "t.test"),
    "^p = ")
  expect_match(
    SEMseeker:::.plot_comparison_pvalue_label(df, "GRP", "VALUE", "kruskal.test"),
    "^p = ")
})

test_that(".plot_comparison_pvalue_label returns NA for an unsupported family", {
  df <- data.frame(VALUE = c(1, 2, 3, 4), GRP = c("a", "a", "b", "b"))
  expect_true(is.na(SEMseeker:::.plot_comparison_pvalue_label(df, "GRP", "VALUE", "nope")))
})

# ---------------------------------------------------------------------------
# .assoc_volcano_marker_from_name
#
# The volcano chart finds which markers to draw by reading the inference file
# names, and it read the part before "_DEPTH_" to do it. The names carried that
# token while the granularity of an artefact was an integer; SCOPE, AREA and
# SUBAREA say it now, so the token is gone, strsplit() returned the whole name, and
# the marker became the file name. The lookup that follows then matched nothing, so
# the function drew nothing and reported that it had found no match - a failure
# that looks like an empty result.
# ---------------------------------------------------------------------------

test_that(".assoc_volcano_marker_from_name reads the marker off the leading token", {
  f <- SEMseeker:::.assoc_volcano_marker_from_name

  # The shape io_inference_file_name() writes: marker, independent variable,
  # transformation, family.
  expect_equal(f("MUTATIONS_AGE_GAUSSIAN.csv"), "MUTATIONS")
  expect_equal(f("DELTAS_TCDD_MOTHER_SCALE_WILCOXON.csv"), "DELTAS")

  # A full path is answered by its basename, not by the first directory.
  expect_equal(f(file.path("study", "Inference", "LESIONS_AGE_GAUSSIAN.csv")),
               "LESIONS")

  # And the shape that used to be parsed: the marker is still the leading token,
  # so a name left over from the old scheme reads correctly rather than specially.
  expect_equal(f("MUTATIONS_DEPTH_3_AGE_GAUSSIAN.csv"), "MUTATIONS")

  # No other token: the marker, without the extension stuck to it.
  expect_equal(f("MUTATIONS.csv"), "MUTATIONS")
})
