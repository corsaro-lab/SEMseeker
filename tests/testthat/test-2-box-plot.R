## test-2-box-plot.R
## Tests for plot_box_plot(): assert PNG files are created for the box-plot
## and violin-plot variants. Requires ggpubr (in Suggests).
##
## plot_box_plot() is a side-effect function: it writes two PNG files per call
## (one boxplot, one violin) under the session chart folder. The tests below
## verify that both files exist after the call and that non-dichotomous
## families are silently skipped.

test_that("plot_box_plot: wilcoxon call creates boxplot and violin PNG files", {
  skip_if_not_installed("ggpubr")

  tf  <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  set.seed(1)
  df <- data.frame(
    BURDEN = c(stats::rnorm(15, 1), stats::rnorm(15, 3)),
    GROUP  = factor(c(rep("ctrl", 15), rep("case", 15)))
  )
  key <- list(AREA = "GENE", SUBAREA = "TSS200", MARKER = "MUTATIONS", FIGURE = "K850")

  SEMseeker:::plot_box_plot(
    dataFrameToPlot       = df,
    independent_variable  = "GROUP",
    dependent_variable    = "BURDEN",
    transformation_y      = "",
    family_test           = "wilcoxon",
    samples_sql_condition = "",
    key                   = key
  )

  ssEnv      <- SEMseeker:::core_get_session_info()
  chart_root <- file.path(ssEnv$result_folderChart, "COMPARISON", "")
  png_files  <- list.files(chart_root, pattern = "\\.png$", recursive = TRUE)
  expect_gte(length(png_files), 2)  # boxplot + violin
})

test_that("plot_box_plot: t.test family also creates PNG files", {
  skip_if_not_installed("ggpubr")

  tf  <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  set.seed(2)
  df <- data.frame(
    BURDEN = c(stats::rnorm(15, 0), stats::rnorm(15, 5)),
    GROUP  = factor(c(rep("ctrl", 15), rep("case", 15)))
  )
  key <- list(AREA = "GENE", SUBAREA = "TSS200", MARKER = "MUTATIONS", FIGURE = "K850")

  SEMseeker:::plot_box_plot(
    dataFrameToPlot       = df,
    independent_variable  = "GROUP",
    dependent_variable    = "BURDEN",
    transformation_y      = "",
    family_test           = "t.test",
    samples_sql_condition = "",
    key                   = key
  )

  ssEnv      <- SEMseeker:::core_get_session_info()
  chart_root <- file.path(ssEnv$result_folderChart, "COMPARISON", "")
  png_files  <- list.files(chart_root, pattern = "\\.png$", recursive = TRUE)
  expect_gte(length(png_files), 2)
})

test_that("plot_box_plot: non-dichotomous family (pearson) is skipped silently", {
  skip_if_not_installed("ggpubr")

  tf  <- sem_test_folder()
  SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df  <- data.frame(BURDEN = 1:10, GROUP = factor(rep("ctrl", 10)))
  key <- list(AREA = "GENE", SUBAREA = "TSS200", MARKER = "MUTATIONS", FIGURE = "K850")

  # Should return NULL invisibly without creating any files
  result <- SEMseeker:::plot_box_plot(
    dataFrameToPlot       = df,
    independent_variable  = "GROUP",
    dependent_variable    = "BURDEN",
    transformation_y      = "",
    family_test           = "pearson",
    samples_sql_condition = "",
    key                   = key
  )

  ssEnv      <- SEMseeker:::core_get_session_info()
  chart_root <- file.path(ssEnv$result_folderChart, "COMPARISON", "")
  png_files  <- list.files(chart_root, pattern = "\\.png$", recursive = TRUE)
  expect_equal(length(png_files), 0)
})
