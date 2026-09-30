## Coverage for io_list_bed_files_for_marker_figure() (previously untested).
##
## Lists per-sample BED/bedGraph files under
##   <result_folderData>/<sample_group>/<MARKER>_<FIGURE>/...
## for a given marker+figure.

test_that("io_list_bed_files_for_marker_figure finds only matching marker/figure beds", {
  tf <- sem_test_folder()
  ss <- SEMseeker:::core_init_env(result_folder = tf, start_fresh = TRUE)
  on.exit({ SEMseeker:::core_close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  data_root <- ss$result_folderData

  # matching: Case / MUTATIONS_HYPER
  d_match <- file.path(data_root, "Case", "MUTATIONS_HYPER")
  dir.create(d_match, recursive = TRUE, showWarnings = FALSE)
  writeLines("chr1\t1\t2", file.path(d_match, "s1.bed"))

  # non-matching marker/figure
  d_other <- file.path(data_root, "Case", "DELTAS_HYPO")
  dir.create(d_other, recursive = TRUE, showWarnings = FALSE)
  writeLines("chr1\t1\t2", file.path(d_other, "s2.bed"))

  hits <- SEMseeker:::io_list_bed_files_for_marker_figure("MUTATIONS", "HYPER")
  expect_length(hits, 1L)
  expect_match(hits, "MUTATIONS_HYPER")
  expect_false(any(grepl("DELTAS_HYPO", hits)))

  # a marker/figure with no files returns an empty character vector
  expect_length(
    SEMseeker:::io_list_bed_files_for_marker_figure("MUTATIONS", "BOTH"),
    0L)
})
