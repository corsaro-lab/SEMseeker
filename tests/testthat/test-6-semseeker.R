test_that("semeeker", {

  tempFolder <- sem_test_folder()
  unlink(tempFolder,recursive = TRUE)

  ####################################################################################

  # "sequential" so the test works under devtools::load_all() too. The
  # "multisession" strategy from setup.R only resolves SEMseeker::: against
  # the INSTALLED package, which makes load_all-time test runs fail on
  # post-merge signature changes (e.g. AI-075 io_bed_file_name(skip_dir_create)).
  # Same rationale as test-7-association_analysis.R.
  SEMseeker::semseeker(
    input         = signal_data,
    sample_sheet  = mySampleSheet,
    result_folder = tempFolder,
    parallel_strategy = "sequential"
  )

  # semseeker() closes the session it opened, so reading its artefacts needs the
  # session back. This run named no areas and no markers, so the reopen names none
  # either and both land on the same defaults.
  SEMseeker:::core_init_env(result_folder = tempFolder,
                            parallel_strategy = "sequential",
                            start_fresh = FALSE)

  ssEnv <- SEMseeker:::core_get_session_info()
  keys <- subset(ssEnv$keys_areas_subareas_markers_figures)
  # core_name_cleaning uppercases Sample_ID inside semseeker(); use the same for comparison
  cleaned_sample_ids <- SEMseeker:::core_name_cleaning(mySampleSheet$Sample_ID)

  for (k in seq_len(nrow(keys)))
  {
    # k <- 1
    key <- keys[k,]
    marker <- as.character(key$MARKER)
    figure <- as.character(key$FIGURE)
    area <- as.character(key$AREA)
    subarea <- as.character(key$SUBAREA)

    mutations_pivot_file_name <- SEMseeker:::io_pivot_file_name_parquet("MUTATIONS",figure,area,subarea)
    # NOTE: silent next preserved on purpose. The iteration builds
    # MUTATIONS_<figure> for EVERY figure in keys (including SIGNAL_MEAN),
    # so a hard expect_true(file.exists(...)) here fires false positives for
    # combos that are not meant to exist (e.g. MUTATIONS_MEAN). Fixing the
    # iteration to filter only (HYPER, HYPO) figures + comparable markers
    # is tracked separately — see AI-090.
    if(file.exists(mutations_pivot_file_name))
      mutations_pivot <- as.data.frame(polars::pl$read_parquet(mutations_pivot_file_name))
    else
      next

    io_pivot_file_name <- SEMseeker:::io_pivot_file_name_parquet(marker,figure,area,subarea)
    # derived markers (LESIONS, DELTA*) may not exist when mutations are too
    # sparse for that figure × area × subarea combo — see AI-090 to rework
    # the iteration so we can hard-assert presence on the combos that ARE
    # supposed to exist.
    if(!file.exists(io_pivot_file_name))
      next
    pivot <- as.data.frame(polars::pl$read_parquet(io_pivot_file_name))

    pivot <- pivot[,-c(1:3)]
    mutations_pivot <- mutations_pivot[,-c(1:3)]
    testthat::expect_true(nrow(pivot)<=nprobes)
    testthat::expect_true(nrow(pivot)>0)

    pivot <- pivot[,order(colnames(pivot))]
    mutations_pivot <- mutations_pivot[,order(colnames(mutations_pivot))]

    testthat::expect_true(all(colnames(pivot) %in% cleaned_sample_ids))
    testthat::expect_true(all(colnames(pivot) == colnames(mutations_pivot)))

    testthat::expect_true(ncol(pivot)==ncol(mutations_pivot))
    testthat::expect_true(nrow(pivot)==nrow(mutations_pivot))

    pivot[is.na(pivot)] <- 0
    mutations_pivot[is.na(mutations_pivot)] <- 0

    pivot[pivot > 0] <- 1
    mutations_pivot[mutations_pivot > 0] <- 1

    # Only MUTATIONS pivot is identical to itself; derived markers (DELTA*, LESIONS)
    # may have signed/windowed values that don't match the binary MUTATIONS mask
    if (marker == "MUTATIONS")
    {
      if (!all(as.data.frame(pivot) == as.data.frame(mutations_pivot)))
      {
        print(marker)
        print(figure)
      }
      testthat::expect_true(all(as.data.frame(pivot) == as.data.frame(mutations_pivot)))
    }

    # Per-sample BED presence is NOT a valid invariant for derived markers.
    # By design (long-reads readiness) io_dump_sample_as_bed_file() skips
    # writing when the per-sample data has 0 rows — DELTAS/DELTAR/LESIONS
    # can legitimately be empty for a given (sample, figure) even when
    # MUTATIONS is not. The pivot-level expectations above (cols match,
    # nrow match, sample IDs match) are the authoritative checks.

  }

  ####################################################################################
  SEMseeker:::core_close_env()
  unlink(tempFolder,recursive = TRUE)
})

