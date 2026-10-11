# 2026-06-09: anno_position_pivot_to_probe returns LazyFrame.
#
# Contract change: returns polars_lazy_frame (was R data.frame in the
# legacy implementation). Caller is responsible for any materialization,
# and most callers should NOT materialize — sem_analyze_batch resume path
# now consumes the LazyFrame end-to-end. Prior R-side subset/sort
# materialised 4–5 copies of the input, blowing 80+ GB peak on
# ~367k × 4k inputs and triggering macOS jetsam OOM-kill silently.

test_that("anno_position_pivot_to_probe accepts a polars LazyFrame and returns a LazyFrame", {

  skip_on_cran()

  tempFolder <- tempfile("ppp_lazy_")
  dir.create(file.path(tempFolder, "Data"), recursive = TRUE)
  ssEnv <- SEMseeker:::core_init_env(tempFolder,
                                parallel_strategy = "sequential",
                                tech = "K850",
                                iqrTimes = 3, verbosity = 1)
  on.exit({ SEMseeker:::core_close_env(); unlink(tempFolder, recursive = TRUE) },
          add = TRUE)

  # Tiny synthetic annotation: 5 probes on chr1, two samples worth of values.
  anno <- data.frame(
    PROBE = paste0("cg", sprintf("%08d", 1:5)),
    CHR   = "1",
    START = c(100L, 200L, 300L, 400L, 500L),
    END   = c(100L, 200L, 300L, 400L, 500L),
    K850  = TRUE,
    stringsAsFactors = FALSE
  )

  pivot_df <- data.frame(
    CHR   = anno$CHR,
    START = anno$START,
    END   = anno$END,
    S001  = c(0.10, 0.20, 0.30, 0.40, 0.50),
    S002  = c(0.15, 0.25, 0.35, 0.45, 0.55),
    stringsAsFactors = FALSE
  )
  pivot_lazy <- polars::as_polars_df(pivot_df)$lazy()
  testthat::expect_s3_class(pivot_lazy, "polars_lazy_frame")

  # Stub anno_probe_annotation_build so the test does not require the Bioconductor
  # annotation package (macOS setup.R skips loading it).
  testthat::local_mocked_bindings(
    anno_probe_annotation_build = function(tech) anno,
    .package = "SEMseeker"
  )

  out <- SEMseeker:::anno_position_pivot_to_probe(pivot_lazy)

  # NEW CONTRACT: lazy passthrough: return type is LazyFrame.
  testthat::expect_s3_class(out, "polars_lazy_frame")
  testthat::expect_false(inherits(out, "data.frame"))

  # Schema: PROBE column + sample columns; CHR/START/END/tech dropped.
  schema_cols <- names(out$collect_schema())
  testthat::expect_true("PROBE" %in% schema_cols)
  testthat::expect_false("CHR"   %in% schema_cols)
  testthat::expect_false("START" %in% schema_cols)
  testthat::expect_false("END"   %in% schema_cols)
  testthat::expect_false("K850"  %in% schema_cols)
  testthat::expect_true("S001" %in% schema_cols)
  testthat::expect_true("S002" %in% schema_cols)

  # Values preserved on collect (compare by PROBE-key order-independently).
  collected <- as.data.frame(out$collect())
  rownames(collected) <- collected$PROBE
  testthat::expect_setequal(collected$PROBE, anno$PROBE)
  testthat::expect_equal(collected[anno$PROBE, "S001"], pivot_df$S001)
  testthat::expect_equal(collected[anno$PROBE, "S002"], pivot_df$S002)
})

test_that("anno_position_pivot_to_probe still accepts an R data.frame input", {

  skip_on_cran()

  tempFolder <- tempfile("ppp_df_")
  dir.create(file.path(tempFolder, "Data"), recursive = TRUE)
  ssEnv <- SEMseeker:::core_init_env(tempFolder,
                                parallel_strategy = "sequential",
                                tech = "K850",
                                iqrTimes = 3, verbosity = 1)
  on.exit({ SEMseeker:::core_close_env(); unlink(tempFolder, recursive = TRUE) },
          add = TRUE)

  anno <- data.frame(
    PROBE = paste0("cg", sprintf("%08d", 1:3)),
    CHR   = "1",
    START = c(100L, 200L, 300L),
    END   = c(100L, 200L, 300L),
    K850  = TRUE,
    stringsAsFactors = FALSE
  )
  pivot_df <- data.frame(
    CHR = anno$CHR, START = anno$START, END = anno$END,
    S001 = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(
    anno_probe_annotation_build = function(tech) anno,
    .package = "SEMseeker"
  )

  out <- SEMseeker:::anno_position_pivot_to_probe(pivot_df)
  testthat::expect_s3_class(out, "polars_lazy_frame")
  collected <- as.data.frame(out$collect())
  testthat::expect_setequal(collected$PROBE, anno$PROBE)
  testthat::expect_setequal(setdiff(colnames(collected), "PROBE"), "S001")
})

test_that("anno_position_pivot_to_probe drops probes flagged FALSE for the tech", {

  skip_on_cran()

  tempFolder <- tempfile("ppp_filter_")
  dir.create(file.path(tempFolder, "Data"), recursive = TRUE)
  ssEnv <- SEMseeker:::core_init_env(tempFolder,
                                parallel_strategy = "sequential",
                                tech = "K850",
                                iqrTimes = 3, verbosity = 1)
  on.exit({ SEMseeker:::core_close_env(); unlink(tempFolder, recursive = TRUE) },
          add = TRUE)

  # Two probes ok for K850, one NOT on K850 (should be dropped).
  anno <- data.frame(
    PROBE = paste0("cg", sprintf("%08d", 1:3)),
    CHR   = "1",
    START = c(100L, 200L, 300L),
    END   = c(100L, 200L, 300L),
    K850  = c(TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  pivot_df <- data.frame(
    CHR = anno$CHR, START = anno$START, END = anno$END,
    S001 = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(
    anno_probe_annotation_build = function(tech) anno,
    .package = "SEMseeker"
  )

  out <- SEMseeker:::anno_position_pivot_to_probe(polars::as_polars_df(pivot_df)$lazy())
  collected <- as.data.frame(out$collect())
  testthat::expect_setequal(collected$PROBE, anno$PROBE[anno$K850])
  testthat::expect_equal(nrow(collected), 2L)
})
