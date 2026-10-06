# 2026-06-09: documents the SAFE PATTERN for extracting schema
# and row count when staying lazy through position_pivot_to_probe.
#
# Background — the bug this test guards against:
#
#   Calling `collect_schema()` or `select(pl$len())$collect()` on the
#   LazyFrame returned by `anno_position_pivot_to_probe(signal_pivot)` FORCES
#   Polars to execute the underlying inner-join on (CHR, START, END).
#   Polars can NOT infer the post-join schema or row count without
#   running the join (an inner join may drop or duplicate rows depending
#   on the right-hand side, so cardinality is data-dependent).
#
#   On ewas-scale (~367k probes × ~4k sample columns) the join allocates
#   ~12 GB of Rust heap to materialise the intermediate join product,
#   pushing the R process over the macOS jetsam threshold. This was the
#   silent OOM kill observed in v18, v21, and v25–v30.
#
#   Sample columns are IDENTICAL between the POSITION pivot and the
#   PROBE pivot (only the key column differs: CHR/START/END → PROBE),
#   so taking schema/count from the RAW LazyFrame (which is a direct
#   scan_parquet) is both equivalent and cheap (parquet footer = O(1)
#   metadata read).
#
# Pattern (analyze_batch.R:56-69):
#
#   schema_cols_position <- names(signal_pivot$collect_schema())   # FAST
#   n_probes <- as.integer(signal_pivot$select(pl$len())$collect()$
#                            to_data_frame()$len[1])                # FAST
#   sample_cols <- setdiff(schema_cols_position,
#                          c("CHR","START","END","PROBE"))
#   signal_lazy <- anno_position_pivot_to_probe(signal_pivot)            # lazy

test_that("sample columns extracted from raw POSITION pivot match the PROBE pivot post-join", {

  skip_on_cran()

  tempFolder <- tempfile("ppp_schema_pattern_")
  dir.create(file.path(tempFolder, "Data"), recursive = TRUE)
  ssEnv <- SEMseeker:::core_init_env(tempFolder,
                                parallel_strategy = "sequential",
                                tech = "K850",
                                iqrTimes = 3, verbosity = 1)
  on.exit({ SEMseeker:::core_close_env(); unlink(tempFolder, recursive = TRUE) },
          add = TRUE)

  anno <- data.frame(
    PROBE = paste0("cg", sprintf("%08d", 1:10)),
    CHR   = "1",
    START = seq(100L, 1000L, by = 100L),
    END   = seq(100L, 1000L, by = 100L),
    K850  = TRUE,
    stringsAsFactors = FALSE
  )
  sample_names <- paste0("S", sprintf("%03d", 1:25))
  pivot_df <- data.frame(
    CHR   = anno$CHR,
    START = anno$START,
    END   = anno$END,
    stringsAsFactors = FALSE
  )
  for (s in sample_names) pivot_df[[s]] <- runif(nrow(anno))
  pivot_lazy <- polars::as_polars_df(pivot_df)$lazy()

  testthat::local_mocked_bindings(
    anno_probe_annotation_build = function(tech) anno,
    .package = "SEMseeker"
  )

  # SAFE PATTERN: schema/count taken BEFORE the join.
  schema_cols_position <- names(pivot_lazy$collect_schema())
  n_probes <- as.integer(
    as.data.frame(pivot_lazy$select(polars::pl$len())$collect())$len[1]
  )
  sample_cols_pre <- setdiff(schema_cols_position,
                             c("CHR", "START", "END", "PROBE"))

  signal_lazy <- SEMseeker:::anno_position_pivot_to_probe(pivot_lazy)

  # Equivalence: sample columns from raw pivot match sample columns from
  # the PROBE pivot (post-join).
  schema_cols_probe <- names(signal_lazy$collect_schema())
  sample_cols_post  <- setdiff(schema_cols_probe,
                               c("CHR", "START", "END", "PROBE"))
  testthat::expect_setequal(sample_cols_pre, sample_cols_post)
  testthat::expect_setequal(sample_cols_pre, sample_names)

  # n_probes from raw POSITION pivot = n rows post anno_position_pivot_to_probe
  # (no probes dropped by the K850 filter in this synthetic data).
  collected <- as.data.frame(signal_lazy$collect())
  testthat::expect_equal(n_probes, nrow(collected))
  testthat::expect_equal(n_probes, nrow(anno))
})

test_that("schema extraction from raw POSITION pivot does NOT depend on anno_position_pivot_to_probe call", {

  skip_on_cran()

  tempFolder <- tempfile("ppp_schema_independence_")
  dir.create(file.path(tempFolder, "Data"), recursive = TRUE)
  ssEnv <- SEMseeker:::core_init_env(tempFolder,
                                parallel_strategy = "sequential",
                                tech = "K850",
                                iqrTimes = 3, verbosity = 1)
  on.exit({ SEMseeker:::core_close_env(); unlink(tempFolder, recursive = TRUE) },
          add = TRUE)

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
    S001  = c(0.1, 0.2, 0.3, 0.4, 0.5),
    S002  = c(0.6, 0.7, 0.8, 0.9, 1.0),
    stringsAsFactors = FALSE
  )
  pivot_lazy <- polars::as_polars_df(pivot_df)$lazy()

  testthat::local_mocked_bindings(
    anno_probe_annotation_build = function(tech) anno,
    .package = "SEMseeker"
  )

  # Take schema/count from raw pivot first.
  schema_cols <- names(pivot_lazy$collect_schema())
  testthat::expect_setequal(schema_cols,
                            c("CHR", "START", "END", "S001", "S002"))
  testthat::expect_false("PROBE" %in% schema_cols)

  # Now run anno_position_pivot_to_probe — the resulting LazyFrame should be
  # internally consistent without forcing materialisation of the raw
  # pivot's schema again.
  signal_lazy <- SEMseeker:::anno_position_pivot_to_probe(pivot_lazy)
  testthat::expect_s3_class(signal_lazy, "polars_lazy_frame")

  # Resulting collect_schema reports PROBE as the new key.
  schema_probe <- names(signal_lazy$collect_schema())
  testthat::expect_true("PROBE" %in% schema_probe)
  testthat::expect_false("CHR" %in% schema_probe)
})
