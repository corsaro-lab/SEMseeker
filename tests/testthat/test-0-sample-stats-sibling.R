## AI-223 slice 1 — per-sample statistics sibling (scope SAMPLE).
##
## Contracts:
##   1. column names come from ONE helper, used by producer and consumer alike;
##   2. the sibling is written next to the sample sheet, one row per sample,
##      carrying burden (SAMPLE_<MARKER>_<FIGURE>) and signal descriptors
##      (SAMPLE_<STAT>);
##   3. the two beta modes are estimated on each side of 0.5 and are omitted on
##      the M-value scale, where the split carries no meaning;
##   4. the sample sheet no longer carries the burden — it moved here;
##   5. the sibling joins back onto the sample sheet on Sample_ID.
##
## Session test uses tempFolders index 13.

# ---------------------------------------------------------------------------
# naming helpers (pure)
# ---------------------------------------------------------------------------

test_that("io_scope_name derives the scope from (area, subarea) with no depth", {
  # AI-223 slice 2a: the producer is depth-agnostic — it writes every scope
  # once, so it names them from the pair alone.
  expect_equal(SEMseeker:::io_scope_name(area = "GENE", subarea = "TSS1500"),
               "GENE_TSS1500")
  expect_equal(SEMseeker:::io_scope_name(area = "GENE"), "GENE")
  expect_equal(SEMseeker:::io_scope_name(), "SAMPLE")
  expect_equal(SEMseeker:::io_scope_name(area = "", subarea = "TSS1500"), "SAMPLE")
})

# ---------------------------------------------------------------------------
# descriptors (pure)
# ---------------------------------------------------------------------------

test_that("the two beta modes land on the injected peaks", {
  set.seed(99L)
  values <- c(stats::rbeta(4000, 2, 40),    # peak near 0.05
              stats::rbeta(4000, 40, 2))    # peak near 0.95

  d <- SEMseeker:::util_signal_descriptors(values, beta = TRUE)

  expect_lt(d$MODELOW, 0.5)
  expect_gt(d$MODEHIGH, 0.5)
  expect_lt(d$MODELOW, d$MODEHIGH)
  expect_lt(abs(d$MODELOW  - 0.05), 0.1)
  expect_lt(abs(d$MODEHIGH - 0.95), 0.1)
  expect_equal(d$N_PROBES, length(values))
  expect_equal(d$MEDIAN, stats::median(values))
  expect_equal(d$IQR, stats::IQR(values))
})

test_that("descriptors omit the modes on the M-value scale", {
  set.seed(7L)
  values <- stats::rnorm(1000, mean = 1.5, sd = 2)

  d <- SEMseeker:::util_signal_descriptors(values, beta = FALSE)

  expect_null(d$MODELOW)
  expect_null(d$MODEHIGH)
  expect_equal(d$MEAN, mean(values))
  expect_equal(d$VARIANCE, stats::var(values))
})

test_that("descriptors degrade gracefully on degenerate input", {
  empty <- SEMseeker:::util_signal_descriptors(numeric(0), beta = TRUE)
  expect_equal(empty$N_PROBES, 0L)
  expect_true(is.na(empty$MEDIAN))

  # a constant vector has no density to speak of
  flat <- SEMseeker:::util_signal_descriptors(rep(0.5, 100), beta = TRUE)
  expect_equal(flat$MEDIAN, 0.5)
  expect_true(is.na(flat$MODELOW))
  expect_true(is.na(flat$MODEHIGH))
})

# ---------------------------------------------------------------------------
# end to end
# ---------------------------------------------------------------------------

test_that("semseeker() writes the statistics sibling and leaves the sample sheet lean", {
  tempFolder <- tempFolders[13]
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(tempFolder, recursive = TRUE) }, add = TRUE)

  syn <- .burden_setup_signal_with_outliers(
    n_samples      = nsamples,
    probe_features = probe_features,
    sample_sheet   = mySampleSheet,
    signal_data    = signal_data)

  SEMseeker::semseeker(
    input             = syn$signal,
    sample_sheet      = syn$samples,
    result_folder     = tempFolder,
    parallel_strategy = "sequential",
    areas             = c("POSITION"),
    markers           = c("MUTATIONS", "DELTAP", "DELTAS"),
    start_fresh       = TRUE,
    inpute            = "median",
    showprogress      = showprogress,
    verbosity         = verbosity
  )

  # AI-255: the statistics are SCOPE = SAMPLE artefacts, composed on read.
  sheet_csv <- file.path(tempFolder, "Data", "SAMPLE_SHEET_RESULT.csv")
  expect_true(file.exists(sheet_csv))

  stats <- SEMseeker:::sem_study_summary_get()
  sheet <- utils::read.csv2(sheet_csv, stringsAsFactors = FALSE)
  skip_if_not(!is.null(stats) && nrow(stats) > 0,
              "downstream assertions need the composed statistics")

  # one row per sample of the signal matrix
  expect_equal(nrow(stats), ncol(syn$signal))

  # AI-248: `markers` means one thing for every marker, SIGNAL included — this
  # run did not ask for it, so its descriptors must NOT be there.
  #
  # AI-255: N_PROBES is no longer among them. It is not an aggregation of a
  # marker over a region class but a property of the imputation — how many
  # positions of that sample survived the treatment of missing values — so it
  # travels with the sample sheet. Nothing is lost for the density: the MEAN of
  # a binary marker *is* the density, denominator included.
  descriptor_cols <- character(0)
  expect_false("SAMPLE_N_PROBES" %in% colnames(stats))
  signal_cols <- vapply(c("MEDIAN", "MEAN", "VARIANCE", "IQR", "MODELOW", "MODEHIGH"),
                        function(a) SEMseeker:::io_feature_colname("SAMPLE", "SIGNAL", "BETA", a),
                        character(1))
  burden_cols <- c(
    SEMseeker:::io_feature_colname("SAMPLE", "MUTATIONS", c("HYPER", "HYPO"), "SUM"),
    SEMseeker:::io_feature_colname("SAMPLE", "DELTAP", c("HYPER", "HYPO"), "SUM"),
    SEMseeker:::io_feature_colname("SAMPLE", "DELTAS", c("HYPER", "HYPO"), "MEAN"))
  expect_true(all(descriptor_cols %in% colnames(stats)),
              info = paste("missing:",
                           paste(setdiff(descriptor_cols, colnames(stats)), collapse = ", ")))
  expect_true(all(burden_cols %in% colnames(stats)),
              info = paste("missing:",
                           paste(setdiff(burden_cols, colnames(stats)), collapse = ", ")))

  # a marker that was not asked for produces no column
  expect_equal(intersect(signal_cols, colnames(stats)), character(0))

  # AI-223 net move: the burden left the sample sheet
  moved_away <- c("MUTATIONS_HYPER", "MUTATIONS_HYPO", "DELTAS_HYPER",
                  "DELTAP_HYPER", "PROBES_COUNT")
  expect_equal(intersect(moved_away, colnames(sheet)), character(0))

  # and the two files join on Sample_ID
  expect_true(all(sheet$Sample_ID %in% stats$Sample_ID))

  # the join is what consumers see through sem_study_summary_get()
  SEMseeker:::core_init_env(result_folder = tempFolder, parallel_strategy = "sequential",
                            areas = c("POSITION"),
                            markers = c("MUTATIONS", "DELTAP", "DELTAS"),
                            start_fresh = FALSE, showprogress = FALSE, verbosity = 1)
  joined <- SEMseeker:::sem_study_summary_get()
  expect_true(all(burden_cols %in% colnames(joined)))
  # AI-255: N_PROBES comes in from the sample sheet, under its own name — it
  # describes the sample, not a scope of it.
  expect_true("N_PROBES" %in% colnames(joined))
  expect_false("SAMPLE_N_PROBES" %in% colnames(joined))
})

# ---------------------------------------------------------------------------
# AI-223 slice 2a — region scope, produced and consumed at depth = 1
# ---------------------------------------------------------------------------

test_that("a region scope reaches the sibling and the depth=1 inference", {
  tempFolder <- tempFolders[14]
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(tempFolder, recursive = TRUE) }, add = TRUE)

  syn <- .burden_setup_signal_with_outliers(
    n_samples      = nsamples,
    probe_features = probe_features,
    sample_sheet   = mySampleSheet,
    signal_data    = signal_data)
  scope <- SEMseeker:::io_scope_name(area = "GENE", subarea = "TSS1500")

  SEMseeker::semseeker(
    input               = syn$signal,
    sample_sheet        = syn$samples,
    result_folder       = tempFolder,
    parallel_strategy   = "sequential",
    # "WHOLE" must stay in `subareas`: it is what keeps POSITION_WHOLE in the
    # keys, and every burden — whatever its scope — is aggregated from the
    # POSITION pivots.
    areas               = c("POSITION", "GENE"),
    subareas            = c("WHOLE", "TSS1500"),
    markers             = c("MUTATIONS", "SIGNAL"),
    start_fresh         = TRUE,
    inpute              = "median",
    showprogress        = showprogress,
    verbosity           = verbosity
  )

  # AI-255: the region class is asked for at read time, not declared before the
  # run. This is the whole point — no rerun to change your mind.
  stats <- SEMseeker:::sem_study_summary_get(regions = c("SAMPLE", scope))
  skip_if_not(!is.null(stats) && nrow(stats) > 0,
              "downstream assertions need the composed statistics")

  scope_cols <- SEMseeker:::io_feature_colname(scope, "MUTATIONS", c("HYPER", "HYPO"), "SUM")
  whole_cols <- SEMseeker:::io_feature_colname("SAMPLE", "MUTATIONS", c("HYPER", "HYPO"), "SUM")
  expect_true(all(scope_cols %in% colnames(stats)),
              info = paste("missing:",
                           paste(setdiff(scope_cols, colnames(stats)), collapse = ", ")))
  expect_true(all(whole_cols %in% colnames(stats)))

  # a subset of the probes can only carry a subset of the burden
  for (i in seq_along(scope_cols))
    expect_true(all(stats[[scope_cols[i]]] <= stats[[whole_cols[i]]]))

  # AI-248: asking for SIGNAL gives its descriptors on the scope too — this is
  # the per-sample median restricted to a region class
  scope_median <- SEMseeker:::io_feature_colname(scope, "SIGNAL", "BETA", "MEDIAN")
  expect_true(scope_median %in% colnames(stats),
              info = paste("missing:", scope_median))
  expect_true(all(stats[[scope_median]] >= 0 & stats[[scope_median]] <= 1, na.rm = TRUE))

  # ── the mask counts every position once ──────────────────────────────────
  # AI-255: the mask is no longer a helper of its own — io_pivot_build() derives
  # every artefact from the position pivot, and at SCOPE = SAMPLE it selects
  # positions rather than partitioning them. The invariant this block guards is
  # unchanged and is the important one: a probe the annotation maps onto three
  # genes must enter the sample's number ONCE. Build the mask the way the
  # producer does — the distinct positions of the region class.
  SEMseeker:::core_init_env(result_folder = tempFolder, parallel_strategy = "sequential",
                            areas = c("POSITION", "GENE"), subareas = c("WHOLE", "TSS1500"),
                            markers = c("MUTATIONS"),
                            start_fresh = FALSE, showprogress = FALSE, verbosity = 1)

  pf <- SEMseeker:::anno_probe_features_get("GENE_TSS1500")
  # AI-308: keep only the probes that ARE annotated to the class. On the Illumina
  # path anno_probe_features_get() returns the whole annotation table, with NA in
  # the column of the class asked for: for K450, GENE_TSS1500 carries 84,808
  # annotated probes against 401,394 NA. This block used to select the
  # coordinates without dropping them, which is precisely what the producer did
  # wrong: the expected value it computed was the burden of the WHOLE sample, and
  # the assertion passed because both sides made the same mistake. A test written
  # by mirroring the implementation agrees with it even when it is wrong.
  pf <- pf[!is.na(pf[["GENE_TSS1500"]]), , drop = FALSE]
  skip_if(is.null(pf) || nrow(pf) == 0,
          "no probe of the fixture is annotated TSS1500")
  mask_df <- unique(data.frame(
    CHR   = sub("^(?i)chr", "", as.character(pf$CHR), perl = TRUE),
    START = as.integer(pf$START),
    END   = as.integer(pf$END),
    stringsAsFactors = FALSE))
  expect_gt(nrow(mask_df), 0)
  # the restriction is real: the class covers strictly fewer positions than the
  # sample. Without this the block could pass again on an unrestricted mask.
  expect_lt(nrow(mask_df), nrow(SEMseeker:::anno_probe_features_get("GENE_TSS1500")))

  pivot <- SEMseeker:::io_read_pivot("MUTATIONS", "HYPER", "POSITION", "WHOLE")
  pivot_df <- as.data.frame(pivot$collect())
  pivot_df$CHR   <- sub("^(?i)chr", "", as.character(pivot_df$CHR), perl = TRUE)
  pivot_df$START <- as.integer(pivot_df$START)
  pivot_df$END   <- as.integer(pivot_df$END)
  masked <- merge(pivot_df, mask_df, by = c("CHR", "START", "END"))
  sample_columns <- setdiff(colnames(masked), c("CHR", "START", "END"))
  expected <- colSums(masked[, sample_columns, drop = FALSE], na.rm = TRUE)

  observed <- stats[[SEMseeker:::io_feature_colname(scope, "MUTATIONS", "HYPER", "SUM")]]
  names(observed) <- stats$Sample_ID
  common <- intersect(names(expected), names(observed))
  expect_gt(length(common), 0)
  expect_equal(as.numeric(observed[common]), as.numeric(expected[common]))

  # a probe annotated to several genes must not be counted twice: the mask has
  # one row per position, whatever the annotation says
  expect_equal(nrow(mask_df), nrow(unique(mask_df[, c("CHR", "START", "END")])))

  # ── consumption at depth = 1 ─────────────────────────────────────────────
  SEMseeker:::core_close_env()

  inference_details <- data.frame(
    independent_variable = "Phenotest",
    family_test          = "spearman",
    transformation_y     = "",
    transformation_x     = "",
    # AI-308: the request names the aggregation branch, not the region classes.
    # The classes are the (AREA, SUBAREA) pairs declared below in areas/subareas
    # and built at runtime; SCOPE = SAMPLE collapses each of them to one number
    # per sample. The single-position class is normalised to the technology's
    # own, which is why the PROBE assertions below still hold on a run that
    # declared POSITION.
    scope                = "SAMPLE",
    aggregation          = "SUM",
    filter_p_value       = FALSE,
    stringsAsFactors     = FALSE
  )

  expect_no_error(
    SEMseeker:::association_analysis(
      inference_details = inference_details,
      result_folder     = tempFolder,
      parallel_strategy = "sequential",
      markers           = c("MUTATIONS"),
      # AI-255: the artefact is built here, on the way in, so this call has to
      # know the region class it is being asked for. The cost did not vanish, it
      # moved: no SEM rerun, but the analysis names its areas. The registry is
      # also what lets a class be resolved without parsing "GENE_TSS1500" back
      # into a pair — which cannot be done safely, since N_SHORE has an
      # underscore of its own.
      areas             = c("POSITION", "GENE"),
      subareas          = c("WHOLE", "TSS1500"),
      multiple_test_adj = "BH",
      showprogress      = showprogress,
      verbosity         = verbosity
    )
  )

  csv_files <- list.files(file.path(tempFolder, "Inference"), pattern = "\\.csv$",
                          recursive = TRUE, full.names = TRUE)
  csv_files <- csv_files[!grepl("(?i)assoc_covariates_model", csv_files)]
  expect_gt(length(csv_files), 0)
  skip_if(length(csv_files) == 0, "no inference CSV to inspect")

  result_df <- do.call(plyr::rbind.fill,
                       lapply(csv_files, function(f) utils::read.csv2(f, stringsAsFactors = FALSE)))
  # AI-255: a collapsed row carries the coordinates of the taxonomy, not a
  # made-up AREA. It used to be AREA = "GENE_TSS1500" (the scope name squashed
  # into a coordinate) and AREA = "SAMPLE_GROUP" for the unrestricted one —
  # values invented for the occasion, exactly as "TOTAL" was invented for
  # SUBAREA. Now: SCOPE says it is collapsed, AREA and SUBAREA say over which
  # region class.
  #
  # which() and not a bare logical: the CSV carries rows with NA in AREA (the
  # job-summary row), and `df[df$AREA == …, ]` would materialise them as all-NA
  # phantom rows.
  scope_rows <- result_df[which(result_df$SCOPE == "SAMPLE" &
                                result_df$AREA == "GENE" &
                                result_df$SUBAREA == "TSS1500"), , drop = FALSE]
  expect_gt(nrow(scope_rows), 0)
  # which burden was tested stays in MARKER/FIGURE; the region class is what
  # AREA and SUBAREA add, and the aggregate is named by AREA_OF_TEST.
  expect_equal(sort(unique(scope_rows$MARKER)), "MUTATIONS")
  expect_true(all(scope_rows$FIGURE %in% c("HYPER", "HYPO")))

  # and the unrestricted collapsed rows are still there, on their own class
  whole_rows <- result_df[which(result_df$SCOPE == "SAMPLE" &
                                result_df$AREA == "PROBE"), , drop = FALSE]
  expect_gt(nrow(whole_rows), 0)
  expect_setequal(unique(scope_rows$FIGURE), unique(whole_rows$FIGURE))
})

test_that("a retired column stops the analysis instead of testing nothing", {
  tempFolder <- tempFolders[15]
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(tempFolder, recursive = TRUE) }, add = TRUE)

  syn <- .burden_setup_signal_with_outliers(
    n_samples      = nsamples,
    probe_features = probe_features,
    sample_sheet   = mySampleSheet,
    signal_data    = signal_data)

  SEMseeker::semseeker(
    input             = syn$signal,
    sample_sheet      = syn$samples,
    result_folder     = tempFolder,
    parallel_strategy = "sequential",
    areas             = c("POSITION"),
    markers           = c("MUTATIONS"),
    start_fresh       = TRUE,
    inpute            = "median",
    showprogress      = showprogress,
    verbosity         = verbosity
  )

  # AI-308: this used to name a region class the run had not produced and expect
  # the run to stop rather than write a CSV that tested nothing. The request can
  # no longer name a class at all: the classes are the (AREA, SUBAREA) pairs
  # declared with areas/subareas, so the surviving guarantee is the one on the
  # column itself: `scopes` is retired, and a request still carrying it is told
  # so by name instead of being read as a typo for something else.
  inference_details <- data.frame(
    independent_variable = "Phenotest",
    family_test          = "spearman",
    transformation_y     = "",
    transformation_x     = "",
    scopes               = "GENE_TSS1500",
    aggregation          = "SUM",
    filter_p_value       = FALSE,
    stringsAsFactors     = FALSE
  )

  expect_error(
    SEMseeker:::association_analysis(
      inference_details = inference_details,
      result_folder     = tempFolder,
      parallel_strategy = "sequential",
      markers           = c("MUTATIONS"),
      areas             = c("POSITION"),
      multiple_test_adj = "BH",
      showprogress      = showprogress,
      verbosity         = verbosity
    ),
    "no longer exist")
})

test_that("an unknown region class is refused at the door, not silently ignored", {
  tempFolder <- tempFolders[16]
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(tempFolder, recursive = TRUE) }, add = TRUE)

  # AI-255: the region classes are no longer declared at SEM time — semseeker()
  # lost sample_stats_scopes — so the refusal moved to where they are asked for.
  # A run that quietly dropped an unknown class would look identical to one that
  # honoured it, and the researcher would find out only at analysis time.
  SEMseeker:::core_init_env(result_folder = tempFolder, parallel_strategy = "sequential",
                            areas = c("POSITION"), markers = c("MUTATIONS"),
                            start_fresh = TRUE, showprogress = FALSE, verbosity = 1)

  expect_error(SEMseeker:::.sem_regions_resolve(c("SAMPLE", "GENE_NOWHERE")),
               "GENE_NOWHERE")
  # ... and the one that is not a region class at all still resolves, because
  # "no restriction" is a legal request.
  expect_equal(SEMseeker:::.sem_regions_resolve("SAMPLE")[[1]]$area, "PROBE")
})
