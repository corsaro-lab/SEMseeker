## AI-248 — taxonomy SCOPE x MARKER x FIGURE x AGGREGATION.
##
## Contracts:
##   1. exactly two compositors — one for column names (io_feature_colname),
##      one for file names (io_pivot_file_name*); nothing else pastes names;
##   2. the aggregation vocabulary is declared in one place and computed in one
##      place, so a column called MEDIAN cannot hold a mean;
##   3. the FIGURE of SIGNAL is the scale, not an aggregation.

# ---------------------------------------------------------------------------
# column names
# ---------------------------------------------------------------------------

test_that("io_feature_colname carries the four axes in order", {
  expect_equal(
    SEMseeker:::io_feature_colname("SAMPLE", "SIGNAL", "BETA", "MEDIAN"),
    "SAMPLE_SIGNAL_BETA_MEDIAN")
  expect_equal(
    SEMseeker:::io_feature_colname("GENE_TSS1500", "MUTATIONS", "HYPER", "SUM"),
    "GENE_TSS1500_MUTATIONS_HYPER_SUM")
  expect_equal(
    SEMseeker:::io_feature_colname("SAMPLE", "SIGNAL", "MVALUE", "IQR"),
    "SAMPLE_SIGNAL_MVALUE_IQR")
})

test_that("io_feature_colname drops the axes that do not apply", {
  # the number of positions is a property of the scope, not an aggregation of a
  # marker: no marker, no figure
  expect_equal(
    SEMseeker:::io_feature_colname("SAMPLE", aggregation = "N_PROBES"),
    "SAMPLE_N_PROBES")
  expect_equal(
    SEMseeker:::io_feature_colname("GENE_TSS1500", aggregation = "N_PROBES"),
    "GENE_TSS1500_N_PROBES")
  # a scope is mandatory: an unscoped column would be ambiguous in the file
  expect_error(SEMseeker:::io_feature_colname(""), "scope")
  expect_error(SEMseeker:::io_feature_colname(NULL), "scope")
})

# ---------------------------------------------------------------------------
# the aggregation vocabulary
# ---------------------------------------------------------------------------

test_that("the two modes are admissible only on the bounded scale", {
  admissible_beta <- SEMseeker:::util_aggregations_allowed("SIGNAL", "BETA",
                                                           default = FALSE)
  admissible_mval <- SEMseeker:::util_aggregations_allowed("SIGNAL", "MVALUE",
                                                           default = FALSE)
  expect_true(all(c("MODELOW", "MODEHIGH") %in% admissible_beta))
  expect_false(any(c("MODELOW", "MODEHIGH") %in% admissible_mval))
  # an unbounded marker has no bimodal split to speak of either
  expect_false(any(c("MODELOW", "MODEHIGH") %in%
                     SEMseeker:::util_aggregations_allowed("MUTATIONS", "HYPER",
                                                           default = FALSE)))
})

test_that("a 0/1 marker admits only the two aggregations that carry information", {
  admissible <- SEMseeker:::util_aggregations_allowed("MUTATIONS", "HYPER",
                                                      discrete = TRUE, default = FALSE)
  # the burden, and the burden per position — the density, which is what makes
  # regions of different size comparable
  expect_setequal(admissible, c("SUM", "MEAN"))
  # median / variance / IQR of a vector of zeros and ones are degenerate
  expect_false(any(c("MEDIAN", "VARIANCE", "IQR") %in% admissible))
})

test_that("the continuous markers are one class: DELTAS and DELTAR match SIGNAL", {
  signal <- SEMseeker:::util_aggregations_allowed("SIGNAL", "MVALUE",
                                                  discrete = FALSE, default = FALSE)
  deltas <- SEMseeker:::util_aggregations_allowed("DELTAS", "HYPER",
                                                  discrete = FALSE, default = FALSE)
  deltar <- SEMseeker:::util_aggregations_allowed("DELTAR", "HYPO",
                                                  discrete = FALSE, default = FALSE)
  expect_setequal(deltas, signal)
  expect_setequal(deltar, signal)
  expect_true(all(c("MEAN", "MEDIAN", "VARIANCE", "IQR") %in% deltas))
})

test_that("what is produced by default is narrower than what is admissible", {
  # a count of epimutations is a burden, and the burden is its sum
  expect_equal(SEMseeker:::util_aggregations_allowed("MUTATIONS", "HYPER"), "SUM")
  # every continuous marker carries the whole descriptor set: adding them does
  # not change the historical mean, it only gives it a suffix
  for (marker in c("SIGNAL", "DELTAS", "DELTAR")) {
    produced <- SEMseeker:::util_aggregations_allowed(marker, "HYPER",
                                                      discrete = FALSE)
    expect_true(all(c("MEAN", "MEDIAN", "VARIANCE", "IQR") %in% produced),
                info = marker)
    expect_false("SUM" %in% produced, info = marker)
  }
})

# ---------------------------------------------------------------------------
# what each aggregation computes
# ---------------------------------------------------------------------------

test_that("util_aggregate_values computes what the name says", {
  set.seed(3L)
  values <- stats::rnorm(500, mean = 2, sd = 3)

  expect_equal(SEMseeker:::util_aggregate_values(values, "SUM"), sum(values))
  expect_equal(SEMseeker:::util_aggregate_values(values, "MEAN"), mean(values))
  expect_equal(SEMseeker:::util_aggregate_values(values, "MEDIAN"),
               stats::median(values))
  expect_equal(SEMseeker:::util_aggregate_values(values, "VARIANCE"),
               stats::var(values))
  expect_equal(SEMseeker:::util_aggregate_values(values, "IQR"), stats::IQR(values))
  # AI-255: VALUE is the identity, and it is only meaningful on a single
  # position — asking it of a set is a request that contradicts itself.
  expect_equal(SEMseeker:::util_aggregate_values(0.42, "VALUE"), 0.42)
  expect_error(SEMseeker:::util_aggregate_values(values, "VALUE"),
               "single position")
  # N_PROBES left the taxonomy: it is a property of the imputation, not an
  # aggregation of a marker.
  expect_error(SEMseeker:::util_aggregate_values(values, "N_PROBES"),
               "unknown aggregation")
})

test_that("util_aggregate_values finds the two peaks of a bimodal sample", {
  set.seed(11L)
  values <- c(stats::rbeta(3000, 2, 40), stats::rbeta(3000, 40, 2))
  low  <- SEMseeker:::util_aggregate_values(values, "MODELOW")
  high <- SEMseeker:::util_aggregate_values(values, "MODEHIGH")
  expect_lt(low, 0.5)
  expect_gt(high, 0.5)
})

test_that("util_aggregate_values degrades instead of inventing a number", {
  expect_true(is.na(SEMseeker:::util_aggregate_values(numeric(0), "MEDIAN")))
  expect_true(is.na(SEMseeker:::util_aggregate_values(numeric(0), "MEAN")))
  # non-finite values are not data: the mean of c(1, NA, Inf, 3) is the mean of
  # the two numbers that are.
  expect_equal(SEMseeker:::util_aggregate_values(c(1, NA, Inf, 3), "MEAN"), 2)
  expect_error(SEMseeker:::util_aggregate_values(1:10, "AVERAGE"), "unknown aggregation")
})

# ---------------------------------------------------------------------------
# the figure of SIGNAL is the scale
# ---------------------------------------------------------------------------

test_that("io_signal_figure reports the scale, not an aggregation", {
  expect_equal(SEMseeker:::io_signal_figure(beta = TRUE), "BETA")
  expect_equal(SEMseeker:::io_signal_figure(beta = FALSE), "MVALUE")
  # a session that has not seen the data yet defaults to the bounded scale
  expect_equal(SEMseeker:::io_signal_figure(beta = NA), "BETA")
})

# ---------------------------------------------------------------------------
# file names
# ---------------------------------------------------------------------------

test_that("the pivot name carries the aggregation, and the scale for SIGNAL", {
  tempFolder <- tempFolders[17]
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(tempFolder, recursive = TRUE) }, add = TRUE)

  SEMseeker:::core_init_env(result_folder = tempFolder, parallel_strategy = "sequential",
                            areas = c("POSITION"), markers = c("MUTATIONS"),
                            start_fresh = TRUE, showprogress = FALSE, verbosity = 1)

  burden <- SEMseeker:::io_pivot_file_name_parquet("MUTATIONS", "HYPER", "CHR",
                                                   "CYTOBAND", aggregation = "SUM")
  expect_match(basename(burden),
               "^MUTATIONS_HYPER_INSTANCE_CHR_CYTOBAND_SUM_HG19\\.parquet$",
               ignore.case = TRUE)

  beta_mean <- SEMseeker:::io_pivot_file_name_parquet("SIGNAL", "BETA", "CHR",
                                                      "CYTOBAND", aggregation = "MEAN")
  mval_mean <- SEMseeker:::io_pivot_file_name_parquet("SIGNAL", "MVALUE", "CHR",
                                                      "CYTOBAND", aggregation = "MEAN")
  # the two scales must not overwrite each other in the same folder
  expect_false(identical(beta_mean, mval_mean))

  # AI-255: omitting the aggregation resolves to VALUE where the block is a
  # single position — the identity now has a name instead of being an absence.
  position <- SEMseeker:::io_pivot_file_name_parquet("MUTATIONS", "HYPER",
                                                     "POSITION", "WHOLE")
  expect_match(basename(position),
               "^MUTATIONS_HYPER_INSTANCE_POSITION_WHOLE_VALUE_HG19\\.parquet$",
               ignore.case = TRUE)

  # ... and stays an error where it is not, because a block holding many
  # positions was reduced somehow and the name has to say how.
  expect_error(
    SEMseeker:::io_pivot_file_name_parquet("MUTATIONS", "HYPER", "GENE", "WHOLE"),
    "aggregation is required")

  # the same region class collapsed to one number per sample is a different
  # artefact, and says so in one token
  collapsed <- SEMseeker:::io_pivot_file_name_parquet(
    "MUTATIONS", "HYPER", "GENE", "WHOLE", aggregation = "SUM", scope = "SAMPLE")
  expect_match(basename(collapsed),
               "^MUTATIONS_HYPER_SAMPLE_GENE_WHOLE_SUM_HG19\\.parquet$",
               ignore.case = TRUE)
})

test_that("the artefact key is injective over the closed vocabularies", {
  tempFolder <- file.path(tempdir(), "test_key_injective")
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(tempFolder, recursive = TRUE) }, add = TRUE)

  SEMseeker:::core_init_env(result_folder = tempFolder, parallel_strategy = "sequential",
                            start_fresh = TRUE, showprogress = FALSE, verbosity = 1)

  # The key is the identity of a row, so two different coordinate tuples must
  # never compose to the same string. It is not obvious that they cannot:
  # core_name_cleaning() maps every separator to "_", and two vocabularies carry
  # an underscore inside their values (N_SHORE, and CHR_CYTOBAND read as three
  # tokens for two coordinates). This is the test that fails loudly the day
  # someone adds a subarea called LOW or a marker called MODE.
  grid <- expand.grid(
    marker      = c("SIGNAL", "MUTATIONS", "DELTAS", "DELTAR", "DELTARP"),
    figure      = c("HYPER", "HYPO", "BETA", "MVALUE"),
    scope       = SEMseeker:::io_scope_vocabulary(),
    area        = c("PROBE", "POSITION", "GENE", "ISLAND", "CHR", "DMR"),
    subarea     = c("WHOLE", "TSS200", "TSS1500", "BODY", "OPENSEA", "CYTOBAND",
                    "N_SHORE", "S_SHORE", "N_SHELF", "S_SHELF"),
    aggregation = SEMseeker:::util_aggregation_vocabulary(),
    stringsAsFactors = FALSE)

  keys <- vapply(seq_len(nrow(grid)), function(i)
    SEMseeker:::io_artefact_key(grid$marker[i], grid$figure[i], grid$scope[i],
                                grid$area[i], grid$subarea[i],
                                grid$aggregation[i]),
    character(1))

  expect_equal(length(unique(keys)), nrow(grid))
})

test_that("every coordinate of the key is required", {
  tempFolder <- file.path(tempdir(), "test_key_required")
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(tempFolder, recursive = TRUE) }, add = TRUE)

  SEMseeker:::core_init_env(result_folder = tempFolder, parallel_strategy = "sequential",
                            start_fresh = TRUE, showprogress = FALSE, verbosity = 1)

  expect_error(
    SEMseeker:::io_artefact_key("MUTATIONS", "HYPER", "INSTANCE", "GENE", "WHOLE", ""),
    "every coordinate is required")
  expect_error(SEMseeker:::io_scope_validate("DEPTH2"), "is not a scope")
})

# ---------------------------------------------------------------------------
# coherence of the request, checked at the door
# ---------------------------------------------------------------------------

.tax_keys <- function(...) {
  k <- data.frame(..., stringsAsFactors = FALSE)
  k
}

test_that("a request that names no aggregation is refused, in either scope", {
  # AI-308: the loop used to run over depth 1, 2, 3: three values of a column
  # that had already been retired, so it ran the same assertion three times. The
  # axis that does vary is the scope, and it varies in two.
  keys <- .tax_keys(MARKER = "MUTATIONS", FIGURE = "HYPER", DISCRETE = TRUE)
  for (scope in c("SAMPLE", "INSTANCE")) {
    details <- data.frame(independent_variable = "Phenotest", family_test = "spearman",
                          aggregation = NA, scope = scope,
                          stringsAsFactors = FALSE)
    expect_error(SEMseeker:::assoc_validate_aggregation(details, keys), "required")
  }
})

test_that("an aggregation outside the taxonomy is refused", {
  keys <- .tax_keys(MARKER = "MUTATIONS", FIGURE = "HYPER", DISCRETE = TRUE)
  details <- data.frame(independent_variable = "Phenotest", family_test = "spearman",
                        aggregation = "AVERAGE", scope = "SAMPLE",
                        stringsAsFactors = FALSE)
  expect_error(SEMseeker:::assoc_validate_aggregation(details, keys),
               "not an aggregation")
})

test_that("an impossible request drops its row instead of stopping the batch", {
  keys <- .tax_keys(MARKER = c("MUTATIONS", "SIGNAL"),
                    FIGURE = c("HYPER", "BETA"),
                    DISCRETE = c(TRUE, FALSE))
  details <- data.frame(
    independent_variable = c("Phenotest", "Phenotest"),
    family_test          = c("spearman", "spearman"),
    # MODELOW exists only for SIGNAL/BETA, so the first row is possible;
    # a run with only counts would make it impossible
    aggregation          = c("MEDIAN", "SUM"),
    scope                = c("SAMPLE", "SAMPLE"),
    stringsAsFactors = FALSE)

  expect_warning(kept <- SEMseeker:::assoc_validate_aggregation(details, keys),
                 "not admissible for MUTATIONS/HYPER")
  # partially admissible: the row survives, the impossible pairs are skipped
  expect_equal(nrow(kept), 2L)

  # now the same median against counts only: nothing can be computed
  counts_only <- .tax_keys(MARKER = "MUTATIONS", FIGURE = "HYPER", DISCRETE = TRUE)
  one_row <- details[1, , drop = FALSE]
  expect_error(
    suppressWarnings(SEMseeker:::assoc_validate_aggregation(one_row, counts_only)),
    "nothing left to analyse")
})

test_that("a batch keeps the rows it can run and drops the ones it cannot", {
  keys <- .tax_keys(MARKER = "MUTATIONS", FIGURE = "HYPER", DISCRETE = TRUE)
  details <- data.frame(
    independent_variable = c("Phenotest", "Phenotest"),
    family_test          = c("spearman", "spearman"),
    aggregation          = c("MEDIAN", "SUM"),
    scope                = c("SAMPLE", "SAMPLE"),
    stringsAsFactors = FALSE)

  expect_warning(kept <- SEMseeker:::assoc_validate_aggregation(details, keys),
                 "row 1")
  expect_equal(nrow(kept), 1L)
  expect_equal(kept$aggregation, "SUM")
})

test_that("a possible request passes silently", {
  keys <- .tax_keys(MARKER = "MUTATIONS", FIGURE = "HYPER", DISCRETE = TRUE)
  details <- data.frame(independent_variable = "Phenotest", family_test = "spearman",
                        aggregation = "SUM", scope = "SAMPLE",
                        stringsAsFactors = FALSE)
  expect_silent(kept <- SEMseeker:::assoc_validate_aggregation(details, keys))
  expect_equal(nrow(kept), 1L)
})

# ---------------------------------------------------------------------------
# the numbers, not just the names
#
# Everything above checks that the taxonomy is COHERENT: the right names, the
# right admissibility, the right refusals. None of it checks that a number
# coming out of the pipeline is the number it claims to be — and that is the gap
# the original defect lived in. `MEDIAN` on GENE_TSS1500 returned the mean, and
# every structural check passed while it did, because the mean is a perfectly
# well-formed number under a column called MEDIAN.
#
# So: build a position pivot whose aggregates are known by hand, ask for two
# different aggregations, and compare against an independent calculation.
# ---------------------------------------------------------------------------

test_that("two aggregations of the same artefact are two different numbers", {
  tempFolder <- file.path(tempdir(), "test_two_aggregations")
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(tempFolder, recursive = TRUE) }, add = TRUE)

  SEMseeker:::core_init_env(result_folder = tempFolder, parallel_strategy = "sequential",
                            start_fresh = TRUE, showprogress = FALSE, verbosity = 1)

  # Chosen so mean and median differ, and so neither coincides with the other
  # sample's: S1 mean 0.45 / median 0.40, S2 mean 0.50 / median 0.50.
  s1 <- c(0.10, 0.90, 0.50, 0.30)
  s2 <- c(0.20, 0.80, 0.40, 0.60)
  df <- data.frame(CHR   = c("1", "1", "1", "2"),
                   START = c(100L, 200L, 300L, 400L),
                   END   = c(101L, 201L, 301L, 401L),
                   S1 = s1, S2 = s2)

  base <- SEMseeker:::io_pivot_file_name_parquet("SIGNAL", "BETA", "POSITION", "WHOLE")
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)
  polars::as_polars_df(df)$write_parquet(base)

  written <- SEMseeker:::io_pivot_build("SIGNAL", "BETA", scope = "SAMPLE",
                                        area = "PROBE", subarea = "WHOLE",
                                        aggregations = c("MEAN", "MEDIAN"))
  expect_length(written, 2L)

  read_one <- function(agg) {
    p <- SEMseeker:::io_pivot_file_name_parquet("SIGNAL", "BETA", "PROBE", "WHOLE",
                                                aggregation = agg, scope = "SAMPLE")
    expect_true(file.exists(p))
    as.data.frame(polars::pl$read_parquet(p))
  }

  got_mean   <- read_one("MEAN")
  got_median <- read_one("MEDIAN")

  # A collapsed artefact is one row tall.
  expect_equal(nrow(got_mean), 1L)
  expect_equal(nrow(got_median), 1L)

  # Each value equals the aggregation computed independently, in plain R.
  expect_equal(got_mean$S1,   mean(s1))
  expect_equal(got_mean$S2,   mean(s2))
  expect_equal(got_median$S1, stats::median(s1))
  expect_equal(got_median$S2, stats::median(s2))

  # And the two are NOT the same number: this is the assertion that fails if an
  # aggregation is ever silently answered with another one.
  expect_false(isTRUE(all.equal(got_mean$S1, got_median$S1)))

  # AREA_OF_TEST names the aggregate, since the region class is already said by
  # AREA and SUBAREA.
  expect_equal(got_mean$AREA,   "MEAN")
  expect_equal(got_median$AREA, "MEDIAN")
})

test_that("the FDR family is the key, and the numbers say so", {
  # AI-257. The family over which an FDR is controlled has to be a statistical
  # choice, and twice it had been something else: adjusted per CHUNK inside
  # assoc_apply_stat_model() (so the family was the memory split), then adjusted
  # across the whole file (so it grew with how many aggregations were asked for).
  # It is now the key minus the instance. This test compares against p.adjust()
  # computed independently on each family, which is the only way to tell the
  # three apart: they all produce well-formed p-values.
  p_mean   <- c(0.001, 0.02, 0.30, 0.60)
  p_median <- c(0.005, 0.05)

  results <- data.frame(
    MARKER      = "DELTAS",
    FIGURE      = "HYPO",
    SCOPE       = "INSTANCE",
    AREA        = "GENE",
    SUBAREA     = "WHOLE",
    AGGREGATION = c(rep("MEAN", length(p_mean)), rep("MEDIAN", length(p_median))),
    AREA_OF_TEST = c(paste0("G", seq_along(p_mean)), paste0("G", seq_along(p_median))),
    PVALUE      = c(p_mean, p_median),
    stringsAsFactors = FALSE)

  out <- SEMseeker:::.assoc_adjust_levels(results, method = "BH",
                                          method_label = "BH", alpha = 0.05)

  expect_equal(out$PVALUE_ADJ_KEY_BH[out$AGGREGATION == "MEAN"],
               stats::p.adjust(p_mean, method = "BH"))
  expect_equal(out$PVALUE_ADJ_KEY_BH[out$AGGREGATION == "MEDIAN"],
               stats::p.adjust(p_median, method = "BH"))

  # The point of the family: asking for a second aggregation must not change the
  # correction of the first. Adjusting across the file would have made every
  # adjusted p of MEAN depend on how many MEDIAN rows happened to be there.
  alone <- results[results$AGGREGATION == "MEAN", ]
  expect_equal(SEMseeker:::.assoc_adjust_levels(alone, method = "BH",
                                                method_label = "BH",
                                                alpha = 0.05)$PVALUE_ADJ_KEY_BH,
               out$PVALUE_ADJ_KEY_BH[out$AGGREGATION == "MEAN"])

  # And it is NOT the global adjustment, which is a different column on purpose.
  expect_false(isTRUE(all.equal(
    out$PVALUE_ADJ_KEY_BH, stats::p.adjust(c(p_mean, p_median), method = "BH"))))
})

test_that("the SCOPE level pools what the key separates, and says which is which", {
  # AI-257. The middle family. Every row here shares SCOPE = INSTANCE, so the
  # SCOPE level adjusts across both aggregations at once while the KEY level
  # keeps them apart. Two columns, two families, one set of p-values: if the two
  # ever came out equal the middle level would be buying nothing.
  p_mean   <- c(0.001, 0.02, 0.30, 0.60)
  p_median <- c(0.005, 0.05)

  results <- data.frame(
    MARKER      = "DELTAS",
    FIGURE      = "HYPO",
    SCOPE       = "INSTANCE",
    AREA        = "GENE",
    SUBAREA     = "WHOLE",
    AGGREGATION = c(rep("MEAN", length(p_mean)), rep("MEDIAN", length(p_median))),
    AREA_OF_TEST = c(paste0("G", seq_along(p_mean)), paste0("H", seq_along(p_median))),
    PVALUE      = c(p_mean, p_median),
    stringsAsFactors = FALSE)

  out <- SEMseeker:::.assoc_adjust_levels(results, method = "BH",
                                          method_label = "BH", alpha = 0.05)

  expect_equal(out$PVALUE_ADJ_SCOPE_BH,
               stats::p.adjust(c(p_mean, p_median), method = "BH"))
  expect_false(isTRUE(all.equal(out$PVALUE_ADJ_SCOPE_BH, out$PVALUE_ADJ_KEY_BH)))
})

test_that("at SCOPE = SAMPLE the key holds one row, and the scope level is what answers", {
  # AI-257. The reason there are three levels and not two. A collapsed artefact
  # is one number per sample, so its key holds a single row and BH on n = 1 is
  # the identity: PVALUE_ADJ_KEY_BH == PVALUE, exactly, and a column named for an
  # adjustment that did not happen is how multiplicity goes unreported. The scope
  # level is the one with members to count.
  p <- c(0.001, 0.02, 0.30, 0.60)

  results <- data.frame(
    MARKER      = "MUTATIONS",
    FIGURE      = c("HYPER", "HYPO", "HYPER", "HYPO"),
    SCOPE       = "SAMPLE",
    AREA        = "GENE",
    SUBAREA     = c("TSS1500", "TSS1500", "BODY", "BODY"),
    AGGREGATION = "SUM",
    AREA_OF_TEST = "GENE",
    PVALUE      = p,
    stringsAsFactors = FALSE)

  out <- SEMseeker:::.assoc_adjust_levels(results, method = "BH",
                                          method_label = "BH", alpha = 0.05)

  # Every key is degenerate: four distinct (FIGURE, SUBAREA) pairs, one row each.
  expect_equal(out$PVALUE_ADJ_KEY_BH, p)
  # The scope level counts all four.
  expect_equal(out$PVALUE_ADJ_SCOPE_BH, stats::p.adjust(p, method = "BH"))
  expect_false(isTRUE(all.equal(out$PVALUE_ADJ_SCOPE_BH, out$PVALUE_ADJ_KEY_BH)))
})

test_that("the level a flag answers for is named, not matched by the method string", {
  # AI-257. `grepl(method, colnames)` caught every adjusted column at once, so
  # adding a second level silently turned the significance flag into an AND
  # across levels. The selector names the level.
  results <- data.frame(
    PVALUE                = c(0.001, 0.20),
    PVALUE_ADJ_KEY_BH     = c(0.900, 0.001),
    PVALUE_ADJ_SCOPE_BH   = c(0.900, 0.001),
    PVALUE_ADJ_ALL_BH     = c(0.001, 0.900),
    I_AGE_PVALUE_ADJ_ALL_BH = c(0.001, 0.900),
    stringsAsFactors = FALSE)

  all_cols <- SEMseeker:::.assoc_level_columns(results, "ALL", "BH")
  expect_setequal(all_cols, c("PVALUE_ADJ_ALL_BH", "I_AGE_PVALUE_ADJ_ALL_BH"))
  expect_equal(SEMseeker:::.assoc_level_columns(results, "KEY", "BH"),
               "PVALUE_ADJ_KEY_BH")

  # The flag follows the ALL level alone: row 1 is significant there and not at
  # the two narrow levels, and it must still read TRUE.
  expect_equal(SEMseeker:::.assoc_all_below(results, all_cols, 0.05),
               c(TRUE, FALSE))

  # No column to answer with is not an answer of TRUE.
  expect_true(all(is.na(SEMseeker:::.assoc_all_below(results, "NOSUCH", 0.05))))
})

test_that("one estimator serves the three levels, and an unnameable one yields NA", {
  # AI-257. The narrow level used to have method = "BH" written into it while
  # the global level honoured the run's setting, so a run asking for something
  # else got a column named for one estimator and computed with another.
  p <- c(0.001, 0.02, 0.30, 0.60)

  expect_equal(SEMseeker:::.assoc_adjust_vector(p, method = "bonferroni",
                                                alpha = 0.05),
               stats::p.adjust(p, method = "bonferroni"))
  expect_equal(SEMseeker:::.assoc_adjust_vector(p, method = "BH", alpha = 0.05),
               stats::p.adjust(p, method = "BH"))

  # qvalue cannot estimate pi0 on a family this small. NA is the honest answer:
  # substituting another estimator would leave the column name wrong.
  skip_if_not_installed("qvalue")
  # suppressWarnings: qvalue's own internals partial-match an argument, and
  # testthat turns that into a warning of ours. Not our call site to fix.
  expect_true(all(is.na(suppressWarnings(
    SEMseeker:::.assoc_adjust_vector(c(0.01, 0.4), method = "q", alpha = 0.05)))))
})

test_that("VALUE carries the per-position values, unreduced", {
  tempFolder <- file.path(tempdir(), "test_value_endtoend")
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(tempFolder, recursive = TRUE) }, add = TRUE)

  SEMseeker:::core_init_env(result_folder = tempFolder, parallel_strategy = "sequential",
                            start_fresh = TRUE, showprogress = FALSE, verbosity = 1)

  s1 <- c(0.10, 0.90, 0.50)
  s2 <- c(0.20, 0.80, 0.40)
  df <- data.frame(CHR   = c("1", "1", "2"),
                   START = c(100L, 200L, 300L),
                   END   = c(101L, 201L, 301L),
                   S1 = s1, S2 = s2)

  p <- SEMseeker:::io_pivot_file_name_parquet("SIGNAL", "BETA", "POSITION", "WHOLE")
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  polars::as_polars_df(df)$write_parquet(p)

  # The name resolves to VALUE on its own: at a single position there is nothing
  # to reduce, so the identity is the only meaningful operator.
  expect_match(basename(p), "_VALUE_", fixed = TRUE)

  masked <- SEMseeker:::.io_pivot_masked_lazy("SIGNAL", "BETA", "INSTANCE",
                                              "POSITION", "WHOLE")
  got <- as.data.frame(masked$lazy$collect())

  # One row per position, nothing aggregated, and the three coordinates composed
  # into ONE identifier - the form io_probe_id_to_coord() splits back.
  expect_equal(nrow(got), length(s1))
  expect_equal(got$AREA, paste(df$CHR, df$START, sep = "_"))
  expect_equal(got$S1, s1)
  expect_equal(got$S2, s2)
  # CHR/START/END must not survive as if they were samples.
  expect_false(any(c("CHR", "START", "END") %in% colnames(got)))
})
