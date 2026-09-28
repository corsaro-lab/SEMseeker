## AI-310: the semantics of a request, and of the identity of what it produces.
##
## These are conceptual invariants rather than behaviours of one function. Each
## of them has already been violated at least once, silently, and in every case
## the run still wrote a file that looked complete. That is the shape of defect
## this file exists to catch.

# ---------------------------------------------------------------------------
# 1. The documentation and the accepted vocabulary cannot drift apart
# ---------------------------------------------------------------------------

test_that("every column the documentation names is a column the package accepts", {
  # association_analysis() documented a `marker` column of inference_details for
  # a long time. It was never in the vocabulary, so a request written from the
  # documentation was refused as carrying an unknown column: the documentation
  # described a package that did not exist. Nothing could have noticed: the two
  # lists live in different files and neither reads the other.
  rd <- system.file("../man/association_analysis.Rd", package = "SEMseeker")
  if (!nzchar(rd) || !file.exists(rd)) rd <- "../../man/association_analysis.Rd"
  skip_if_not(file.exists(rd), "association_analysis.Rd not reachable from here")

  lines <- readLines(rd, warn = FALSE)
  # The columns of inference_details are the NESTED \item entries: they sit
  # inside a \describe{} and are indented, while the arguments of the function
  # itself start at column zero.
  nested <- grep("^\\s+\\\\item\\{", lines, value = TRUE)
  documented <- sub("^\\s+\\\\item\\{([^}]+)\\}.*$", "\\1", nested)
  documented <- unique(documented[nzchar(documented)])
  expect_gt(length(documented), 0)

  accepted <- colnames(SEMseeker:::assoc_validate_inference_schema(
    data.frame(independent_variable = "x", stringsAsFactors = FALSE)))

  undocumented_but_accepted <- setdiff(accepted, documented)
  documented_but_refused    <- setdiff(documented, accepted)

  expect_equal(documented_but_refused, character(0),
               info = paste("documented but not accepted:",
                            paste(documented_but_refused, collapse = ", ")))

  # The other direction is a documentation gap, not a lie: these columns work,
  # they are simply not described in the argument block. Pinned so it can only
  # shrink: a new accepted column that nobody documents will fail here.
  known_gap <- c("covariates", "covariates_dummy", "covariates_pca",
                 "collinearity_check", "transformation_x", "filter_p_value",
                 "samples_sql_condition", "areas_sql_condition",
                 "association_results_sql_condition")
  expect_setequal(undocumented_but_accepted, known_gap)
})

test_that("a retired column is named as retired, not guessed at as a typo", {
  # Both of these were real columns. Telling their author to "register it as
  # legal" or suggesting a near-miss would send them the wrong way; the message
  # has to say the column is gone and what replaced it.
  for (col in c("scopes", "depth_analysis")) {
    d <- data.frame(independent_variable = "x", stringsAsFactors = FALSE)
    d[[col]] <- "whatever"
    err <- expect_error(SEMseeker:::assoc_validate_inference_schema(d))
    expect_match(conditionMessage(err), "no longer exist")
    expect_match(conditionMessage(err), col, fixed = TRUE)
  }
})

# ---------------------------------------------------------------------------
# 2. Identity: two different requests cannot name the same artefact
# ---------------------------------------------------------------------------

test_that("distinct coordinate tuples never compose the same artefact key", {
  # io_artefact_key() is the ONLY compositor of the six-coordinate name, and the
  # name is what a pivot file is called. If two distinct tuples composed the
  # same string, one artefact would silently overwrite the other and the second
  # analysis would read the first one's numbers.
  #
  # The risk is not hypothetical: the key is a flat underscore join and two
  # vocabularies carry an underscore inside their values (N_SHORE, S_SHELF), so
  # nothing about the format prevents a collision, only the vocabularies do.
  # That is exactly the kind of guarantee that has to be measured.
  tempFolder <- tempFolders[20]
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(tempFolder, recursive = TRUE) }, add = TRUE)
  SEMseeker:::core_init_env(result_folder = tempFolder, parallel_strategy = "sequential",
                            areas = "ALL", subareas = "ALL", start_fresh = TRUE,
                            showprogress = FALSE, verbosity = 1)
  ssEnv <- SEMseeker:::core_get_session_info()

  # Taken from the package's own registries, so the test tracks the vocabulary
  # instead of carrying a copy of it that can go stale.
  regions <- unique(ssEnv$default$keys_areas_subareas_default[, c("AREA", "SUBAREA")])
  mf      <- unique(ssEnv$default$keys_markers_figures_default[, c("MARKER", "FIGURE")])
  expect_gt(nrow(regions), 1)
  expect_gt(nrow(mf), 1)

  grid <- expand.grid(row_r = seq_len(nrow(regions)), row_m = seq_len(nrow(mf)),
                      SCOPE = SEMseeker:::io_scope_vocabulary(),
                      AGGREGATION = SEMseeker:::util_aggregation_vocabulary(),
                      stringsAsFactors = FALSE)
  tuples <- data.frame(
    MARKER      = mf$MARKER[grid$row_m],
    FIGURE      = mf$FIGURE[grid$row_m],
    SCOPE       = grid$SCOPE,
    AREA        = regions$AREA[grid$row_r],
    SUBAREA     = regions$SUBAREA[grid$row_r],
    AGGREGATION = grid$AGGREGATION,
    stringsAsFactors = FALSE)

  keys <- vapply(seq_len(nrow(tuples)), function(i)
    SEMseeker:::io_artefact_key(tuples$MARKER[i], tuples$FIGURE[i], tuples$SCOPE[i],
                                tuples$AREA[i], tuples$SUBAREA[i], tuples$AGGREGATION[i]),
    character(1))

  dup <- keys[duplicated(keys)]
  expect_equal(unique(dup), character(0),
               info = paste("colliding keys:", paste(unique(dup), collapse = ", ")))
  expect_equal(length(unique(keys)), nrow(tuples))
})

test_that("every coordinate is load-bearing: dropping one changes the key", {
  # A coordinate that never changes the name is a coordinate the name does not
  # carry, and an artefact whose name omits one cannot be told apart from an
  # artefact that omits a different one.
  base <- list(marker = "MUTATIONS", figure = "HYPER", scope = "SAMPLE",
               area = "GENE", subarea = "TSS1500", aggregation = "SUM")
  ref  <- do.call(SEMseeker:::io_artefact_key, base)

  alt <- list(marker = "LESIONS", figure = "HYPO", scope = "INSTANCE",
              area = "ISLAND", subarea = "WHOLE", aggregation = "MEAN")
  for (coord in names(base)) {
    changed <- base
    changed[[coord]] <- alt[[coord]]
    expect_false(identical(ref, do.call(SEMseeker:::io_artefact_key, changed)),
                 info = paste("changing", coord, "left the key unchanged"))
  }
})

# ---------------------------------------------------------------------------
# 3. Identity of the VALUES: two combinations must not answer with one number
# ---------------------------------------------------------------------------

test_that("the taxonomy forbids the one case where aggregations would coincide", {
  # A block holding a single position has no reduction to speak of: its sum, its
  # mean and its median are the block itself. Rather than let three names return
  # one number, which would invite the reader to believe three things were
  # computed: the taxonomy admits exactly one name there, VALUE.
  #
  # So the constraint "two combinations never give the same answer" is kept at
  # the level where it can be kept: the combinations that would collide are not
  # expressible.
  for (area in c("PROBE", "POSITION")) {
    admissible <- SEMseeker:::util_aggregations_allowed(
      "MUTATIONS", "HYPER", discrete = TRUE, default = FALSE,
      scope = "INSTANCE", area = area)
    expect_equal(admissible, "VALUE",
                 info = paste("single-position area", area,
                              "admits more than one aggregation"))
  }

  # At SCOPE = SAMPLE the same area is the whole sample, thousands of positions,
  # so the reduction is real and the names are not interchangeable.
  at_sample <- SEMseeker:::util_aggregations_allowed(
    "MUTATIONS", "HYPER", discrete = TRUE, default = FALSE,
    scope = "SAMPLE", area = "PROBE")
  expect_true(length(at_sample) > 1)
  expect_false("VALUE" %in% at_sample)

  # And the two peaks of a bimodal density exist only where there is a
  # distribution to find them in.
  per_instance <- SEMseeker:::util_aggregations_allowed(
    "SIGNAL", "BETA", discrete = FALSE, default = FALSE, scope = "INSTANCE",
    area = "GENE")
  expect_false(any(c("MODELOW", "MODEHIGH") %in% per_instance))
})

test_that("no two (region class, aggregation) pairs answer with the same numbers", {
  # The constraint stated on real artefacts. Every collapsed artefact of the run
  # is one row of numbers, one per sample; two different questions must not
  # return the same row.
  #
  # This is where the unrestricted-mask defect showed: every region class
  # returned the burden of the whole sample, so the whole matrix below collapsed
  # to a single distinct row while every name in it promised something else.
  tempFolder <- tempFolders[21]
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE) }, add = TRUE)

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
    areas             = c("GENE", "PROBE"),
    subareas          = c("WHOLE", "TSS1500"),
    markers           = c("MUTATIONS"),
    start_fresh       = TRUE,
    inpute            = "median",
    showprogress      = showprogress,
    verbosity         = verbosity)

  classes <- list(c("PROBE", "WHOLE"), c("GENE", "WHOLE"), c("GENE", "TSS1500"))
  aggs    <- SEMseeker:::util_aggregations_allowed("MUTATIONS", "HYPER",
                                                   discrete = TRUE, default = FALSE,
                                                   scope = "SAMPLE", area = "GENE")
  expect_true(all(c("SUM", "MEAN") %in% aggs))

  answers <- list()
  labels  <- character(0)
  for (cl in classes) for (ag in aggs) {
    p <- SEMseeker:::io_read_pivot("MUTATIONS", "HYPER", cl[1], cl[2],
                                   aggregation = ag, scope = "SAMPLE")
    skip_if(is.null(p), paste("artefact unavailable:", cl[1], cl[2], ag))
    df <- as.data.frame(p$collect())
    v  <- unlist(df[1, setdiff(colnames(df), "AREA")], use.names = TRUE)
    answers[[length(answers) + 1L]] <- unname(v[order(names(v))])
    labels <- c(labels, paste(cl[1], cl[2], ag, sep = "_"))
  }
  names(answers) <- labels
  expect_equal(length(answers), length(classes) * length(aggs))

  # every pair, named, so a failure says which two questions gave one answer
  for (i in seq_along(answers)) for (j in seq_along(answers)) if (i < j)
    expect_false(isTRUE(all.equal(answers[[i]], answers[[j]])),
                 info = paste(labels[i], "and", labels[j],
                              "returned the same numbers"))
})

# ---------------------------------------------------------------------------
# 4. What an enrichment can be about
# ---------------------------------------------------------------------------

test_that("the enrichment input invariant is declared once, not five times", {
  # It used to be two literals written out by hand in every backend, so nothing
  # stated the invariant and nothing enforced it: a sixth backend asking for
  # ISLAND would have been served, because assoc_results_get() accepts any
  # region class: the cross-study overlaps iterate over all of them.
  want <- SEMseeker:::enrich_input_invariant()
  expect_equal(want$scope, "INSTANCE")
  expect_equal(want$area,  "GENE")

  backends <- c("enrich_ctdR", "enrich_PathfindR", "enrich_phenotype_phenolyzer",
                "enrich_STRINGdb", "enrich_WebGestalt")
  for (b in backends) {
    f <- system.file(paste0("../R/", b, ".R"), package = "SEMseeker")
    if (!nzchar(f) || !file.exists(f)) f <- file.path("../../R", paste0(b, ".R"))
    skip_if_not(file.exists(f), paste("source of", b, "not reachable"))
    src <- paste(readLines(f, warn = FALSE), collapse = "\n")
    expect_false(grepl('area\\s*=\\s*"GENE"', src),
                 info = paste(b, "still writes the region class out by hand"))
    expect_false(grepl('scope\\s*=\\s*"INSTANCE"', src),
                 info = paste(b, "still writes the scope out by hand"))
    expect_true(grepl("enrich_input_invariant\\(\\)", src),
                info = paste(b, "does not read the declared invariant"))
  }
})

test_that("an enrichment on a folder without gene rows is refused, not answered emptily", {
  # The case AI-308 made easy to reach: a run at scope = "SAMPLE" is a complete,
  # legitimate analysis that cannot feed a pathway enrichment. Before this, every
  # backend read zero rows and wrote nothing, and the folder looked like a study
  # where no pathway was significant.
  tempFolder <- tempFolders[22]
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE) }, add = TRUE)

  # a run whose registry has no GENE class at all: the first of the two refusals
  SEMseeker:::core_init_env(result_folder = tempFolder, parallel_strategy = "sequential",
                            areas = c("POSITION"), markers = c("MUTATIONS"),
                            start_fresh = TRUE, showprogress = FALSE, verbosity = 1)
  details <- data.frame(independent_variable = "Phenotest", family_test = "spearman",
                        scope = "INSTANCE", aggregation = "SUM",
                        stringsAsFactors = FALSE)
  err <- expect_error(SEMseeker:::enrich_input_assert(details))
  expect_match(conditionMessage(err), "a pathway is a set of genes")
  expect_match(conditionMessage(err), "association_analysis", fixed = TRUE)
})
