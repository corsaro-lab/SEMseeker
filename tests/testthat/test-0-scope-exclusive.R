## AI-308: SCOPE picks one branch of aggregation, and that branch restricts.
##
## Two properties, both of which fail *silently*: a result file is produced
## either way, and it looks complete either way.
##
##   1. EXCLUSIVITY. SAMPLE and INSTANCE are alternatives, not addends. They used
##      to be built together and written into one file, so every result was the
##      union of two different questions and no column said which rows answered
##      which. A request now names its branch and gets that branch only.
##
##   2. RESTRICTION. A collapsed burden on GENE_TSS1500 must be summed over the
##      positions of GENE_TSS1500 and no others. It was not: on the Illumina path
##      anno_probe_features_get() hands back the whole annotation table with NA
##      in the column of the class asked for, and the SAMPLE branch selected its
##      coordinates without dropping those rows, so the mask was the entire
##      array whatever the class, and every "restricted" burden came out equal to
##      the burden of the whole sample.
##
## Session tests use tempFolders indices 18 and 19.

.scope_details <- function(scope, family_test = "spearman", aggregation = "SUM") {
  data.frame(
    independent_variable = "Phenotest",
    family_test          = family_test,
    transformation_y     = "",
    transformation_x     = "",
    scope                = scope,
    aggregation          = aggregation,
    filter_p_value       = FALSE,
    stringsAsFactors     = FALSE)
}

# The region classes of the run, and what each of them is expected to cover.
# PROBE_WHOLE is the whole sample (no mask), GENE_WHOLE the probes annotated to
# any gene, GENE_TSS1500 the promoter windows only: a strict chain of subsets,
# which is what makes the restriction testable without trusting a single number.
.scope_areas    <- c("GENE", "PROBE")
.scope_subareas <- c("WHOLE", "TSS1500")

.scope_run_sem <- function(tempFolder) {
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
    areas             = .scope_areas,
    subareas          = .scope_subareas,
    markers           = c("MUTATIONS"),
    start_fresh       = TRUE,
    inpute            = "median",
    showprogress      = showprogress,
    verbosity         = verbosity)
  invisible(syn)
}

.scope_inference_rows <- function(tempFolder) {
  csv_files <- list.files(file.path(tempFolder, "Inference"), pattern = "\\.csv$",
                          recursive = TRUE, full.names = TRUE)
  csv_files <- csv_files[!grepl("(?i)assoc_covariates_model", csv_files)]
  if (length(csv_files) == 0) return(data.frame())
  df <- do.call(plyr::rbind.fill,
                lapply(csv_files, function(f)
                  utils::read.csv2(f, stringsAsFactors = FALSE)))
  # which() and not a bare logical: the CSV carries a job-summary row with NA in
  # SCOPE, and `df[df$SCOPE == x, ]` would materialise it as an all-NA phantom.
  df[which(!is.na(df$SCOPE)), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# the door: what is refused before any result exists
# ---------------------------------------------------------------------------

test_that("a request that does not name its scope is refused", {
  expect_error(SEMseeker:::assoc_validate_scope(.scope_details(NA)),
               "'scope' is required")
  expect_error(SEMseeker:::assoc_validate_scope(.scope_details("")),
               "'scope' is required")
})

test_that("a scope outside the taxonomy is refused, and the row is named", {
  err <- expect_error(SEMseeker:::assoc_validate_scope(.scope_details("AREA")))
  expect_match(conditionMessage(err), "row 1")
  expect_match(conditionMessage(err), "SAMPLE")
  expect_match(conditionMessage(err), "INSTANCE")
})

test_that("a batch family cannot be fitted on a collapsed artefact", {
  # limma/voom estimate a prior variance across the instances they are handed;
  # a SCOPE = SAMPLE artefact holds one row. While the two branches were
  # additive this could be an INFO that skipped the collapsed keys and carried
  # on with the per-instance ones: with one branch per request there is nothing
  # left to carry on with, so it is an error and not an empty file.
  expect_error(
    SEMseeker:::assoc_validate_scope(.scope_details("SAMPLE", family_test = "limma_trend")),
    "scope = \"SAMPLE\"")
  expect_silent(
    SEMseeker:::assoc_validate_scope(.scope_details("INSTANCE", family_test = "limma_trend")))
})

test_that("the scope is returned in canonical spelling", {
  kept <- SEMseeker:::assoc_validate_scope(.scope_details("instance"))
  expect_equal(kept$scope, "INSTANCE")
})

test_that("a request whose model cannot be fitted stops the run", {
  # AI-309: the row used to be dropped and the run continued, so the result file
  # was indistinguishable from one where every requested test had been fitted.
  # A reader cannot tell a model that found nothing from a model never run.
  expect_error(SEMseeker:::assoc_validate_family(.scope_details("SAMPLE", family_test = NA)),
               "'family_test' is required")
  expect_error(SEMseeker:::assoc_validate_family(.scope_details("SAMPLE", family_test = "")),
               "'family_test' is required")
  expect_error(SEMseeker:::assoc_validate_family(.scope_details("SAMPLE", family_test = "not_a_test")),
               "not a model this package can fit")
  # the row number travels with the refusal, so a batch says which request failed
  batch <- rbind(.scope_details("SAMPLE"), .scope_details("SAMPLE", family_test = "nope"))
  expect_error(SEMseeker:::assoc_validate_family(batch), "row 2")
  # and a legal one passes without noise
  expect_silent(SEMseeker:::assoc_validate_family(.scope_details("INSTANCE")))
})

# ---------------------------------------------------------------------------
# the two branches, on a real run
# ---------------------------------------------------------------------------

test_that("scope = SAMPLE tests the collapsed artefacts and leaves the instances alone", {
  tempFolder <- tempFolders[18]
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE) }, add = TRUE)

  .scope_run_sem(tempFolder)
  SEMseeker:::core_close_env()

  # Everything the SEM run left behind, before the analysis adds anything. The
  # position pivot is among these: it is the SOURCE every aggregate is built
  # from, not an aggregate, so it must not be counted against the branch.
  before <- basename(list.files(tempFolder, pattern = "\\.parquet$",
                                recursive = TRUE, full.names = TRUE))

  expect_no_error(
    SEMseeker:::association_analysis(
      inference_details = .scope_details("SAMPLE"),
      result_folder     = tempFolder,
      parallel_strategy = "sequential",
      markers           = c("MUTATIONS"),
      areas             = .scope_areas,
      subareas          = .scope_subareas,
      multiple_test_adj = "BH",
      start_fresh       = TRUE,
      showprogress      = showprogress,
      verbosity         = verbosity))

  rows <- .scope_inference_rows(tempFolder)
  expect_gt(nrow(rows), 0)

  # (1) exclusivity, in the result
  expect_equal(unique(rows$SCOPE), "SAMPLE")

  # the classes of the run, all of them, and the single-position one normalised
  # to the technology's own: the run declared POSITION implicitly and PROBE
  # explicitly, and collapsed they are the same number
  expect_setequal(unique(paste(rows$AREA, rows$SUBAREA, sep = "_")),
                  c("GENE_WHOLE", "GENE_TSS1500", "PROBE_WHOLE"))

  # (2) exclusivity, in what was computed: no per-instance aggregate was built.
  # A file whose key carries INSTANCE and an aggregation other than VALUE is an
  # aggregate; INSTANCE + VALUE is the position pivot, the source.
  added <- setdiff(basename(list.files(tempFolder, pattern = "\\.parquet$",
                                       recursive = TRUE, full.names = TRUE)), before)
  instance_aggregates <- added[grepl("_INSTANCE_", added) & !grepl("_VALUE_", added)]
  expect_equal(instance_aggregates, character(0))
})

test_that("the collapsed burden is restricted to its own region class", {
  tempFolder <- tempFolders[18]
  skip_if(!dir.exists(file.path(tempFolder, "Inference")),
          "scope = SAMPLE run did not complete")
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE) }, add = TRUE)
  SEMseeker:::core_init_env(result_folder = tempFolder, parallel_strategy = "sequential",
                            areas = .scope_areas, subareas = .scope_subareas,
                            markers = c("MUTATIONS"), start_fresh = FALSE,
                            showprogress = showprogress, verbosity = verbosity)

  # The mask itself: how many positions each class contributes. This is the
  # measurement that was wrong, it returned the whole array for every class,
  # and it is deterministic, so it can be asserted strictly.
  mask_size <- function(area, subarea) {
    m <- SEMseeker:::.io_pivot_masked_lazy("MUTATIONS", "HYPER", "SAMPLE", area, subarea)
    if (is.null(m)) return(NA_integer_)
    nrow(as.data.frame(m$lazy$collect()))
  }
  whole    <- mask_size("PROBE", "WHOLE")
  gene     <- mask_size("GENE",  "WHOLE")
  promoter <- mask_size("GENE",  "TSS1500")

  expect_gt(promoter, 0)
  # a strict chain of subsets: promoter windows inside genes inside the sample
  expect_lt(promoter, gene)
  expect_lt(gene, whole)

  # and the burdens that follow from those masks, per sample
  burden <- function(area, subarea) {
    p <- SEMseeker:::io_read_pivot("MUTATIONS", "HYPER", area, subarea,
                                   aggregation = "SUM", scope = "SAMPLE")
    skip_if(is.null(p), paste("collapsed artefact missing for", area, subarea))
    df <- as.data.frame(p$collect())
    v <- unlist(df[1, setdiff(colnames(df), "AREA")], use.names = TRUE)
    v[order(names(v))]
  }
  b_whole    <- burden("PROBE", "WHOLE")
  b_gene     <- burden("GENE",  "WHOLE")
  b_promoter <- burden("GENE",  "TSS1500")

  # containment holds sample by sample: summing a subset of non-negative counts
  # cannot exceed summing the superset
  expect_true(all(b_promoter <= b_gene))
  expect_true(all(b_gene     <= b_whole))
  # and the three are genuinely different aggregates, not one number relabelled
  expect_false(identical(unname(b_promoter), unname(b_gene)))
  expect_false(identical(unname(b_gene),     unname(b_whole)))
})

test_that("scope = INSTANCE tests the per-instance artefacts and nothing else", {
  tempFolder <- tempFolders[19]
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE) }, add = TRUE)

  .scope_run_sem(tempFolder)
  SEMseeker:::core_close_env()

  expect_no_error(
    SEMseeker:::association_analysis(
      inference_details = .scope_details("INSTANCE"),
      result_folder     = tempFolder,
      parallel_strategy = "sequential",
      markers           = c("MUTATIONS"),
      areas             = .scope_areas,
      subareas          = .scope_subareas,
      multiple_test_adj = "BH",
      start_fresh       = TRUE,
      showprogress      = showprogress,
      verbosity         = verbosity))

  rows <- .scope_inference_rows(tempFolder)
  expect_gt(nrow(rows), 0)
  expect_equal(unique(rows$SCOPE), "INSTANCE")
  # one row per instance, so the same class carries many AREA_OF_TEST values,
  # which is exactly what a collapsed row does not have
  gene_rows <- rows[which(rows$AREA == "GENE" & rows$SUBAREA == "WHOLE"), , drop = FALSE]
  expect_gt(length(unique(gene_rows$AREA_OF_TEST)), 1)
})

test_that("a request that wants both branches writes two rows, and gets both", {
  tempFolder <- tempFolders[19]
  skip_if(!dir.exists(file.path(tempFolder, "Data")), "INSTANCE run did not complete")
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE) }, add = TRUE)

  both <- rbind(.scope_details("SAMPLE"), .scope_details("INSTANCE"))

  expect_no_error(
    SEMseeker:::association_analysis(
      inference_details = both,
      result_folder     = tempFolder,
      parallel_strategy = "sequential",
      markers           = c("MUTATIONS"),
      areas             = .scope_areas,
      subareas          = .scope_subareas,
      multiple_test_adj = "BH",
      start_fresh       = TRUE,
      showprogress      = showprogress,
      verbosity         = verbosity))

  rows <- .scope_inference_rows(tempFolder)
  # the file name carries neither the scope nor the aggregation, so both rows
  # land in the same CSV, and the coordinates are what tells them apart
  expect_setequal(unique(rows$SCOPE), c("SAMPLE", "INSTANCE"))
  sample_rows   <- rows[which(rows$SCOPE == "SAMPLE"), , drop = FALSE]
  instance_rows <- rows[which(rows$SCOPE == "INSTANCE"), , drop = FALSE]
  expect_gt(nrow(sample_rows), 0)
  expect_gt(nrow(instance_rows), 0)
  # disjoint: a collapsed row names its class, an instance row names an instance
  expect_equal(length(intersect(sample_rows$AREA_OF_TEST, instance_rows$AREA_OF_TEST)), 0)
})
