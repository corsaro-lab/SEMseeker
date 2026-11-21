# The DMR area and the CYTOBAND subarea of CHR are not built any more.
#
# Both came from data bundled in the package and fixed on hg19 whatever the
# assembly of the run. They are removed, and asking for them is refused at the
# door of semseeker(), before the result folder exists, like any other region
# class SEMseeker does not build.
#
# Expected values are the names of the removed classes and datasets, written
# here by hand, and the absence of the result folder after a refusal.

.removed_probes <- function() {
  utils::data("test_signal_gse133774", "test_samplesheet_gse133774", "test_master_features",
              package = "SEMseeker", envir = environment())
  demo_probes <- subset(test_master_features, DMR_LABEL %in% c("KCNQ1OT1:TSS-DMR", "H19_IGF2:IG-DMR") &
                          PROBE %in% rownames(test_signal_gse133774))
  list(signal = test_signal_gse133774[demo_probes$PROBE, , drop = FALSE], sample_sheet = test_samplesheet_gse133774)
}

test_that("asking for the DMR area is refused before anything is written", {
  probes <- .removed_probes()
  result_folder <- file.path(tempdir(), "regions_refused_dmr")
  unlink(result_folder, recursive = TRUE)

  expect_error(SEMseeker::semseeker(input = probes$signal, sample_sheet = probes$sample_sheet, result_folder = result_folder,
                                    areas = c("GENE", "DMR")),
               "areas: SEMseeker builds no area 'DMR'. Allowed: GENE, ISLAND, CHR, PROBE, or 'ALL'")
  expect_false(dir.exists(result_folder))
})

test_that("asking for the CYTOBAND subarea is refused before anything is written", {
  probes <- .removed_probes()
  result_folder <- file.path(tempdir(), "regions_refused_cytoband")
  unlink(result_folder, recursive = TRUE)

  expect_error(SEMseeker::semseeker(input = probes$signal, sample_sheet = probes$sample_sheet, result_folder = result_folder,
                                    areas = "CHR", subareas = "CYTOBAND"),
               "subareas: SEMseeker builds no subarea 'CYTOBAND'")
  expect_false(dir.exists(result_folder))
})

test_that("a run on every area writes no DMR or CYTOBAND result, and the annotation declares no such column", {
  skip_on_cran()
  probes <- .removed_probes()
  result_folder <- sem_test_folder()
  on.exit({ SEMseeker:::core_close_env(); unlink(result_folder, recursive = TRUE) }, add = TRUE)

  SEMseeker::semseeker(input = probes$signal, sample_sheet = probes$sample_sheet, result_folder = result_folder,
                       parallel_strategy = "sequential", areas = "ALL", subareas = "WHOLE")

  written <- list.files(result_folder, recursive = TRUE)
  expect_gt(length(grep("_GENE_", written)), 0L)
  expect_length(grep("DMR|CYTOBAND", written), 0L)
  # The schema is what the builder keeps and what its cache is checked against,
  # so a column it declares is a column every annotation carries.
  expect_length(grep("DMR|CYTOBAND", SEMseeker:::.anno_probe_schema("K850")), 0L)
})

test_that("the package ships neither dmr_annotation nor cytoband_hg19", {
  shipped <- utils::data(package = "SEMseeker")$results[, "Item"]
  expect_false("dmr_annotation" %in% shipped)
  expect_false("cytoband_hg19" %in% shipped)
})
