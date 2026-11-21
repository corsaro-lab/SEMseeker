# One assembly per run, fixed before anything is written.
#
# Every result is keyed by coordinates, so the assembly is part of the identity
# of each position. An Illumina manifest states its own (hg19 for the 27k, 450k
# and EPIC manifests); a coordinate input states none, so genome_build has no
# default and whoever aligned the reads declares it.
#
# Expected values are the assembly of the manifest as named by the package
# (hg19), the assembly declared in the call, and the absence of the result
# folder after a refusal; none of them is computed by the code under test.

.build_probes <- function() {
  utils::data("test_signal_gse133774", "test_samplesheet_gse133774", "test_master_features",
              package = "SEMseeker", envir = environment())
  demo_probes <- subset(test_master_features, DMR_LABEL %in% c("KCNQ1OT1:TSS-DMR", "H19_IGF2:IG-DMR") &
                          PROBE %in% rownames(test_signal_gse133774))
  list(signal = test_signal_gse133774[demo_probes$PROBE, , drop = FALSE], features = demo_probes,
       sample_sheet = test_samplesheet_gse133774)
}

# The same values as a coordinate data frame, keyed by CHR/START/END.
.build_coordinates <- function(probes) {
  positions <- probes$features[match(rownames(probes$signal), probes$features$PROBE), c("CHR", "START", "END")]
  cbind(positions, as.data.frame(probes$signal), row.names = NULL)
}

.pivot_builds <- function(result_folder) {
  pivots <- list.files(file.path(result_folder, "Data", "Pivots"), "\\.parquet$", recursive = TRUE)
  unique(sub("^.*_(HG[0-9]+|MM[0-9]+|LEGACY)\\.parquet$", "\\1", pivots))
}

test_that("an Illumina matrix takes the assembly of its manifest", {
  skip_on_cran()
  probes <- .build_probes()
  result_folder <- sem_test_folder()
  on.exit({ SEMseeker:::core_close_env(); unlink(result_folder, recursive = TRUE) }, add = TRUE)

  SEMseeker::semseeker(input = probes$signal, sample_sheet = probes$sample_sheet, result_folder = result_folder,
                       parallel_strategy = "sequential", areas = "PROBE")

  expect_identical(readRDS(file.path(result_folder, "Log", "session_info.rds"))$genome_build, "hg19")
  expect_identical(.pivot_builds(result_folder), "HG19")
})

test_that("an Illumina matrix declared on another assembly is refused, and nothing is written", {
  probes <- .build_probes()
  result_folder <- file.path(tempdir(), "build_refused_illumina")
  unlink(result_folder, recursive = TRUE)

  expect_error(SEMseeker::semseeker(input = probes$signal, sample_sheet = probes$sample_sheet, result_folder = result_folder,
                                    genome_build = "hg38"),
               "does not match the Illumina manifest.*hg19")
  expect_false(dir.exists(result_folder))
})

test_that("coordinate input without an assembly is refused, and nothing is written", {
  probes <- .build_probes()
  result_folder <- file.path(tempdir(), "build_refused_coordinates")
  unlink(result_folder, recursive = TRUE)

  expect_error(SEMseeker::semseeker(input = .build_coordinates(probes), sample_sheet = probes$sample_sheet,
                                    result_folder = result_folder),
               "genome_build is required for coordinate input")
  expect_error(SEMseeker::semseeker(input = c("a.bed", "b.bed"), sample_sheet = probes$sample_sheet,
                                    result_folder = result_folder, tech = "LONGREAD"),
               "genome_build is required for coordinate input")
  expect_false(dir.exists(result_folder))
})

test_that("coordinate input runs on the assembly declared", {
  skip_on_cran()
  probes <- .build_probes()
  result_folder <- sem_test_folder()
  on.exit({ SEMseeker:::core_close_env(); unlink(result_folder, recursive = TRUE) }, add = TRUE)

  SEMseeker::semseeker(input = .build_coordinates(probes), sample_sheet = probes$sample_sheet, result_folder = result_folder,
                       genome_build = "hg38", parallel_strategy = "sequential", areas = "PROBE")

  expect_identical(readRDS(file.path(result_folder, "Log", "session_info.rds"))$genome_build, "hg38")
  expect_identical(.pivot_builds(result_folder), "HG38")
})

test_that("an assembly SEMseeker does not know is refused with the values allowed", {
  probes <- .build_probes()
  result_folder <- file.path(tempdir(), "build_refused_unknown")
  unlink(result_folder, recursive = TRUE)

  expect_error(SEMseeker::semseeker(input = .build_coordinates(probes), sample_sheet = probes$sample_sheet,
                                    result_folder = result_folder, genome_build = "hg37"),
               "genome_build must be one of 'hg19', 'hg38', 'mm10', 'legacy'; got 'hg37'")
  expect_false(dir.exists(result_folder))
})

test_that("a result folder keeps the assembly of its first run for good", {
  skip_on_cran()
  probes <- .build_probes()
  result_folder <- sem_test_folder()
  on.exit({ SEMseeker:::core_close_env(); unlink(result_folder, recursive = TRUE) }, add = TRUE)
  run_on <- function(...) SEMseeker::semseeker(input = .build_coordinates(probes), sample_sheet = probes$sample_sheet,
                                               result_folder = result_folder, parallel_strategy = "sequential", areas = "PROBE", ...)
  recorded_build <- function() readRDS(file.path(result_folder, "Log", "session_info.rds"))$genome_build

  run_on(genome_build = "hg19")
  SEMseeker:::core_close_env()
  pivots_before <- list.files(file.path(result_folder, "Data", "Pivots"), recursive = TRUE)

  expect_error(run_on(genome_build = "hg38"), "holds results on hg19 and keeps that assembly for good")
  expect_error(run_on(genome_build = "hg38", start_fresh = TRUE), "holds results on hg19 and keeps that assembly for good")
  expect_identical(recorded_build(), "hg19")
  expect_identical(list.files(file.path(result_folder, "Data", "Pivots"), recursive = TRUE), pivots_before)
  expect_identical(.pivot_builds(result_folder), "HG19")

  # A run that declares no assembly takes the recorded one. Checked on the
  # resolver: rerunning a coordinate run on an existing folder stops in the
  # resume path for a reason unrelated to the assembly, tracked separately.
  expect_identical(SEMseeker:::.sem_genome_build_resolve(.build_coordinates(probes), "auto", NULL, NULL, result_folder), "hg19")
})
