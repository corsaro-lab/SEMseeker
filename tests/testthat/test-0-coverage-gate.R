## AI-074 — coverage is a mandatory pre-step of every SEM analysis.
##
## Three contracts, in decreasing order of importance:
##   1. the coverage charts are produced on EVERY run, including the runs the
##      gate rejects — they are the artefact the analyst reads;
##   2. the run STOPS when the input barely overlaps the reference annotation,
##      with a message naming the numbers;
##   3. coverage_minimum lowers the bar explicitly for deliberate
##      cross-technology runs.
##

.coverage_chart_files <- function(tf) {
  list.files(file.path(tf, "Chart", "COVERAGE"), full.names = TRUE)
}

.coverage_sidecar <- function(tf) {
  file.path(tf, "Data", "COVERAGE_GATE.json")
}

test_that("coverage gate passes and writes its sidecar when the input matches the annotation", {
  tf <- sem_test_folder()
  unlink(tf, recursive = TRUE)
  SEMseeker:::core_init_env(result_folder = tf, tech = "K850",
                            parallel_strategy = "sequential", start_fresh = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE); unlink(tf, recursive = TRUE) },
          add = TRUE)

  # real EPIC probe ids from the bundled fixture
  observed <- probe_features$PROBE[1:500]
  metrics <- SEMseeker:::sem_coverage_gate(observed)

  expect_equal(metrics$gate, "pass")
  expect_equal(metrics$input_positions, length(unique(observed)))
  expect_gte(metrics$coverage_percent, 80)
  expect_equal(metrics$coverage_minimum, 80)

  expect_true(file.exists(.coverage_sidecar(tf)))
  back <- jsonlite::fromJSON(.coverage_sidecar(tf))
  expect_equal(back$gate, "pass")
})

test_that("coverage gate stops the run below threshold AND still leaves the charts behind", {
  tf <- sem_test_folder()
  unlink(tf, recursive = TRUE)
  SEMseeker:::core_init_env(result_folder = tf, tech = "K850",
                            parallel_strategy = "sequential", start_fresh = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE); unlink(tf, recursive = TRUE) },
          add = TRUE)

  # identifiers that exist in no manifest: 0% coverage
  observed <- paste0("not_a_probe_", seq_len(300))

  expect_error(SEMseeker:::sem_coverage_gate(observed), "Coverage gate failed")

  # contract 1: charts are the artefact of the run, produced even when the
  # gate rejects it — this is what the analyst opens to see WHAT was missing.
  expect_gt(length(.coverage_chart_files(tf)), 0)

  # the sidecar records the rejection for later audit
  expect_true(file.exists(.coverage_sidecar(tf)))
  back <- jsonlite::fromJSON(.coverage_sidecar(tf))
  expect_equal(back$gate, "fail")
  expect_equal(back$covered_positions, 0)
})

test_that("coverage_minimum lowers the bar explicitly for cross-technology runs", {
  tf <- sem_test_folder()
  unlink(tf, recursive = TRUE)
  SEMseeker:::core_init_env(result_folder = tf, tech = "K850",
                            parallel_strategy = "sequential", start_fresh = TRUE,
                            coverage_minimum = 0)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE); unlink(tf, recursive = TRUE) },
          add = TRUE)

  observed <- paste0("not_a_probe_", seq_len(300))
  metrics <- SEMseeker:::sem_coverage_gate(observed)

  expect_equal(metrics$gate, "pass")
  expect_equal(metrics$coverage_minimum, 0)
})

test_that("coverage gate is not enforced on coordinate-based technologies", {
  tf <- sem_test_folder()
  unlink(tf, recursive = TRUE)
  SEMseeker:::core_init_env(result_folder = tf, tech = "LONGREAD",
                            genome_build = "hg38",
                            parallel_strategy = "sequential", start_fresh = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE); unlink(tf, recursive = TRUE) },
          add = TRUE)

  observed <- c("1:1000-1001", "1:2000-2001", "2:3000-3001")
  metrics <- SEMseeker:::sem_coverage_gate(observed, tech = "LONGREAD")

  expect_equal(metrics$gate, "skipped")
  expect_equal(metrics$coverage_percent, 100)
  expect_true(file.exists(.coverage_sidecar(tf)))
})
