# test-0-analyze_population_bulk.R
#
# Unit test per sem_analyze_population_bulk (AI-042, 2026-06-08): verifica che
# DELTAS, MUTATIONS, DELTAR pivots siano bit-identici alla formula attesa
# (rank-invariante per MUTATIONS, magnitude per DELTAS/DELTAR) usando un
# input sintetico piccolo. LESIONS pivot si verifica solo a livello shape +
# presenza (la logica binomial è esercitata in test-0-lesions_get.R).

test_that("sem_analyze_population_bulk - synthetic 10x5 produces correct pivots", {

  set.seed(42)

  tempFolder <- tempfile("test_apb_")
  dir.create(tempFolder, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tempFolder, recursive = TRUE, force = TRUE), add = TRUE)

  ssEnv <- SEMseeker:::core_init_env(
    tempFolder,
    parallel_strategy = "sequential",
    inpute            = "median",
    bulk_population   = TRUE,
    LESIONS_BP        = 5000L,
    bonferroni_threshold = 0.05
  )

  # ---- Synthetic SIGNAL pivot --------------------------------------------
  n_probes  <- 10L
  n_samples <- 5L
  sample_ids <- paste0("S", seq_len(n_samples))

  signal_matrix <- matrix(runif(n_probes * n_samples), n_probes, n_samples,
                          dimnames = list(NULL, sample_ids))
  probe_df <- data.frame(
    CHR   = c(rep("1", 5L), rep("2", 5L)),
    START = c(seq(100L, 500L, by = 100L), seq(200L, 600L, by = 100L)),
    END   = c(seq(100L, 500L, by = 100L), seq(200L, 600L, by = 100L)) + 1L,
    signal_matrix,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )

  # The FIGURE of SIGNAL is the scale of the session, resolved by the same
  # function sem_analyze_population_bulk() uses to look this pivot up.
  signal_path <- SEMseeker:::io_pivot_file_name_parquet("SIGNAL",
                                                    SEMseeker:::io_signal_figure(),
                                                    "POSITION", "WHOLE")
  dir.create(dirname(signal_path), recursive = TRUE, showWarnings = FALSE)
  polars::as_polars_df(probe_df)$write_parquet(signal_path)

  # ---- Synthetic thresholds (high=0.8, low=0.2 per ogni probe) -----------
  thresholds <- data.frame(
    CHR   = probe_df$CHR,
    START = probe_df$START,
    END   = probe_df$END,
    signal_superior_thresholds = rep(0.8, n_probes),
    signal_inferior_thresholds = rep(0.2, n_probes),
    signal_median_values       = rep(0.5, n_probes),
    stringsAsFactors = FALSE
  )

  sample_sheet <- data.frame(
    Sample_ID    = sample_ids,
    Sample_Group = c("Control", "Control", "Case", "Case", "Case"),
    stringsAsFactors = FALSE
  )

  # ---- Invoke bulk function ----------------------------------------------
  SEMseeker:::sem_analyze_population_bulk(
    signal_data       = NULL,
    sample_sheet      = sample_sheet,
    signal_thresholds = thresholds,
    probe_features    = NULL
  )

  # Helper: align output df rows to input order via (CHR, START, END) key
  # (Polars join non preserva l'ordine delle righe del LHS).
  align_by_coord <- function(actual_df, expected_input_df, sample_ids) {
    key_actual <- paste(actual_df$CHR, actual_df$START, actual_df$END, sep = "|")
    key_input  <- paste(expected_input_df$CHR, expected_input_df$START,
                        expected_input_df$END, sep = "|")
    idx <- match(key_input, key_actual)
    testthat::expect_false(any(is.na(idx)),
                           "Some input rows not found in actual output")
    as.matrix(actual_df[idx, sample_ids])
  }

  # ---- Verify DELTAS_HYPER ------------------------------------------------
  deltas_hyper_path <- SEMseeker:::io_pivot_file_name_parquet("DELTAS", "HYPER",
                                                          "POSITION", "WHOLE")
  testthat::expect_true(file.exists(deltas_hyper_path))

  deltas_hyper_df <- as.data.frame(polars::pl$read_parquet(deltas_hyper_path))
  testthat::expect_equal(nrow(deltas_hyper_df), n_probes)
  testthat::expect_setequal(setdiff(colnames(deltas_hyper_df),
                                    c("CHR", "START", "END")), sample_ids)

  expected_hyper <- pmax(signal_matrix - 0.8, 0)
  actual_hyper   <- align_by_coord(deltas_hyper_df, probe_df, sample_ids)
  testthat::expect_equal(actual_hyper, expected_hyper,
                         tolerance = 1e-10, ignore_attr = TRUE)

  # ---- Verify DELTAS_HYPO ------------------------------------------------
  deltas_hypo_path <- SEMseeker:::io_pivot_file_name_parquet("DELTAS", "HYPO",
                                                         "POSITION", "WHOLE")
  testthat::expect_true(file.exists(deltas_hypo_path))

  deltas_hypo_df <- as.data.frame(polars::pl$read_parquet(deltas_hypo_path))
  expected_hypo  <- pmax(0.2 - signal_matrix, 0)
  actual_hypo    <- align_by_coord(deltas_hypo_df, probe_df, sample_ids)
  testthat::expect_equal(actual_hypo, expected_hypo,
                         tolerance = 1e-10, ignore_attr = TRUE)

  # ---- Verify MUTATIONS_HYPER --------------------------------------------
  mut_hyper_path <- SEMseeker:::io_pivot_file_name_parquet("MUTATIONS", "HYPER",
                                                       "POSITION", "WHOLE")
  testthat::expect_true(file.exists(mut_hyper_path))

  mut_hyper_df <- as.data.frame(polars::pl$read_parquet(mut_hyper_path))
  expected_mut <- (expected_hyper > 0) * 1L
  actual_mut   <- align_by_coord(mut_hyper_df, probe_df, sample_ids)
  testthat::expect_equal(actual_mut, expected_mut, ignore_attr = TRUE)

  # ---- Verify DELTAR_HYPER -----------------------------------------------
  deltar_hyper_path <- SEMseeker:::io_pivot_file_name_parquet("DELTAR", "HYPER",
                                                          "POSITION", "WHOLE")
  testthat::expect_true(file.exists(deltar_hyper_path))

  deltar_hyper_df <- as.data.frame(polars::pl$read_parquet(deltar_hyper_path))
  expected_deltar <- expected_hyper / (0.8 - 0.2)
  actual_deltar   <- align_by_coord(deltar_hyper_df, probe_df, sample_ids)
  testthat::expect_equal(actual_deltar, expected_deltar,
                         tolerance = 1e-9, ignore_attr = TRUE)

  # ---- Verify LESIONS_HYPER existence + shape ----------------------------
  les_hyper_path <- SEMseeker:::io_pivot_file_name_parquet("LESIONS", "HYPER",
                                                       "POSITION", "WHOLE")
  testthat::expect_true(file.exists(les_hyper_path))

  les_hyper_df <- as.data.frame(polars::pl$read_parquet(les_hyper_path))
  testthat::expect_equal(nrow(les_hyper_df), n_probes)
  testthat::expect_setequal(setdiff(colnames(les_hyper_df),
                                    c("CHR", "START", "END")), sample_ids)
  # values must be 0/1 integers
  les_mat <- as.matrix(les_hyper_df[, sample_ids])
  testthat::expect_true(all(les_mat %in% c(0L, 1L)))
})


test_that("sem_analyze_population_bulk fails clean if SIGNAL pivot is missing", {

  tempFolder <- tempfile("test_apb_nopivot_")
  dir.create(tempFolder, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tempFolder, recursive = TRUE, force = TRUE), add = TRUE)

  ssEnv <- SEMseeker:::core_init_env(
    tempFolder,
    parallel_strategy = "sequential",
    inpute            = "median",
    bulk_population   = TRUE,
    LESIONS_BP        = 5000L
  )

  thresholds <- data.frame(
    CHR   = "1", START = 100L, END = 101L,
    signal_superior_thresholds = 0.8,
    signal_inferior_thresholds = 0.2,
    signal_median_values       = 0.5
  )

  sample_sheet <- data.frame(
    Sample_ID    = "S1", Sample_Group = "Control",
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    SEMseeker:::sem_analyze_population_bulk(
      signal_data       = NULL,
      sample_sheet      = sample_sheet,
      signal_thresholds = thresholds,
      probe_features    = NULL
    ),
    regexp = "SIGNAL POSITION pivot mancante"
  )
})
