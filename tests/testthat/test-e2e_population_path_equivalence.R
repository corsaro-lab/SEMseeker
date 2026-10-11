# =============================================================================
# The two population paths must produce the same numbers.
#
# semseeker() has two engines for the per-sample markers: the bulk path
# (bulk_population = TRUE, the default) computes every sample at once and
# writes the POSITION pivots directly; the per-sample path
# (bulk_population = FALSE) loops sample by sample, writes one BED/bedgraph
# per sample and marker, and the pivots are merged from those files.
#
# Nothing compared them until this test. It runs the same real fixture through
# both and compares every parquet the two runs write, file by file and cell by
# cell.
#
# The two paths do NOT share a representation, and the comparison is built
# around that on purpose:
#   - the bulk pivots are dense: every position, every sample, 0 where nothing
#     was found;
#   - the per-sample pivots are sparse: a position appears only if some sample
#     has a value there, and a sample with no value at all (zero lesions, for
#     instance) has no column.
# So a missing row or column is read as 0, and only then are the values
# compared. What must agree is the number, not the layout.
# =============================================================================

.equivalence_markers <- c("SIGNAL", "MUTATIONS", "LESIONS", "DELTAS", "DELTAR",
                          "DELTAQ", "DELTARQ", "DELTAP", "DELTARP")

.run_population_path <- function(result_folder, bulk) {
  SEMseeker::semseeker(
    input             = signal_data,
    sample_sheet      = mySampleSheet,
    result_folder     = result_folder,
    parallel_strategy = "sequential",
    areas             = "POSITION",
    markers           = .equivalence_markers,
    start_fresh       = TRUE,
    inpute            = "median",
    bulk_population   = bulk,
    showprogress      = showprogress,
    verbosity         = verbosity
  )
  file.path(result_folder, "Data")
}

# The parquet files a run wrote, relative to its Data folder.
.parquet_files <- function(data_folder) {
  sort(list.files(data_folder, pattern = "\\.parquet$", recursive = TRUE))
}

# One parquet as a value matrix keyed by its non-value columns. Value columns
# are the numeric ones other than START/END; everything else (CHR, START, END,
# AREA, ...) identifies the row.
.as_keyed_matrix <- function(path) {
  df <- as.data.frame(polars::pl$read_parquet(path))
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  value_cols   <- setdiff(numeric_cols, c("START", "END"))
  key_cols     <- setdiff(names(df), value_cols)
  list(key    = do.call(paste, c(df[key_cols], sep = "|")),
       values = as.matrix(df[value_cols]))
}

# Both files on the union of their rows and value columns, absent read as 0.
.aligned <- function(path_a, path_b) {
  a <- .as_keyed_matrix(path_a)
  b <- .as_keyed_matrix(path_b)
  keys <- union(a$key, b$key)
  cols <- union(colnames(a$values), colnames(b$values))
  fill <- function(m) {
    out <- matrix(0, length(keys), length(cols), dimnames = list(keys, cols))
    rows <- match(rownames(out), m$key)
    for (cl in intersect(cols, colnames(m$values))) {
      out[!is.na(rows), cl] <- m$values[rows[!is.na(rows)], cl]
    }
    out[is.na(out)] <- 0
    out
  }
  list(a = fill(a), b = fill(b))
}

test_that("the bulk and the per-sample population paths write the same values", {
  folder_bulk   <- sem_test_folder()
  folder_legacy <- sem_test_folder()
  on.exit({
    try(SEMseeker:::core_close_env(), silent = TRUE)
    unlink(c(folder_bulk, folder_legacy), recursive = TRUE)
  }, add = TRUE)

  data_bulk   <- .run_population_path(folder_bulk,   bulk = TRUE)
  data_legacy <- .run_population_path(folder_legacy, bulk = FALSE)

  files_bulk   <- .parquet_files(data_bulk)
  files_legacy <- .parquet_files(data_legacy)

  # Same artefacts: a pivot one path writes and the other does not is already
  # a divergence, whatever its values.
  expect_identical(files_bulk, files_legacy)

  # Every marker asked for is present, so the comparison below cannot pass by
  # comparing nothing.
  pivot_markers <- unique(sub("_.*$", "", basename(dirname(
    files_bulk[grepl("^Pivots/", files_bulk)]))))
  expect_setequal(pivot_markers, .equivalence_markers)

  for (f in intersect(files_bulk, files_legacy)) {
    both <- .aligned(file.path(data_bulk, f), file.path(data_legacy, f))
    expect_equal(both$a, both$b, tolerance = 1e-12, info = f)
  }

  # The one thing the per-sample path gives that the bulk path does not: a
  # BED/bedgraph per sample and marker. It is why the bulk path needs an
  # exporter before the per-sample path can go.
  per_sample <- function(d) list.files(d, pattern = "\\.(bed|bedgraph)(\\.gz)?$",
                                       recursive = TRUE)
  expect_length(per_sample(data_bulk), 0L)
  expect_gt(length(per_sample(data_legacy)), 0L)
})
