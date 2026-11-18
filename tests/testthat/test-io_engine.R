# Foundations of the DuckDB engine: the connection of the process, its limits,
# and the three primitives the rest of the package will use in place of polars
# (scan a parquet, register a data frame, sink a query to parquet).
#
# The limits are checked by reading them back from DuckDB, never by trusting
# what was passed: memory_limit and temp_directory given to dbConnect(config = )
# are ignored without a warning, and the benchmark that found it had been
# measuring an unlimited engine while believing it limited.

#' Value of one DuckDB setting as DuckDB reports it.
#'
#' @param connection an open DuckDB connection
#' @param setting_name the setting, for instance "memory_limit" or "threads"
#' @return the setting as text, for instance "1.0 GiB"
engine_setting <- function(connection, setting_name) {
  DBI::dbGetQuery(connection, sprintf("SELECT current_setting('%s') AS setting_value", setting_name))$setting_value
}

#' Bytes of a memory amount printed by DuckDB.
#'
#' @param amount_text text such as "32.0 GiB" or "500.0 MiB"
#' @return the amount in bytes, as a double
amount_bytes <- function(amount_text) {
  parts <- strsplit(trimws(amount_text), " ")[[1]]
  unit_bytes <- c(bytes = 1, KiB = 1024, MiB = 1024^2, GiB = 1024^3, TiB = 1024^4)
  as.numeric(parts[1]) * unit_bytes[[parts[2]]]
}

#' Open a session for these tests, with the engine options under test.
#'
#' @param ... session options passed to core_init_env(), memory_limit among them
#' @return the session list returned by core_init_env()
engine_session <- function(...) {
  SEMseeker:::core_init_env(sem_test_folder(), parallel_strategy = "sequential", start_fresh = TRUE, ...)
}

test_that("the memory limit and the threads asked for are the ones in force", {
  ssEnv <- engine_session(memory_limit = "1GiB")
  on.exit(SEMseeker:::core_close_env(), add = TRUE)
  connection <- SEMseeker:::io_engine_connection()

  expect_identical(engine_setting(connection, "memory_limit"), "1.0 GiB")
  expect_identical(as.integer(engine_setting(connection, "threads")), as.integer(ssEnv$parallel$nCore))
})

test_that("the default limit is half of the physical memory, read by DuckDB on every system", {
  engine_session()
  on.exit(SEMseeker:::core_close_env(), add = TRUE)
  connection <- SEMseeker:::io_engine_connection()

  # DuckDB's own default is 80% of the physical memory it detects, so a plain
  # connection gives the denominator without any code specific to one system.
  plain_connection <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(plain_connection, shutdown = TRUE), add = TRUE)
  ratio <- amount_bytes(engine_setting(connection, "memory_limit")) / amount_bytes(engine_setting(plain_connection, "memory_limit"))

  expect_equal(ratio, 50 / 80, tolerance = 0.01)
})

test_that("a limit passed to dbConnect(config = ) is not applied, so the check above would catch it", {
  # Seeded defect: the way of opening the connection that the benchmark first
  # used. If DuckDB honoured config here this test would fail, and the reading
  # back in io_engine_connection() could be simplified.
  wrong_connection <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE), config = list(memory_limit = "1GiB"))
  on.exit(DBI::dbDisconnect(wrong_connection, shutdown = TRUE), add = TRUE)

  expect_false(identical(engine_setting(wrong_connection, "memory_limit"), "1.0 GiB"))
})

test_that("the spill directory lives inside the session temporary folder", {
  ssEnv <- engine_session()
  on.exit(SEMseeker:::core_close_env(), add = TRUE)
  connection <- SEMseeker:::io_engine_connection()

  spill_directory <- normalizePath(engine_setting(connection, "temp_directory"), mustWork = FALSE)
  session_temp <- normalizePath(ssEnv$temp_folder, mustWork = FALSE)
  expect_true(startsWith(spill_directory, session_temp))
})

test_that("a data frame written to parquet and read back by a third reader is unchanged", {
  skip_if_not_installed("nanoparquet")
  engine_session()
  on.exit(SEMseeker:::core_close_env(), add = TRUE)

  original <- data.frame(
    CHR    = c("1", "2", "X", NA),
    START  = c(10L, 20L, NA, 40L),
    VALUE  = c(0.123456789012345, NA, -1e-12, 3),
    FLAG   = c(TRUE, FALSE, NA, TRUE),
    stringsAsFactors = FALSE
  )
  parquet_file <- tempfile(fileext = ".parquet")
  on.exit(unlink(parquet_file), add = TRUE)

  SEMseeker:::io_sink_parquet(SEMseeker:::io_register_frame(original), parquet_file)
  # nanoparquet is neither the engine that wrote the file nor polars, which the
  # migration removes: the comparison does not depend on the code under test.
  read_back <- as.data.frame(nanoparquet::read_parquet(parquet_file))

  expect_identical(read_back, original)
  expect_identical(unique(nanoparquet::read_parquet_metadata(parquet_file)$column_chunks$codec), "ZSTD")
})

test_that("a parquet file is scanned lazily and collected unchanged", {
  skip_if_not_installed("nanoparquet")
  engine_session()
  on.exit(SEMseeker:::core_close_env(), add = TRUE)

  original <- data.frame(CHR = c("1", "1", "2"), START = c(5L, 9L, 1L), VALUE = c(0.5, NA, 0.25), stringsAsFactors = FALSE)
  parquet_file <- tempfile(fileext = ".parquet")
  on.exit(unlink(parquet_file), add = TRUE)
  nanoparquet::write_parquet(original, parquet_file)

  scanned <- SEMseeker:::io_scan_parquet(parquet_file)
  expect_s3_class(scanned, "tbl_lazy")
  expect_identical(as.data.frame(dplyr::collect(scanned)), original)
})

test_that("closing twice is not an error, and a closed engine opens again", {
  engine_session()
  first_connection <- SEMseeker:::io_engine_connection()

  SEMseeker:::core_close_env()
  expect_false(DBI::dbIsValid(first_connection))
  expect_silent(SEMseeker:::io_engine_disconnect())

  engine_session()
  on.exit(SEMseeker:::core_close_env(), add = TRUE)
  second_connection <- SEMseeker:::io_engine_connection()
  expect_true(DBI::dbIsValid(second_connection))
})

test_that("another process gets its own connection, with its share of the limit and one thread", {
  ssEnv <- engine_session(memory_limit = "1GiB")
  on.exit(SEMseeker:::core_close_env(), add = TRUE)
  parent_connection <- SEMseeker:::io_engine_connection()
  # Once the cache points at the simulated worker, nothing closes the parent's.
  on.exit(if (DBI::dbIsValid(parent_connection)) DBI::dbDisconnect(parent_connection, shutdown = TRUE), add = TRUE)

  # A multisession worker is a different process holding a copy of the session:
  # the connection cached here is not its own, and the parent is the one that
  # opened the session. Both conditions are simulated by changing the recorded
  # process ids.
  assign("engine_process_id", -1L, envir = SEMseeker:::.pkgglobalenv)
  ssEnv$main_process_id <- -1L
  SEMseeker:::core_update_session_info(ssEnv)
  worker_connection <- SEMseeker:::io_engine_connection()

  expect_false(identical(worker_connection, parent_connection))
  expect_identical(as.integer(engine_setting(worker_connection, "threads")), 1L)
  expect_equal(amount_bytes(engine_setting(worker_connection, "memory_limit")), 1024^3 / ssEnv$parallel$nCore, tolerance = 0.01)
})

test_that("the engine writes nothing in the home directory", {
  fake_home <- tempfile("home_")
  dir.create(fake_home)
  withr::local_envvar(HOME = fake_home)
  engine_session()
  on.exit(SEMseeker:::core_close_env(), add = TRUE)

  parquet_file <- tempfile(fileext = ".parquet")
  on.exit(unlink(parquet_file), add = TRUE)
  SEMseeker:::io_sink_parquet(SEMseeker:::io_register_frame(data.frame(VALUE = 1:3)), parquet_file)
  dplyr::collect(SEMseeker:::io_scan_parquet(parquet_file))

  expect_false(dir.exists(file.path(fake_home, ".duckdb")))
})
