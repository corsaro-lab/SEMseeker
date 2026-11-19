# ---------------------------------------------------------------------------
# io_engine - the DuckDB connection of the process and the three primitives
# that replace polars: scan a parquet file, register a data frame, sink a
# query to parquet.
#
# Why DuckDB: polars is on neither CRAN nor Bioconductor, and Bioconductor
# refuses the extra repository that makes it installable. DuckDB is on CRAN,
# reads and writes parquet natively and is queried with the dplyr verbs through
# dbplyr, so parquet stays the format on disk and nothing already written by
# SEMseeker has to be converted.
#
# Four facts measured before this file was written shape it:
#   - memory_limit and temp_directory passed to dbConnect(config = ) are
#     ignored without a warning: they are set with SET and read back;
#   - DuckDB exceeds its memory limit by 70-130%, because the limit governs
#     only part of the memory it uses: the default is half of the physical
#     memory, not DuckDB's own 80%;
#   - without shared_home = FALSE DuckDB keeps extensions and secrets in
#     ~/.duckdb, and a Bioconductor package writes nowhere outside tempdir()
#     and the folders the user names;
#   - the default parquet compression is snappy, about a third larger than
#     the zstd files polars wrote.
#
# One connection per process, kept in .pkgglobalenv and not in the session:
# the session is serialised to multisession workers, which are separate
# processes and cannot use a connection opened by the parent. Each process
# opens its own on first use. A worker takes a share of the memory limit and
# one thread, because the workers run at the same time on the same machine.
# ---------------------------------------------------------------------------

#' The DuckDB connection of the current process, opened on first use.
#'
#' The connection is reused while it belongs to this process, is still valid
#' and was opened for the current session temporary folder; otherwise a new
#' one is opened. The memory limit, the threads and the spill directory are
#' set with SET. The threads and the spill directory are read back and a
#' mismatch stops here; the memory limit cannot be compared with what was
#' asked (DuckDB resolves "50%" into an amount), so the amount in force is
#' written to the log and test-io_engine.R checks it.
#'
#' The settings come from the session: `memory_limit` (session option,
#' default "50%" of the physical memory, as DuckDB detects it on every
#' system), the number of workers `parallel$nCore` for the threads, and
#' `temp_folder` for the spill directory. Without an open session the
#' defaults of the session options apply and the spill directory goes under
#' tempdir().
#'
#' @return an open DBIConnection to an in-memory DuckDB database
#' @noRd
io_engine_connection <- function() {
  ssEnv <- tryCatch(core_get_session_info(), error = function(e) list())
  temp_folder <- if (is.null(ssEnv$temp_folder)) file.path(tempdir(), "semseeker_engine") else ssEnv$temp_folder

  connection <- .pkgglobalenv$engine_connection
  if (!is.null(connection) && identical(.pkgglobalenv$engine_process_id, Sys.getpid()) &&
      identical(.pkgglobalenv$engine_temp_folder, temp_folder) && DBI::dbIsValid(connection))
    return(connection)

  # A connection cached by another process (the parent whose session a worker
  # received) is left alone: it is not this process's to close.
  if (!is.null(connection) && identical(.pkgglobalenv$engine_process_id, Sys.getpid()))
    io_engine_disconnect()

  memory_limit <- if (is.null(ssEnv$memory_limit)) .SS_DEFAULTS$memory_limit$value else ssEnv$memory_limit
  worker_count <- if (is.null(ssEnv$parallel$nCore)) 1L else as.integer(ssEnv$parallel$nCore)
  is_worker <- !is.null(ssEnv$main_process_id) && !identical(as.integer(ssEnv$main_process_id), Sys.getpid())
  thread_count <- if (is_worker) 1L else max(1L, worker_count)

  spill_directory <- file.path(temp_folder, "duckdb")
  dir.create(spill_directory, recursive = TRUE, showWarnings = FALSE)

  connection <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  # SET and current_setting() are the only ways to change and read a setting:
  # dbConnect(config = ) ignores them (see the header of this file).
  # setting_name: a DuckDB setting such as "memory_limit"; setting_value:
  # its value, quoted as text (DuckDB parses "50%", "8GB", "4").
  set_value <- function(setting_name, setting_value)
    DBI::dbExecute(connection, sprintf("SET %s = %s", setting_name, DBI::dbQuoteString(connection, as.character(setting_value))))
  # setting_name: the setting to read back; returns its value as DuckDB
  # reports it ("32.0 GiB", 4, a path).
  read_value <- function(setting_name)
    DBI::dbGetQuery(connection, sprintf("SELECT current_setting('%s') AS setting_value", setting_name))$setting_value

  tryCatch({
    set_value("temp_directory", spill_directory)
    set_value("memory_limit", memory_limit)
    if (is_worker && worker_count > 1L)
      set_value("memory_limit", sprintf("%.0fB", io_engine_amount_bytes(read_value("memory_limit")) / worker_count))
    set_value("threads", thread_count)
  }, error = function(e) {
    DBI::dbDisconnect(connection, shutdown = TRUE)
    stop("SEMseeker: the data engine refused its settings (memory_limit = '", memory_limit, "'): ", conditionMessage(e), call. = FALSE)
  })

  if (!identical(as.integer(read_value("threads")), thread_count) ||
      !identical(normalizePath(read_value("temp_directory"), mustWork = FALSE), normalizePath(spill_directory, mustWork = FALSE))) {
    DBI::dbDisconnect(connection, shutdown = TRUE)
    stop("SEMseeker: the data engine did not apply its settings.", call. = FALSE)
  }

  core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " data engine opened: memory_limit ", read_value("memory_limit"),
                 ", threads ", thread_count, ", spill to ", spill_directory)

  assign("engine_connection", connection, envir = .pkgglobalenv)
  assign("engine_process_id", Sys.getpid(), envir = .pkgglobalenv)
  assign("engine_temp_folder", temp_folder, envir = .pkgglobalenv)
  connection
}

#' Close the DuckDB connection of this process, if there is one.
#'
#' Closing what is not open is not an error, as for core_close_env(), which
#' calls this before it removes the session temporary folder that holds the
#' spill directory. A connection cached by another process is forgotten, not
#' closed.
#'
#' @return NULL, invisibly
#' @noRd
io_engine_disconnect <- function() {
  connection <- .pkgglobalenv$engine_connection
  if (!is.null(connection) && identical(.pkgglobalenv$engine_process_id, Sys.getpid()) && DBI::dbIsValid(connection))
    DBI::dbDisconnect(connection, shutdown = TRUE)
  for (name in c("engine_connection", "engine_process_id", "engine_temp_folder"))
    if (exists(name, envir = .pkgglobalenv, inherits = FALSE)) rm(list = name, envir = .pkgglobalenv)
  invisible(NULL)
}

#' Bytes of a memory amount as DuckDB prints it.
#'
#' DuckDB reports its settings as text ("32.0 GiB"); the share of a worker is
#' computed from the limit DuckDB resolved, so a percentage and an absolute
#' amount are divided the same way.
#'
#' @param amount_text text such as "32.0 GiB", "500.0 MiB" or "0 bytes"
#' @return the amount in bytes, as a double
#' @noRd
io_engine_amount_bytes <- function(amount_text) {
  parts <- strsplit(trimws(amount_text), " +")[[1]]
  unit_bytes <- c(bytes = 1, KiB = 1024, MiB = 1024^2, GiB = 1024^3, TiB = 1024^4, PiB = 1024^5)
  if (length(parts) != 2L || !(parts[2] %in% names(unit_bytes)))
    stop("SEMseeker: unexpected memory amount from the data engine: '", amount_text, "'", call. = FALSE)
  as.numeric(parts[1]) * unit_bytes[[parts[2]]]
}

#' A lazy table over one or more parquet files.
#'
#' Replaces polars' scan_parquet(): nothing is read until the query is
#' collected or sunk, and dplyr verbs applied to the result are pushed down to
#' DuckDB.
#'
#' @param parquet_paths character vector of parquet files read as one table;
#'   they must share their columns
#' @return a tbl_lazy on the connection of this process
#' @noRd
io_scan_parquet <- function(parquet_paths) {
  connection <- io_engine_connection()
  path_list <- paste(DBI::dbQuoteString(connection, normalizePath(parquet_paths, mustWork = TRUE)), collapse = ", ")
  dplyr::tbl(connection, dbplyr::sql(sprintf("SELECT * FROM read_parquet([%s])", path_list)))
}

#' Make an R data frame queryable together with the parquet files.
#'
#' Replaces polars' as_polars_df(frame)$lazy(), used to join annotations and
#' sample sheets with the pivots. The frame is registered as a view, without
#' copying it into the database, for the life of the connection.
#'
#' @param frame a data.frame
#' @return a tbl_lazy on the connection of this process
#' @noRd
io_register_frame <- function(frame) {
  connection <- io_engine_connection()
  view_name <- paste0("frame_", gsub("[^A-Za-z0-9]", "", basename(tempfile(""))))
  duckdb::duckdb_register(connection, view_name, as.data.frame(frame))
  dplyr::tbl(connection, view_name)
}

#' Write a lazy table to a parquet file without bringing it into R.
#'
#' Replaces polars' sink_parquet(). The compression is declared, zstd,
#' because DuckDB's default (snappy) wrote files about a third larger than
#' the ones polars wrote.
#'
#' @param lazy_table a tbl_lazy on the connection of this process
#' @param parquet_path the file to write; an existing file is replaced
#' @return parquet_path, invisibly
#' @noRd
io_sink_parquet <- function(lazy_table, parquet_path) {
  connection <- dbplyr::remote_con(lazy_table)
  query_text <- dbplyr::sql_render(lazy_table)
  DBI::dbExecute(connection, sprintf("COPY (%s) TO %s (FORMAT parquet, COMPRESSION zstd)",
                                     query_text, DBI::dbQuoteString(connection, parquet_path)))
  invisible(parquet_path)
}

#' Write an R data frame to a parquet file through the data engine.
#'
#' Replaces polars' as_polars_df(frame)$write_parquet(path). The frame is
#' registered only for the copy and unregistered afterwards, so the engine
#' does not keep a reference to a frame that may be the whole signal matrix.
#'
#' @param frame a data.frame; its column types are kept (integer as INTEGER,
#'   double as DOUBLE, character as VARCHAR)
#' @param parquet_path the file to write; an existing file is replaced
#' @return parquet_path, invisibly
#' @noRd
io_write_frame_parquet <- function(frame, parquet_path) {
  connection <- io_engine_connection()
  view_name <- paste0("frame_", gsub("[^A-Za-z0-9]", "", basename(tempfile(""))))
  duckdb::duckdb_register(connection, view_name, as.data.frame(frame))
  on.exit(duckdb::duckdb_unregister(connection, view_name), add = TRUE)
  DBI::dbExecute(connection, sprintf("COPY (SELECT * FROM %s) TO %s (FORMAT parquet, COMPRESSION zstd)",
                                     DBI::dbQuoteIdentifier(connection, view_name), DBI::dbQuoteString(connection, parquet_path)))
  invisible(parquet_path)
}

#' Read a whole parquet file into an R data frame through the data engine.
#'
#' Replaces as.data.frame(polars::pl$read_parquet(path)) for the small files
#' read whole (thresholds, cleaned inputs).
#'
#' @param parquet_path the file to read
#' @return a data.frame
#' @noRd
io_read_parquet_frame <- function(parquet_path) {
  connection <- io_engine_connection()
  DBI::dbGetQuery(connection, sprintf("SELECT * FROM read_parquet(%s)", DBI::dbQuoteString(connection, parquet_path)))
}
