# Tests for utility helper functions (no session required)
#
# Covered:
#  - util_boolean_check()          string/logical/NA → logical coercion
#  - io_file_path_build()        path construction helper
#  - io_name_composer()          filename fragment assembler
#  - util_data_frame_add_column()  safe column add/update on data.frames
#  - assoc_filter_sql()             SQL-condition row filter (runs in DuckDB)

# ---------------------------------------------------------------------------
# 1. util_boolean_check
# ---------------------------------------------------------------------------

test_that("util_boolean_check returns FALSE for NA, empty string, NULL", {
  expect_false(SEMseeker:::util_boolean_check(NA))
  expect_false(SEMseeker:::util_boolean_check(""))
  expect_false(SEMseeker:::util_boolean_check(NULL))
})

test_that("util_boolean_check returns TRUE for TRUE, 'TRUE', 'T', 'true', 'True', '1', 1", {
  expect_true(SEMseeker:::util_boolean_check(TRUE))
  expect_true(SEMseeker:::util_boolean_check("TRUE"))
  expect_true(SEMseeker:::util_boolean_check("T"))
  expect_true(SEMseeker:::util_boolean_check("true"))
  expect_true(SEMseeker:::util_boolean_check("True"))
  expect_true(SEMseeker:::util_boolean_check("1"))
  expect_true(SEMseeker:::util_boolean_check(1))
})

test_that("util_boolean_check returns FALSE for FALSE, 'FALSE', '0', 0", {
  expect_false(SEMseeker:::util_boolean_check(FALSE))
  expect_false(SEMseeker:::util_boolean_check("FALSE"))
  expect_false(SEMseeker:::util_boolean_check("0"))
  expect_false(SEMseeker:::util_boolean_check(0))
})

test_that("util_boolean_check returns FALSE for arbitrary strings", {
  expect_false(SEMseeker:::util_boolean_check("yes"))
  expect_false(SEMseeker:::util_boolean_check("no"))
  expect_false(SEMseeker:::util_boolean_check("maybe"))
})

# ---------------------------------------------------------------------------
# 2. io_file_path_build
# ---------------------------------------------------------------------------

test_that("io_file_path_build assembles path correctly", {
  result <- SEMseeker:::io_file_path_build("/base", c("a", "b"), "csv")
  expect_true(grepl("^/base", result))
  expect_true(grepl("A_B\\.csv$", result))
})

test_that("io_file_path_build adds .gz when add_gz = TRUE", {
  result <- SEMseeker:::io_file_path_build("/base", c("sample"), "bed", add_gz = TRUE)
  expect_true(endsWith(result, ".gz"))
  expect_true(grepl("\\.bed\\.gz$", result))
})

test_that("io_file_path_build does not add .gz when add_gz = FALSE", {
  result <- SEMseeker:::io_file_path_build("/base", c("sample"), "csv", add_gz = FALSE)
  expect_false(endsWith(result, ".gz"))
})

test_that("io_file_path_build cleans the filename part (upper, underscore)", {
  result <- SEMseeker:::io_file_path_build("/tmp", c("hello world", "foo.bar"), "txt")
  fname <- basename(result)
  expect_true(grepl("HELLO_WORLD", fname))
  expect_true(grepl("FOO_BAR",     fname))
})

test_that("io_file_path_build collapses multiple filenames with underscore", {
  result <- SEMseeker:::io_file_path_build("/tmp", c("A", "B", "C"), "csv")
  expect_true(grepl("A_B_C", basename(result)))
})

test_that("io_file_path_build avoids double dots", {
  result <- SEMseeker:::io_file_path_build("/tmp", c("X"), "csv")
  expect_false(grepl("\\.\\.", result))
})

# ---------------------------------------------------------------------------
# 3. io_name_composer
# ---------------------------------------------------------------------------

test_that("io_name_composer joins fragments with underscore", {
  result <- SEMseeker:::io_name_composer("a", "b", "c")
  expect_equal(result, "a_b_c")
})

test_that("io_name_composer replaces spaces with underscore", {
  result <- SEMseeker:::io_name_composer("hello world", "foo")
  expect_false(grepl(" ", result))
  expect_true(grepl("_", result))
})

test_that("io_name_composer replaces colons with underscore", {
  result <- SEMseeker:::io_name_composer("key:value")
  expect_false(grepl(":", result))
})

test_that("io_name_composer collapses double underscores", {
  result <- SEMseeker:::io_name_composer("a", "", "b")
  expect_false(grepl("__", result))
})

test_that("io_name_composer avoids double dots", {
  result <- SEMseeker:::io_name_composer("file.name", ".csv")
  expect_false(grepl("\\.\\.", result))
})

test_that("io_name_composer handles a single fragment", {
  result <- SEMseeker:::io_name_composer("single")
  expect_equal(result, "single")
})

# ---------------------------------------------------------------------------
# 4. util_data_frame_add_column
# ---------------------------------------------------------------------------

test_that("util_data_frame_add_column adds a new column when absent", {
  df <- data.frame(x = 1:3)
  result <- SEMseeker:::util_data_frame_add_column(df, "y", c(4, 5, 6))
  expect_true("y" %in% colnames(result))
  expect_equal(result$y, c(4, 5, 6))
})

test_that("util_data_frame_add_column updates existing column", {
  df <- data.frame(x = 1:3, y = c(0, 0, 0))
  result <- SEMseeker:::util_data_frame_add_column(df, "y", c(7, 8, 9))
  expect_equal(result$y, c(7, 8, 9))
  expect_equal(ncol(result), 2)   # no extra column added
})

test_that("util_data_frame_add_column works on empty data.frame", {
  df <- data.frame()
  result <- SEMseeker:::util_data_frame_add_column(df, "col", 42)
  expect_true("col" %in% colnames(result))
})

test_that("util_data_frame_add_column preserves existing columns", {
  df <- data.frame(a = 1:3, b = 4:6)
  result <- SEMseeker:::util_data_frame_add_column(df, "c", 7:9)
  expect_true(all(c("a", "b", "c") %in% colnames(result)))
  expect_equal(result$a, 1:3)
})

# ---------------------------------------------------------------------------
# 5. assoc_filter_sql
# ---------------------------------------------------------------------------
# assoc_filter_sql runs the condition in DuckDB and calls core_log_event (requires
# a session). Tests that reach the query are wrapped with core_init_env / close_env.
#
# The table is registered on the engine connection under a name of its own and
# removed when the query ends, so nothing is written to the caller's frame.

test_that("assoc_filter_sql returns data unchanged for empty conditions", {
  df <- data.frame(x = 1:5, g = c("a","b","a","b","a"))
  result <- SEMseeker:::assoc_filter_sql(c(), df)
  expect_equal(nrow(result), nrow(df))
})

test_that("assoc_filter_sql filters rows by numeric condition", {
  tempFolder <- sem_test_folder()
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  df <- data.frame(x = 1:10)
  result <- SEMseeker:::assoc_filter_sql("x > 5", df)
  expect_true(all(result$x > 5))
  expect_equal(nrow(result), 5)

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})

test_that("assoc_filter_sql filters rows by character condition", {
  tempFolder <- sem_test_folder()
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  df <- data.frame(
    val = c(1, 2, 3),
    grp = c("case", "control", "case"),
    stringsAsFactors = FALSE
  )
  result <- SEMseeker:::assoc_filter_sql("grp = 'case'", df)
  expect_equal(nrow(result), 2)
  expect_true(all(result$grp == "case"))

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})

test_that("assoc_filter_sql applies multiple conditions sequentially", {
  tempFolder <- sem_test_folder()
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  df <- data.frame(x = 1:10, y = c(rep("A",5), rep("B",5)),
                   stringsAsFactors = FALSE)
  result <- SEMseeker:::assoc_filter_sql(c("x > 3", "x < 8"), df)
  expect_true(all(result$x > 3 & result$x < 8))

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})

test_that("assoc_filter_sql returns empty data.frame when no rows match", {
  tempFolder <- sem_test_folder()
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  df <- data.frame(x = 1:5)
  result <- SEMseeker:::assoc_filter_sql("x > 100", df)
  expect_equal(nrow(result), 0)

  SEMseeker:::core_close_env()
  unlink(tempFolder, recursive = TRUE)
})

# ---------------------------------------------------------------------------
# assoc_filter_sql and the user's workspace
#
# The table used to be registered in the global environment under a fixed name and
# removed afterwards. A package may not write to the workspace of whoever is using
# it, and the cost was not only formal: an object of that name belonging to the
# user was overwritten and then deleted, and the removal was not protected against
# an error, so the one case that left the workspace dirty was the case where the
# query had already failed.
# ---------------------------------------------------------------------------

test_that("assoc_filter_sql leaves an object of the user's alone", {
  tempFolder <- sem_test_folder()
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                            showprogress = showprogress, verbosity = verbosity)
  on.exit({ SEMseeker:::core_close_env(); unlink(tempFolder, recursive = TRUE) },
          add = TRUE)

  # Under the very name the function used to borrow.
  sentinel <- "this belongs to the caller"
  assign("ssdf_tmp", sentinel, envir = globalenv())
  on.exit(suppressWarnings(rm("ssdf_tmp", envir = globalenv())), add = TRUE)

  result <- SEMseeker:::assoc_filter_sql("x > 7", data.frame(x = 1:10))

  # The filter did its job...
  expect_equal(nrow(result), 3)
  # ...and the caller's object is still there, with its value.
  expect_true(exists("ssdf_tmp", envir = globalenv(), inherits = FALSE))
  expect_identical(get("ssdf_tmp", envir = globalenv(), inherits = FALSE), sentinel)
})

test_that("assoc_filter_sql writes nothing to the workspace, failing or not", {
  tempFolder <- sem_test_folder()
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                            showprogress = showprogress, verbosity = verbosity)
  on.exit({ SEMseeker:::core_close_env(); unlink(tempFolder, recursive = TRUE) },
          add = TRUE)

  skip_if(exists("ssdf_tmp", envir = globalenv(), inherits = FALSE),
          "the workspace already holds that name")

  # A query that succeeds.
  SEMseeker:::assoc_filter_sql("x > 2", data.frame(x = 1:5))
  expect_false(exists("ssdf_tmp", envir = globalenv(), inherits = FALSE))

  # And one that cannot: the column is not there. Reached only because nothing is
  # written in the first place, which is why there is nothing to clean up.
  try(suppressWarnings(
    SEMseeker:::assoc_filter_sql("no_such_column > 1", data.frame(x = 1:5))),
    silent = TRUE)
  expect_false(exists("ssdf_tmp", envir = globalenv(), inherits = FALSE))
})

test_that("a condition carrying a sub-statement finds the table", {
  tempFolder <- sem_test_folder()
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                            showprogress = showprogress, verbosity = verbosity)
  on.exit({ SEMseeker:::core_close_env(); unlink(tempFolder, recursive = TRUE) },
          add = TRUE)

  # TABLE is the placeholder the function documents by refusing anything else. It
  # was replaced with the name of a local variable rather than the name the table
  # is registered under, so this branch asked for a table nobody had registered.
  df <- data.frame(x = 1:10)
  result <- SEMseeker:::assoc_filter_sql("x IN (SELECT x FROM TABLE WHERE x > 6)", df)

  expect_equal(sort(result$x), 7:10)
})

# ---------------------------------------------------------------------------
# assoc_filter_sql in the DuckDB dialect
#
# Expected rows are written by hand from the input. Cases 'LIKE' and 'division'
# give a different answer in SQLite, the engine sqldf used, so they also tell
# which engine ran the query.
# ---------------------------------------------------------------------------

filter_sql_session <- function(env = parent.frame()) {
  tempFolder <- sem_test_folder()
  SEMseeker:::core_init_env(tempFolder, parallel_strategy = parallel_strategy,
                            showprogress = showprogress, verbosity = verbosity)
  withr::defer({ SEMseeker:::core_close_env(); unlink(tempFolder, recursive = TRUE) },
               envir = env)
}

test_that("AND and OR combine comparisons", {
  filter_sql_session()
  df <- data.frame(x = 1:10, g = rep(c("a", "b"), 5), stringsAsFactors = FALSE)
  expect_identical(SEMseeker:::assoc_filter_sql("x > 2 AND g = 'a'", df)$x, c(3L, 5L, 7L, 9L))
  expect_identical(SEMseeker:::assoc_filter_sql("x < 2 OR x > 9", df)$x, c(1L, 10L))
})

test_that("IN keeps the rows whose value is in the list", {
  filter_sql_session()
  df <- data.frame(id = 1:4, g = c("a", "b", "c", "a"), stringsAsFactors = FALSE)
  expect_identical(SEMseeker:::assoc_filter_sql("g IN ('a', 'c')", df)$id, c(1L, 3L, 4L))
})

test_that("LIKE is case-sensitive", {
  filter_sql_session()
  # SQLite's LIKE ignores case and kept both "case" and "Case".
  df <- data.frame(g = c("case", "Case", "control"), stringsAsFactors = FALSE)
  expect_identical(SEMseeker:::assoc_filter_sql("g LIKE 'ca%'", df)$g, "case")
})

test_that("IS NULL matches the missing values", {
  filter_sql_session()
  df <- data.frame(id = 1:4, y = c(1, NA, 3, NA))
  expect_identical(SEMseeker:::assoc_filter_sql("y IS NULL", df)$id, c(2L, 4L))
  expect_identical(SEMseeker:::assoc_filter_sql("y IS NOT NULL", df)$id, c(1L, 3L))
})

test_that("dividing two integers keeps the fraction", {
  filter_sql_session()
  # 5 / 2 is 2.5 > 2 here; SQLite truncated it to 2 and dropped the row.
  df <- data.frame(x = 1:10)
  expect_identical(SEMseeker:::assoc_filter_sql("x / 2 > 2", df)$x, 5:10)
})

test_that("a condition with double quotes is refused before it runs", {
  filter_sql_session()
  df <- data.frame(grp = c("case", "control"), stringsAsFactors = FALSE)
  expect_error(SEMseeker:::assoc_filter_sql('grp = "case"', df), "single quotes")
})

test_that("rows come back in input order, sub-statement included", {
  filter_sql_session()
  df <- data.frame(x = c(5L, 1L, 9L, 3L, 7L))
  expect_identical(SEMseeker:::assoc_filter_sql("x IN (SELECT x FROM TABLE WHERE x > 2)", df)$x,
                   c(5L, 9L, 3L, 7L))
  expect_identical(SEMseeker:::assoc_filter_sql("x > 2", df)$x, c(5L, 9L, 3L, 7L))
})

test_that("column types survive the filter", {
  filter_sql_session()
  df <- data.frame(chr = c("a", "b"), int = 1:2, dbl = c(0.5, 1.5), lgl = c(TRUE, FALSE),
                   stringsAsFactors = FALSE)
  result <- SEMseeker:::assoc_filter_sql("int > 0", df)
  expect_identical(vapply(result, class, character(1)), vapply(df, class, character(1)))
  expect_identical(result, df)
})

test_that("the filter does not load tcltk", {
  filter_sql_session()
  skip_if("tcltk" %in% loadedNamespaces(), "tcltk was loaded before the filter ran")
  SEMseeker:::assoc_filter_sql("x > 1", data.frame(x = 1:3))
  expect_false("tcltk" %in% loadedNamespaces())
})

test_that("the filter leaves no table registered on the engine", {
  filter_sql_session()
  connection <- SEMseeker:::io_engine_connection()
  tables_before <- DBI::dbListTables(connection)
  SEMseeker:::assoc_filter_sql("x > 1", data.frame(x = 1:3))
  try(SEMseeker:::assoc_filter_sql("no_such_column > 1", data.frame(x = 1:3)), silent = TRUE)
  expect_identical(DBI::dbListTables(connection), tables_before)
})
