#' Keep the rows of a data frame that satisfy SQL WHERE conditions.
#'
#' Users select samples, areas and results with a WHERE fragment written in
#' their own inference details, so the filter speaks SQL rather than R. The
#' fragment runs in DuckDB, the data engine of the package. It ran in SQLite
#' through sqldf until 0.99.7; sqldf brought tcltk, and with it XQuartz on
#' macOS, for nothing else.
#'
#' The dialect is DuckDB's, and three things differ from SQLite: double quotes
#' name a column instead of delimiting a string, LIKE is case-sensitive (ILIKE
#' is not), and dividing two integers gives a fractional result. The first is
#' refused at the door, because a string in double quotes would otherwise be
#' read as a column and fail with a message about a missing column.
#'
#' @param sql_conditions character vector of WHERE fragments, applied one after
#'   the other; an empty vector, NA or "" leaves the data unchanged. A fragment
#'   may hold a sub-statement on the same data, written FROM TABLE.
#' @param data_frame data.frame to filter. Columns whose names differ only by
#'   case are reduced to the first, because SQL does not tell them apart.
#' @return the rows of data_frame satisfying every condition, in their input
#'   order, with the column types of the input.
#' @noRd
assoc_filter_sql <- function(sql_conditions, data_frame)
{

  # check if sql_conditions is am empty string
  if (length(sql_conditions) == 0 )
    return(data_frame)

  for ( s in seq_along(sql_conditions))
  {
    sql_condition <- as.character(sql_conditions[s])
    sql_condition <- trimws(sql_condition)
    if(!is.na(sql_condition))
      if(!is.null(sql_condition))
        if (sql_condition != "")
        {
          # sql_condition <- toupper(sql_condition)
          sql_condition <- gsub("  ", " ", sql_condition)
          # In DuckDB a double-quoted word is a column name. SQLite accepted it
          # as a string, so conditions written for sqldf would now fail, or
          # compare two columns, without saying why.
          if (grepl("\"", sql_condition, fixed = TRUE))
            stop("The SQL condition '", sql_condition, "' contains double quotes. ",
                 "Write strings in single quotes, e.g. Sample_Group = 'Case'.", call. = FALSE)
          if (grepl("FROM", sql_condition))
          {
            if (!grepl("FROM TABLE", sql_condition))
            {
              core_log_event("ERROR: ",format(Sys.time(), "%a %b %d %X %Y")," the only allowed table placeholder with sub statement is TABLE!")
              stop()
            }
          }
          columns <- toupper(colnames(data_frame))
          # drop duplicates columns
          data_frame <- data_frame[, !duplicated(columns), drop = FALSE]
          data_frame <- assoc_filter_sql_run(sql_condition, data_frame)
          core_log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y"), " Executed sql: " , sql_condition)
        }
  }
  return(data_frame)
}

#' Run one WHERE fragment on a data frame in the data engine.
#'
#' Kept apart from assoc_filter_sql() so that the registration of the table
#' and its removal sit in one frame, with on.exit() guaranteeing the removal
#' when the query fails.
#'
#' The table gets a row ordinal of its own, and the query sorts on it, because
#' a sub-statement is planned as a join and joins do not keep row order. The
#' ordinal is dropped from the result with EXCLUDE.
#'
#' @param sql_condition one WHERE fragment, already checked by the caller;
#'   TABLE stands for the data itself.
#' @param data_frame data.frame without columns differing only by case.
#' @return the matching rows as a data.frame, in input order.
#' @noRd
assoc_filter_sql_run <- function(sql_condition, data_frame)
{
  connection <- io_engine_connection()
  view_name <- paste0("filter_", gsub("[^A-Za-z0-9]", "", basename(tempfile(""))))
  ordinal_column <- paste0(view_name, "_row")

  registered_frame <- as.data.frame(data_frame)
  registered_frame[[ordinal_column]] <- seq_len(nrow(registered_frame))
  duckdb::duckdb_register(connection, view_name, registered_frame)
  on.exit(duckdb::duckdb_unregister(connection, view_name), add = TRUE)

  # TABLE is the placeholder users write for a sub-statement on the same data.
  sql_condition <- gsub("FROM TABLE", paste("FROM", view_name), sql_condition, fixed = TRUE)
  sql <- sprintf("SELECT * EXCLUDE (%s) FROM %s WHERE %s ORDER BY %s",
                 ordinal_column, view_name, sql_condition, ordinal_column)
  result <- DBI::dbGetQuery(connection, sql)
  rownames(result) <- NULL
  result
}
