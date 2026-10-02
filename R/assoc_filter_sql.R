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
          if (grepl("FROM", sql_condition))
          {
            if (!grepl("FROM TABLE", sql_condition))
            {
              core_log_event("ERROR: ",format(Sys.time(), "%a %b %d %X %Y")," the only allowed table placeholder with sub statement is TABLE!")
              stop()
            }
            else
              # The same name the table is registered under below. It used to be
              # replaced with "data_frame", which is the name of the local
              # variable and of no table sqldf can see, so a condition carrying a
              # sub-statement asked for a table that does not exist. Nothing
              # exercises that branch - no caller passes one and no test covers
              # it - which is how it stayed broken.
              sql_condition <- gsub("TABLE", "ssdf_tmp", sql_condition)
          }
          columns <- toupper(colnames(data_frame))
          # drop duplicates columns
          data_frame <- data_frame[, !duplicated(columns), drop = FALSE]
          columns <- colnames(data_frame)
          # sqldf looks its tables up by name in the environment it is given, so
          # it is given one of its own, holding nothing else and living only as
          # long as this query.
          #
          # It used to be handed the user's workspace: the table was assigned
          # there under a fixed name and removed afterwards. A package may not
          # write to the workspace of whoever is using it - R CMD check reports it
          # and Bioconductor refuses it - and the cost was not only formal. An
          # object of that name belonging to the user was overwritten and then
          # deleted, and the removal was not protected, so the one case that left
          # the workspace dirty was the case where the query had failed.
          #
          # Empty parent on purpose: the only name this environment resolves is
          # the table, so nothing of the caller's can be mistaken for one.
          query_env <- new.env(parent = emptyenv())
          assign("ssdf_tmp", data_frame, envir = query_env)
          sql <- paste("select * from ssdf_tmp where ", sql_condition)
          data_frame <- sqldf::sqldf(sql, envir = query_env)
          core_log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y"), " Executed sql: " , sql_condition)
        }
  }
  return(data_frame)
}
