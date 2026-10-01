core_set_env_variable <- function(arguments, var_name, var_default_value, check_values=NULL)
{


  ssEnv <- core_get_session_info()
  if(is.null(ssEnv[[var_name]]))
    ssEnv[[var_name]]  <- var_default_value

  if(!is.null(arguments[[var_name]]))
  {
    if(any(arguments[[var_name]]=="DEFAULT"))
      ssEnv[[var_name]]  <- var_default_value
    else
      ssEnv[[var_name]]  <- arguments[[var_name]]
  }
  # In memory only. This runs once per option, about twenty times per
  # core_init_env(), and core_update_session_info() writes the WHOLE session
  # twice per call. That was ~40 writes of the entire ssEnv per initialisation,
  # and it is what made core_init_env() take 287 s once the probe annotation had
  # been built, against 0.8 s before it. The annotation no longer lives in ssEnv,
  # so the file is small again, but writing it per option was never the right
  # shape: the caller flushes once, when it is done setting options.
  core_update_session_info(ssEnv, save_to_disk = FALSE)
  core_log_event("INFO:",format(Sys.time(), "%a %b %d %X %Y"), " " , var_name ," ",  ssEnv[[var_name]])

  if(!is.null(check_values))
  {
    values <- arguments[[var_name]]
    unknown_options <- values[!values %in% check_values]
  }
  else
  {
    unknown_options <- NULL
  }

  # remove from arguments var_name
  arguments[[var_name]] <- unknown_options

  return(arguments)
}
