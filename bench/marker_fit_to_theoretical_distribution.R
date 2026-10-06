#' @importFrom doRNG %dorng%
#' @importFrom doFuture %dofuture%
marker_fit_to_theoretical_distribution <- function()
{


  # result_folder, maxResources = 90, parallel_strategy  = "multisession", ...

  # ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)
  ssEnv <- get_session_info()

  keys <- ssEnv$keys_areas_subareas_markers_figures
  keys <- keys[keys$AREA == "PROBE",]

  remove_zeros <- c(TRUE,FALSE)
  to_export <- c("keys")
  nkeys <- nrow(keys)

  # sort key desc for marker
  keys <- keys[order(keys$MARKER, decreasing = TRUE),]

  res <- data.frame()
  for (rz in 1:2)
  {
    if(ssEnv$showprogress)
      progress_bar <- progressr::progressor(along = 1:nkeys)
    else
      progress_bar <- ""

    result_temp <- foreach::foreach(k = 1:nkeys, .combine =  plyr::rbind.fill, .export = to_export) %dorng%
      # for (k in 1:nkeys)
      {
        # in-memory only (workers need .pkgglobalenv populated,
        # but no need to hit disk per key).
        update_session_info(ssEnv, save_to_disk = FALSE)
        # k <- 1
        key <- keys [k,]
        log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y"), "Processing key: ", key$MARKER," ", key$FIGURE," ", key$AREA," ", key$SUBAREA)
        if(key$EXT == "bed")
          distributions <- list(
            Poisson = list(name = "pois"),
            NegativeBinomial = list(name = "nbinom"),
            Binomial = list(name = "binom"),
            Geometric = list(name = "geom")
          )
        else
          distributions <- list(
            Normal = list(name = "norm"),
            Exponential = list(name = "exp"),
            Lognormal = list(name = "lnorm"),
            Gamma = list(name = "gamma"),
            Weibull = list(name = "weibull"),
            Beta = list(name = "beta", start = list(shape1 = 2, shape2 = 2))
          )

        res_temp <- data.frame("MARKER" = key$MARKER, "FIGURE" = key$FIGURE, "AREA" = key$AREA, "SUBAREA" = key$SUBAREA)
        res_temp$ZERO_REMOVED <- remove_zeros[rz]

        pivot_subfolder <- dir_check_and_create(result_folderPivot, key$MARKER)
        fname <- file_path_build( pivot_subfolder ,c(key$MARKER, key$FIGURE, key$AREA,key$SUBAREA),"csv", add_gz=TRUE)

        data <- utils::read.csv2(gzfile(fname), header  =  TRUE)
        data <- data[-1, -1]
        data <- as.numeric(as.vector(as.matrix(data)))
        if(remove_zeros[rz])
          data <- data[data!=0]

        table(data==0)

        fits <- lapply(distributions, function(dist) {
          fit_distribution(dist$name, data, start = dist$start)
        })


        # Remove NULLs from results (failed fits)
        fits <- Filter(Negate(is.null), fits)

        if(length(fits) == 0)
        {
          aic_values <- c()
          log_event("WARNING: ",format(Sys.time(), "%a %b %d %X %Y")," No distribution fits found for marker ", key$MARKER, " figure ", key$FIGURE, " area ", key$AREA, " subarea ", key$SUBAREA)
          res_temp$best_fit <- "none"
        }
        else
        {
          # Compute AIC for each distribution
          aic_values <- sapply(fits, function(fit) { fit$distname = round(fit$aic,0)})
          aic_values <- t(as.data.frame(aic_values))

          # store AIC values
          if(length(aic_values) > 0)
            res_temp <- cbind (res_temp, aic_values)

          # Select the distribution with the lowest AIC
          best_fit <- names(which.min(aic_values))
          res_temp$best_fit <- best_fit

          # Print estimated parameters for the best-fit distribution
          best_fit_dist <- fits[[which.min(aic_values)]]
          res_temp$best_fit_dist <- best_fit_dist$distname

          # Obtain theoretical quantiles for the best-fit distribution
          # theoretical_quantiles <- get_theoretical_quantiles(best_fit_dist, data, best_fit)

          # Check for NA values in theoretical quantiles
          # if (any(is.na(theoretical_quantiles))) {
          #   cat("Warning: theoretical quantiles contain NA values.\n")
          #   table(is.na(theoretical_quantiles))
          # } else
          # {
          #   # Crea il QQ-plot
          #   qqplot(theoretical_quantiles, data, main = paste("QQ-plot: Distribuzione", best_fit), xlab = "Quantili Teorici", ylab = "Quantili Campionari")
          #   qqline(data, distribution = function(p) get_theoretical_quantiles(best_fit_dist, data, best_fit), col = "red")
          # }
        }

        colors_dist <- c("magenta", "red", "green", "purple", "orange", "brown")
        chartFolder <- dir_check_and_create(ssEnv$result_folderChart,c("Distributions"))
        if(remove_zeros[rz])
          filename  =  file_path_build(chartFolder,c(key$MARKER,key$FIGURE,key$AREA,key$SUBAREA, "HISTOGRAM","ZERO","REMOVED"),ssEnv$plot_format)
        else
          filename  =  file_path_build(chartFolder,c(key$MARKER,key$FIGURE,key$AREA,key$SUBAREA, "HISTOGRAM"),ssEnv$plot_format)
        # filename  =  file_path_build(chartFolder,c(key$MARKER,key$FIGURE,key$AREA,key$SUBAREA, "HISTOGRAM"),"png")

        # data <- data.frame(VALUE = data)
        # # Plotting KDE with overlayed histogram
        # p <- ggplot2::ggplot(data, ggplot2::aes(x = VALUE)) +
        #   ggplot2::geom_histogram(ggplot2::aes(y = ..density..), bins = 30, fill = "grey", alpha = 0.5, color = "black") +
        #   ggplot2::geom_density(adjust = 0.5, fill = "blue", alpha = 0.3) +
        #   ggplot2::labs(title = "Histogram with theoretical distribution curves", x = "Value", y = "Density") +
        #   ggplot2::theme_minimal()
        #
        #
        # if(length(fits) > 1)
        # {
        #   # Adding theoretical density curves
        #   for (i in seq_along(fits)) {
        #     p <- p + add_density_curve(data, fits[[i]]$distname, colors_dist[i])
        #   }
        #
        #   p <- p + ggplot2::scale_color_manual(
        #     name = "Theoretical Densities",
        #     values = setNames(colors_dist[seq_along(fits)], sapply(fits, `[[`, "distname"))
        #   )
        # }

        # Save the plot to a file
        # ggplot2::ggsave(filename, plot = p, width = 2480 / ssEnv$plot_resolution, height = 2480 / ssEnv$plot_resolution, dpi = as.numeric(ssEnv$plot_resolution_ppi))

        if(ssEnv$plot_format == "png")
          grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15)
        if(ssEnv$plot_format == "eps")
          grDevices::postscript(file =  filename, width = 2480,height = 2480, pointsize  =  15)
        graphics::hist(data, breaks = 100, probability = FALSE, col = "gray", main = "Histogram with theoretical distribution curves")
        # Add theoretical density curves to the histogram

        legend_labels <- character(0)
        if(length(fits)>0)
          for (fd in seq_along(fits))
          {
            # fd <-2
            if(length(aic_values)>0)
            {
              add_density_curve(fits[[fd]], fits[[fd]]$distname, colors_dist[fd])
              # create labels with dist name and aic
              legend_labels <- c(legend_labels, paste(fits[[fd]]$distname, ifelse(fd==which.min(aic_values),"*","")))
            }
            else
              legend_labels <- c(legend_labels, fits[[fd]]$distname)
            # Legend for the plot
            graphics::legend("topright", legend = legend_labels, col = colors_dist, lwd = 2)
          }
        grDevices::dev.off()

        # x <- seq(min(data), max(data), length = 100)
        # Aggiungi etichette AIC evitando sovrapposizioni
        # for (i in seq_along(fits)) {
        #   if (!is.null(fits[[i]])) {
        #     offsets <- seq(-0.1, 0.1, length.out = length(fits))
        #     text(x = 0.5, y = 0.5 + offsets[i],
        #       labels = paste(names(distributions)[i], "AIC:", round(fits[[i]]$aic, 2)),
        #       col = colors_dist[i], pos = 4)
        #   }
        # }

        if(ssEnv$showprogress)
          progress_bar(sprintf("Doing fitting: %s", stringr::str_pad(key$MARKER, 10, pad = " ")))

        # result_temp <- plyr::rbind.fill(result_temp, res_temp)
        res_temp
      }
  }
  # end-of-foreach disk snapshot.
  update_session_info(ssEnv, save_to_disk = TRUE)

  # res <- plyr::rbind.fill(res, res_temp)
  dataFolder <- dir_check_and_create(ssEnv$result_folderData,c("Distributions"))
  # close_env()
  filename  =  file_path_build(dataFolder,c("SUMMARY", "HISTOGRAM"),"csv")
  utils::write.csv2(result_temp, file = filename, row.names = FALSE)
}

# Helper function to add theoretical distribution curves
add_density_curve <- function(fit, dist_name, col) {
  graphics::curve(do.call(paste0("d", dist_name), c(list(x = x), as.list(fit$estimate))),col = col, lwd = 2, add = TRUE)
}

#' Fit a parametric distribution to data
#'
#' Wrapper around \code{\link[fitdistrplus]{fitdist}} that silently returns
#' \code{NULL} on failure instead of stopping.  Used internally to try
#' multiple candidate distributions and keep only the ones that converge.
#'
#' @param dist_name Character scalar: distribution name as accepted by
#'   \code{fitdistrplus::fitdist} (e.g. \code{"norm"}, \code{"beta"},
#'   \code{"gamma"}).
#' @param data Numeric vector of observations to fit.
#' @param start Named list of initial parameter values passed to
#'   \code{fitdistrplus::fitdist}.  \code{NULL} lets the function choose
#'   defaults.
#'
#' @return A \code{fitdist} object on success, or \code{NULL} if fitting
#'   fails.
#'
# Fit a distribution with exception handling
fit_distribution <- function(dist_name, data, start = NULL) {
  tryCatch({
    fit <- fitdistrplus::fitdist(data, dist_name, start = start)
    return(fit)
  }, error = function(e) {
    return(NULL)
  })
}

# Get theoretical quantiles for the fitted distribution
get_theoretical_quantiles <- function(fit, data, dist_name) {
  # fit <- best_fit_dist
  # dist_name <- best_fit
  dist_name <- tolower(dist_name)
  probs <- ppoints(length(data))
  if (dist_name == "beta") {
    qtheo <- stats::qbeta(probs, shape1 = fit$estimate["shape1"], shape2 = fit$estimate["shape2"])
  } else if (dist_name == "gamma") {
    qtheo <- stats::qgamma(probs, shape = fit$estimate["shape"], rate = fit$estimate["rate"])
  } else if (dist_name == "weibull") {
    qtheo <- stats::qweibull(probs, shape = fit$estimate["shape"], scale = fit$estimate["scale"])
  } else if (dist_name == "lnorm") {
    qtheo <- stats::qlnorm(probs, meanlog = fit$estimate["meanlog"], sdlog = fit$estimate["sdlog"])
  } else if (dist_name == "exp") {
    qtheo <- stats::qexp(probs, rate = fit$estimate["rate"])
  } else {
    qtheo <- stats::qnorm(probs, mean = fit$estimate["mean"], sd = fit$estimate["sd"])
  }
  return(qtheo)
}


# add_density_curve <- function(fit, dist_name, col) {
#   dist_fun <- match.fun(paste0("d", dist_name))
#   stat_function(
#     fun = dist_fun,
#     args = as.list(fit$estimate),
#     ggplot2::aes(color = dist_name),
#     size = 1
#   )
# }
