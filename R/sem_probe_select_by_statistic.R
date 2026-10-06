#' The probe at a given point of the cohort burden distribution
#'
#' Sums a marker across the samples of a probe pivot, giving one total per probe,
#' and returns the name of the probe sitting at a requested point of that
#' distribution: its minimum, first quartile, median, third quartile or maximum.
#'
#' It exists so that a chart can be asked for a representative probe without the
#' chart having to decide what representative means. Through 0.99.6 the two were
#' one function: asking for a chart scanned every marker, computed all five
#' statistics, wrote them to `PROBES_STAT.csv` and drew ten charts, so there was
#' no way to ask for one chart of one probe. Separating them costs a second call
#' and buys a chart whose subject is an argument.
#'
#' @section Why a quartile returns a probe and not a value:
#' The quantile of a distribution of totals is a number, and almost never the
#' total of any actual probe. The probe returned is the one whose total is
#' nearest that number, which is the question a chart is asking: show me a probe
#' that looks typical, or extreme. The value is returned alongside so the caller
#' can see how near the nearest was - a `MEDIAN` whose probe sits far from the
#' median total says the distribution is sparse there, and that is worth knowing
#' before reading the chart as representative.
#'
#' @section Ties:
#' The first probe in the order of the pivot wins. The alternative, failing on a
#' tie, would fail constantly: a marker that counts events has many probes with
#' a total of zero, so `MIN` is a tie by construction on any real study.
#'
#' @param marker character. Marker name, e.g. `"MUTATIONS"`, `"DELTAS"`.
#' @param figure character. Figure to sum: `"HYPER"`, `"HYPO"` or `"BOTH"`.
#' @param statistics character vector. Any of `"MIN"`, `"Q1"`, `"MEDIAN"`,
#'   `"Q3"`, `"MAX"`. All five by default.
#'
#' @return A data.frame with `STATISTIC`, `PROBE` and `VALUE`, one row per
#'   statistic asked for, or `NULL` when the pivot does not exist.
#'
#' @seealso [sem_marker_value_per_sample_plot()], which draws the probe chosen
#'   here.
#'
#' @examples
#' # Reads the probe pivots written by a previous semseeker() run; sample_sheet
#' # is the one that run was given.
#' \dontrun{
#' picked <- sem_probe_select_by_statistic("DELTAS", "BOTH", c("MEDIAN", "MAX"))
#' picked
#' sem_marker_value_per_sample_plot(
#'   marker = "DELTAS", figure = "BOTH",
#'   probe = picked$PROBE[picked$STATISTIC == "MAX"],
#'   sample_sheet = sample_sheet
#' )
#' }
#' @export
sem_probe_select_by_statistic <- function(marker, figure,
                                         statistics = c("MIN", "Q1", "MEDIAN",
                                                        "Q3", "MAX")) {

  known <- c("MIN", "Q1", "MEDIAN", "Q3", "MAX")
  unknown <- setdiff(toupper(statistics), known)
  if (length(unknown) > 0L)
    stop("unknown statistic: ", paste(unknown, collapse = ", "),
         ". Known: ", paste(known, collapse = ", "))
  statistics <- toupper(statistics)

  pivot_lazy <- io_read_pivot(marker, figure, area = "PROBE", subarea = "WHOLE")
  if (is.null(pivot_lazy))
    return(NULL)

  pivot <- as.data.frame(pivot_lazy$collect())
  sample_columns <- setdiff(colnames(pivot), "AREA")
  if (nrow(pivot) == 0L || length(sample_columns) == 0L)
    return(NULL)

  probe_total <- rowSums(
    as.matrix(vapply(pivot[sample_columns], as.numeric, numeric(nrow(pivot)))),
    na.rm = TRUE)
  probe_name <- as.character(pivot$AREA)

  quantile_of <- c(MIN = 0, Q1 = 0.25, MEDIAN = 0.5, Q3 = 0.75, MAX = 1)

  do.call(rbind, lapply(statistics, function(statistic) {
    target <- stats::quantile(probe_total, quantile_of[[statistic]],
                              names = FALSE, na.rm = TRUE)
    # which.min() takes the first minimum, so a tie resolves to the first probe
    # in the order of the pivot.
    nearest <- which.min(abs(probe_total - target))
    data.frame(STATISTIC = statistic, PROBE = probe_name[nearest],
               VALUE = probe_total[nearest], stringsAsFactors = FALSE)
  }))
}
