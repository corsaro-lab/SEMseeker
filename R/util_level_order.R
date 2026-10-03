#' The order of the levels of an ordinal independent variable
#'
#' Pure helper: no I/O except a log line, no side effects. Answers one question -
#' in what order do the levels of this variable stand - and answers it in one
#' place, because the package converts the independent variable to a factor at
#' two points with two different family lists
#' (\code{io_data_preparation()} and \code{sem_prepare_study_for_analysis()}),
#' and an order applied at one of them only would diverge from the other.
#'
#' @section Why the order has to be declarable:
#' Without a declared order the levels are whatever sorting the labels into a
#' sequence produces, and for an ordinal variable that is rarely the order the
#' question is about. A tumour stage reading reference, in situ, I, II, III, IV
#' sorts to I, II, III, in situ, IV, reference - the Roman numerals happen to
#' land correctly and the words do not - so the two comparisons that matter,
#' reference against in situ and in situ against I, are not adjacent at all, and
#' two that mean nothing are. Nothing in the result says which order was used,
#' so a wrong one and a right one read the same.
#'
#' The consequences are not confined to figures. The levels decide which
#' comparisons are consecutive, they decide the sign of a difference between
#' them, and \code{levels()[1]} is the reference category of every regression
#' fitted on the variable.
#'
#' @section The fallback, and what it is not:
#' With no declared order the levels are sorted with \code{gtools::mixedsort},
#' which reads embedded numbers as numbers: \code{"1" < "2" < "3" < "10"} where
#' plain sorting gives \code{"1" < "10" < "2" < "3"}. That is a real improvement
#' for levels that are numbers written as text and it does nothing at all for
#' labels that are words - \code{mixedsort} puts in situ between III and IV. The
#' fallback keeps a request working; it does not make it right. Only the
#' declaration does that.
#'
#' @param values The observed values of the independent variable.
#' @param declared_order Optional. The levels in the order the request intends,
#'   as a \code{+}-separated string or a character vector. \code{NULL}, \code{NA}
#'   or an empty string mean no declaration.
#' @return A character vector of the observed levels, in order.
#'
#' @keywords internal
#' @noRd
util_level_order <- function(values, declared_order = NULL) {

  observed <- unique(as.character(values))
  observed <- observed[!is.na(observed)]

  declared <- util_split_and_clean(declared_order)
  declared <- as.character(declared)
  declared <- declared[!is.na(declared) & nzchar(declared)]

  if (!length(declared))
    return(gtools::mixedsort(observed))

  # An observed level the declaration does not mention is a request that does not
  # say where that level goes. Refusing is the only answer that cannot be wrong:
  # appending it would invent a position, and dropping it would silently remove
  # samples from the analysis.
  unplaced <- setdiff(observed, declared)
  if (length(unplaced))
    stop("independent_variable_order does not place every level present in the data.\n",
         "  Not placed: ", paste(unplaced, collapse = ", "), "\n",
         "  Declared:   ", paste(declared, collapse = ", "), "\n",
         "  Every observed level has to appear, because the order decides which ",
         "comparisons are consecutive, the sign of the difference between them, ",
         "and the reference category of a regression.",
         call. = FALSE)

  # The other direction is legitimate: a samples filter can remove a level the
  # request still names, and the request should not have to be edited per subset.
  unobserved <- setdiff(declared, observed)
  if (length(unobserved))
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " independent_variable_order names ", paste(unobserved, collapse = ", "),
              ", absent from this subset: dropped from the order, the rest keeps its sequence.")

  declared[declared %in% observed]
}
