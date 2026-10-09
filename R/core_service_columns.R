#--- core_service_columns ---------------------------------------------------------
# The one list of columns that are not data.
#
# A pivot or a sample table carries, next to its values, columns that locate or
# identify a row: genomic coordinates, probe and area names, array membership,
# sample identity. They are never a sample, an area or a burden, and they can
# sit first or last depending on how the table was built (group_by puts AREA
# first, with_columns puts it last). Reading the data by position therefore
# picks the wrong columns as soon as the layout changes; reading it as "every
# column but these" does not.
#
# Every place that separates data from service columns reads this list, so a
# new service column (an Is_Reference flag replacing Sample_Group, for
# instance) is declared once.
#
# @return character vector of service column names.
# @keywords internal
core_service_columns <- function() {
  c(
    # genomic coordinates
    "CHR", "START", "END",
    # probe and area identifiers
    "PROBE", "PROBE_WHOLE", "AREA",
    # array membership flags, in both spellings the manifests use
    "K27", "K450", "K850", "k27", "k450", "k850",
    # sample identity
    "Sample_ID", "Sample_Group"
  )
}

# The data columns of a table: every column that is not a service column.
#
# @param x a data.frame, a polars frame or a character vector of column names.
# @param also further names to exclude (the independent variable and the
#   covariates of a model, for instance).
# @return character vector of column names, in the table's order.
# @keywords internal
core_data_columns <- function(x, also = character(0)) {
  nm <- if (is.character(x)) x else names(x)
  setdiff(nm, c(core_service_columns(), also))
}
