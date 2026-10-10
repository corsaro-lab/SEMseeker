# =============================================================================
# Golden master: the whole output of a fixed run, and how to compare it.
#
# The unit tests say that each function does what it claims. They did not say
# that the output as a whole stays the same, and two defects lived behind a
# green suite for that reason: a SAMPLE-scope test that lost its first sample,
# and a transformation_x that never reached a model.
#
# .golden_run() runs semseeker() and association_analysis() on the bundled
# GSE133774 fixture with a fixed request. .golden_fingerprint() reduces every
# parquet and CSV the run wrote to something that can be stored and compared:
#   - small tables keep their values, so a failure names the cell that moved;
#   - large tables keep per-column summaries (n, non-zero, sum, sum of squares,
#     min, max), which change when any value changes and do not depend on the
#     last bits a different platform may produce.
# Columns are sorted by name and rows by their content, so neither order counts.
#
# Every name starts with a dot: setup.R runs rm(list = ls()) after the helpers
# are sourced, and only hidden names survive it.
#
# It says whether a result CHANGED, not whether it is right. When a change is
# intended, dev/golden-update.R rewrites the reference and the diff of the
# reference in the pull request shows which results moved.
# =============================================================================

.golden_dir <- function() testthat::test_path("golden")

# Cells kept by value up to this size; above it, per-column summaries.
.golden_small <- 50000L

# Columns that describe the run rather than the result.
.golden_volatile <- c("TIME", "TIME_TAKEN", "DURATION", "START_TIME", "END_TIME",
                      "ELAPSED", "SESSION", "SESSION_ID", "PATH", "FILE", "HOST",
                      "CREATED", "DATE", "X")

.golden_input <- function() {
  env <- new.env()
  utils::data("test_master_features", package = "SEMseeker", envir = env)
  utils::data("test_signal_gse133774", package = "SEMseeker", envir = env)
  utils::data("test_samplesheet_gse133774", package = "SEMseeker", envir = env)
  features <- as.data.frame(env$test_master_features, stringsAsFactors = FALSE)
  common <- intersect(features$PROBE, rownames(env$test_signal_gse133774))
  signal <- as.data.frame(env$test_signal_gse133774[common, , drop = FALSE])

  # Deterministic phenotypes: no random number reaches the reference.
  samples <- env$test_samplesheet_gse133774
  k <- seq_len(nrow(samples))
  samples$Phenotest   <- k * 1.5 + (k %% 3)
  samples$Covariates1 <- 40 + (k * 7) %% 23
  samples$Group2      <- ifelse(k %% 2 == 0, "A", "B")
  list(signal = signal, samples = samples)
}

.golden_requests <- function() {
  base <- data.frame(independent_variable = "Phenotest", family_test = "spearman",
                     transformation_y = "none", transformation_x = "none",
                     covariates = "", aggregation = "SUM", scope = "SAMPLE",
                     filter_p_value = FALSE, stringsAsFactors = FALSE)
  # a group family on a two-level variable
  groups <- base
  groups$independent_variable <- "Group2"; groups$family_test <- "wilcoxon"
  # a regression family, with a covariate and a transformation of x
  regression <- base
  regression$family_test <- "gaussian"; regression$covariates <- "Covariates1"
  regression$transformation_x <- "log10"
  # one instance per region
  instance <- base; instance$scope <- "INSTANCE"
  rbind(base, groups, regression, instance)
}

.golden_run <- function(result_folder) {
  input <- .golden_input()
  SEMseeker::semseeker(
    input             = input$signal,
    sample_sheet      = input$samples,
    result_folder     = result_folder,
    parallel_strategy = "sequential",
    areas             = "POSITION",
    markers           = c("SIGNAL", "MUTATIONS", "LESIONS", "DELTAS", "DELTAR",
                          "DELTAQ", "DELTARQ", "DELTAP", "DELTARP"),
    start_fresh       = TRUE,
    inpute            = "median",
    showprogress      = FALSE,
    verbosity         = 0
  )
  requests <- .golden_requests()
  for (r in seq_len(nrow(requests))) {
    instance <- identical(requests$scope[r], "INSTANCE")
    SEMseeker::association_analysis(
      inference_details = requests[r, , drop = FALSE],
      result_folder     = result_folder,
      parallel_strategy = "sequential",
      markers           = c("MUTATIONS", "LESIONS", "DELTAS"),
      figures           = c("HYPER", "HYPO"),
      areas             = if (instance) "GENE" else c("POSITION", "GENE"),
      subareas          = if (instance) "WHOLE" else c("WHOLE", "TSS1500"),
      multiple_test_adj = "BH",
      showprogress      = FALSE,
      verbosity         = 0
    )
  }
  invisible(result_folder)
}

.golden_read <- function(path) {
  df <- if (grepl("\\.parquet$", path)) {
    as.data.frame(polars::pl$read_parquet(path))
  } else {
    utils::read.csv2(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  # unnamed columns are row numbers written by write.csv, not results
  keep <- nzchar(colnames(df)) & !(colnames(df) %in% .golden_volatile)
  df <- df[, keep, drop = FALSE]
  # radix: the C order whatever the locale. testthat runs with LC_COLLATE=C, a
  # shell does not, and the default sort() would order the columns differently.
  df <- df[, sort(colnames(df), method = "radix"), drop = FALSE]
  for (cl in colnames(df)) {
    if (is.factor(df[[cl]])) df[[cl]] <- as.character(df[[cl]])
    if (is.numeric(df[[cl]])) df[[cl]] <- signif(df[[cl]], 10)
  }
  if (nrow(df) > 1L && ncol(df) > 0L)
    df <- df[do.call(order, c(unname(as.list(df)), list(method = "radix"))), , drop = FALSE]
  rownames(df) <- NULL
  df
}

.golden_summary <- function(df) {
  num <- colnames(df)[vapply(df, is.numeric, logical(1))]
  data.frame(
    column   = num,
    n        = vapply(num, function(c) sum(!is.na(df[[c]])), numeric(1)),
    non_zero = vapply(num, function(c) sum(df[[c]] != 0, na.rm = TRUE), numeric(1)),
    sum      = vapply(num, function(c) sum(df[[c]], na.rm = TRUE), numeric(1)),
    sum_sq   = vapply(num, function(c) sum(df[[c]]^2, na.rm = TRUE), numeric(1)),
    min      = vapply(num, function(c) suppressWarnings(min(df[[c]], na.rm = TRUE)), numeric(1)),
    max      = vapply(num, function(c) suppressWarnings(max(df[[c]], na.rm = TRUE)), numeric(1)),
    stringsAsFactors = FALSE, row.names = NULL)
}

# Every parquet and CSV under Data/ and Inference/, by path relative to the run.
.golden_files <- function(result_folder) {
  files <- list.files(result_folder, pattern = "\\.(parquet|csv)$",
                      recursive = TRUE, full.names = FALSE)
  sort(files[grepl("^(Data|Inference)/", files)], method = "radix")
}

.golden_fingerprint <- function(result_folder) {
  out <- list()
  for (f in .golden_files(result_folder)) {
    df <- .golden_read(file.path(result_folder, f))
    entry <- list(nrow = nrow(df), columns = colnames(df))
    if (nrow(df) * max(1L, ncol(df)) <= .golden_small) entry$values <- df
    else entry$summary <- .golden_summary(df)
    out[[f]] <- entry
  }
  out
}

.golden_write <- function(fingerprint, dir = .golden_dir()) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(fingerprint, file.path(dir, "fingerprint.rds"), version = 2)
  # A readable index, so a pull request shows which artefacts exist and how big.
  index <- data.frame(
    file    = names(fingerprint),
    rows    = vapply(fingerprint, function(e) e$nrow, numeric(1)),
    columns = vapply(fingerprint, function(e) length(e$columns), numeric(1)),
    stored  = vapply(fingerprint, function(e) if (is.null(e$values)) "summary" else "values",
                     character(1)),
    stringsAsFactors = FALSE)
  utils::write.csv(index, file.path(dir, "index.csv"), row.names = FALSE)
  invisible(index)
}

.golden_read_reference <- function(dir = .golden_dir()) {
  readRDS(file.path(dir, "fingerprint.rds"))
}

# Differences between two fingerprints, one line each; empty when they agree.
.golden_compare <- function(reference, current, tolerance = 1e-8) {
  msg <- character(0)
  missing <- setdiff(names(reference), names(current))
  added   <- setdiff(names(current), names(reference))
  if (length(missing)) msg <- c(msg, paste("no longer written:", missing))
  if (length(added))   msg <- c(msg, paste("newly written:", added))

  close <- function(a, b) {
    if (is.numeric(a) && is.numeric(b))
      return(isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE)))
    identical(as.character(a), as.character(b))
  }

  for (f in intersect(names(reference), names(current))) {
    r <- reference[[f]]; c <- current[[f]]
    if (!identical(r$columns, c$columns)) {
      msg <- c(msg, sprintf("%s: columns differ (-%s +%s)", f,
                            paste(setdiff(r$columns, c$columns), collapse = ","),
                            paste(setdiff(c$columns, r$columns), collapse = ",")))
      next
    }
    if (!identical(r$nrow, c$nrow)) {
      msg <- c(msg, sprintf("%s: %d rows, reference %d", f, c$nrow, r$nrow))
      next
    }
    if (!is.null(r$values) && !is.null(c$values)) {
      for (cl in r$columns) {
        a <- r$values[[cl]]; b <- c$values[[cl]]
        if (!close(a, b)) {
          i <- which(!(mapply(function(x, y) close(x, y), a, b)))[1]
          msg <- c(msg, sprintf("%s: column %s row %d is %s, reference %s", f, cl, i,
                                format(b[i]), format(a[i])))
        }
      }
    } else if (!is.null(r$summary) && !is.null(c$summary)) {
      for (k in seq_len(nrow(r$summary))) {
        a <- unlist(r$summary[k, -1]); b <- unlist(c$summary[k, -1])
        if (!close(a, b))
          msg <- c(msg, sprintf("%s: column %s summary moved (sum %s, reference %s)", f,
                                r$summary$column[k], format(b[["sum"]]), format(a[["sum"]])))
      }
    } else {
      msg <- c(msg, sprintf("%s: stored as values in one and as summary in the other", f))
    }
  }
  msg
}
