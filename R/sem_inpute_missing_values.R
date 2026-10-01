# AI-096 Phase 2 (2026-06-09): KNN imputation memory gate.
#
# Estimates the peak RAM needed for KNN on the input matrix and compares
# it to the available system RAM (scaled by SEMSEEKER_KNN_MEM_FRACTION,
# default 0.6 to leave room for Polars cache + R working set + OS + apps).
# Fail-fast with an actionable error before the kmeans_knn() call so the
# user knows to switch to median/mean or run KNN upstream of semseeker().
#
# Memory model:
#   matrix    = n_probes × n_samples × 8 bytes        (double)
#   distance  = n_samples × n_samples × 8 bytes       (sample-sample)
#   working   = 2 × matrix                            (KNN intermediate buffers)
#   needed    = matrix + distance + working
#
# All in bytes; converted to GB at log time.
.sem_knn_memory_gate <- function(mat) {
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    return(invisible(NULL))   # only gate matrices/frames
  }
  n_probes  <- nrow(mat)
  n_samples <- ncol(mat)
  byte_per_double <- 8
  matrix_GB   <- n_probes  * n_samples  * byte_per_double / (1024^3)
  distance_GB <- n_samples * n_samples  * byte_per_double / (1024^3)
  working_GB  <- 2 * matrix_GB
  needed_GB   <- matrix_GB + distance_GB + working_GB

  # SEMSEEKER_KNN_MEM_FRACTION (env var) overrides the default 0.6.
  mem_frac <- suppressWarnings(as.numeric(Sys.getenv("SEMSEEKER_KNN_MEM_FRACTION", "0.6")))
  if (is.na(mem_frac) || mem_frac <= 0 || mem_frac > 1) mem_frac <- 0.6

  # Available RAM: prefer memuse if installed, else parse /proc/meminfo
  # on linux, fall back to `sysctl hw.memsize` on macOS. If everything
  # fails, skip the gate (don't block the user on a portability issue).
  total_GB <- .sem_total_ram_GB()
  if (is.na(total_GB) || total_GB <= 0) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " sem_inpute_missing_values: could not detect total RAM - skipping KNN memory gate.")
    return(invisible(NULL))
  }
  available_GB <- total_GB * mem_frac

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " sem_inpute_missing_values: KNN memory estimate matrix=",
            round(matrix_GB, 2), " GB + distance=", round(distance_GB, 3),
            " GB + buffers=", round(working_GB, 2), " GB = ",
            round(needed_GB, 2), " GB; available=", round(available_GB, 2),
            " GB (", round(mem_frac * 100), "% of ", round(total_GB), " GB).")

  if (needed_GB > available_GB) {
    msg <- sprintf(
      "KNN imputation needs ~%.1f GB (matrix=%.1f, distance=%.2f, buffers=%.1f), available ~%.1f GB (%.0f%% of total %.0f GB). Options:\n  - use inpute='median' (matrixStats row-median, C-level, no full materialization)\n  - use inpute='mean' (same shape, mean instead of median)\n  - run KNN upstream and pass the imputed matrix to semseeker()\n  - raise SEMSEEKER_KNN_MEM_FRACTION (current %.2f) if you have headroom",
      needed_GB, matrix_GB, distance_GB, working_GB,
      available_GB, mem_frac * 100, total_GB, mem_frac
    )
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
              " sem_inpute_missing_values: KNN gate triggered. ", msg)
    stop(msg, call. = FALSE)
  }
  invisible(NULL)
}

# Best-effort total RAM detection. Returns GB or NA.
.sem_total_ram_GB <- function() {
  # 1) memuse (cross-platform, optional dependency)
  if (requireNamespace("memuse", quietly = TRUE)) {
    bytes <- tryCatch(as.numeric(memuse::Sys.meminfo()$totalram),
                       error = function(e) NA_real_)
    if (!is.na(bytes) && bytes > 0) return(bytes / (1024^3))
  }
  # 2) macOS: sysctl hw.memsize → bytes
  os <- Sys.info()[["sysname"]]
  if (os == "Darwin") {
    out <- tryCatch(suppressWarnings(system2("sysctl", c("-n", "hw.memsize"),
                                              stdout = TRUE, stderr = FALSE)),
                     error = function(e) NULL)
    if (!is.null(out) && length(out) == 1L) {
      bytes <- suppressWarnings(as.numeric(out))
      if (!is.na(bytes) && bytes > 0) return(bytes / (1024^3))
    }
  }
  # 3) Linux: /proc/meminfo MemTotal (kB)
  if (os == "Linux" && file.exists("/proc/meminfo")) {
    lines <- tryCatch(readLines("/proc/meminfo", n = 1L), error = function(e) NULL)
    if (!is.null(lines)) {
      m <- regmatches(lines, regexec("MemTotal:\\s+(\\d+)\\s+kB", lines))[[1]]
      if (length(m) >= 2L) {
        kB <- suppressWarnings(as.numeric(m[2]))
        if (!is.na(kB) && kB > 0) return(kB * 1024 / (1024^3))
      }
    }
  }
  # 4) Windows: PowerShell CIM query → bytes
  if (os == "Windows") {
    out <- tryCatch(
      suppressWarnings(system2(
        "powershell",
        c("-Command", "(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory"),
        stdout = TRUE, stderr = FALSE)),
      error = function(e) NULL)
    if (!is.null(out) && length(out) >= 1L) {
      bytes <- suppressWarnings(as.numeric(trimws(out[[1L]])))
      if (!is.na(bytes) && bytes > 0) return(bytes / (1024^3))
    }
  }
  NA_real_
}

sem_inpute_missing_values <- function(signal_data){

  ssEnv <- core_get_session_info()
  # count number of na
  n_na <- sum(is.na(signal_data))

  if (n_na==0)
    return (signal_data)

  nrow_ex_ante <- nrow(signal_data)

  # count missed per rows
  n_missed_per_row <- rowSums(is.na(signal_data))
  if (any(n_missed_per_row > 0.1 * ncol(signal_data))) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " There are rows with missing values more than 10%. I will remove them.")
    core_log_event("JOURNAL: ", format(Sys.time(), "%a %b %d %X %Y"), " Imputing missing values using ", ssEnv$inpute , " method. Number of missing values: ", n_na, " corresponding to: ", round(n_na/(nrow(signal_data)*ncol(signal_data)), 2), " % of the data.")
    signal_data <- signal_data[n_missed_per_row  < 0.1 * ncol(signal_data), ]
  }

  n_missed_per_col <- colSums(is.na(signal_data))
  if (any(n_missed_per_col > 0.1 * nrow(signal_data))) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " There are columns with missing values more than 10%. I will remove them.")
    core_log_event("JOURNAL: ", format(Sys.time(), "%a %b %d %X %Y"), " There are columns with missing values more than 10%. I will remove them.")
    signal_data <- signal_data[, n_missed_per_col  < 0.1 * nrow(signal_data)]
  }

  chunk_size <- 10000  # Define a chunk size
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = seq_along(seq(1, nrow(signal_data), by = chunk_size)))
  else
    progress_bar <- ""

  n_item <- nrow(signal_data)*ncol(signal_data)/100
  core_log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y") ," Imputing missing values using ", ssEnv$inpute , " method. Number of missing values: ", n_na, " corresponding to: ", round(n_na/n_item, 2), " % of the data.")
  core_log_event("JOURNAL: Imputing missing values using ", ssEnv$inpute , " method. Number of missing values: ", n_na, " corresponding to: ", round(n_na/n_item, 2), " % of the data.")
  mat <- as.matrix(signal_data)

  if (ssEnv$inpute == "median")
  {
    # Vectorized: matrixStats::rowMedians is C-level, ~10x faster than apply
    row_med <- matrixStats::rowMedians(mat, na.rm = TRUE)
    na_idx  <- which(is.na(mat), arr.ind = TRUE)
    mat[na_idx] <- row_med[na_idx[, 1L]]
    signal_data <- as.data.frame(mat)
  }
  else if (ssEnv$inpute == "mean")
  {
    row_mn <- matrixStats::rowMeans2(mat, na.rm = TRUE)
    na_idx <- which(is.na(mat), arr.ind = TRUE)
    mat[na_idx] <- row_mn[na_idx[, 1L]]
    signal_data <- as.data.frame(mat)
  }
  else if (grepl("knn", ssEnv$inpute))
  {
    if (length(strsplit(ssEnv$inpute, ";")[[1]]) != 3)
    {
      core_log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y") ," Invalid inpute value. Please provide the number of centers and the number of clusters separated by a semicolon.")
      stop()
    }

    # AI-096 Phase 2 (2026-06-09): memory gate for KNN imputation.
    # KNN requires the full matrix + a sample×sample distance matrix
    # + working buffers. On ewas-scale (367k × 4013 ≈ 12 GB) the
    # 64 GB Mac can absorb it; on long-reads (10⁶ × variable) it
    # blows up to 30+ GB and macOS jetsam kills R silently. Fail-fast
    # before the call so the user gets an actionable error.
    .sem_knn_memory_gate(mat)

    centers <- strsplit(ssEnv$inpute, ";")[[1]][2]
    k <- strsplit(ssEnv$inpute, ";")[[1]][3]
    imputed_matrix <- KMEANS.KNN::kmeans_knn(mat, centers = centers, k = k)
    signal_data <- as.data.frame(imputed_matrix)
  }
  else if (ssEnv$inpute == "none")
  {
    # Vectorized: rowAnys(is.na) avoids creating full logical matrix copy
    has_na <- matrixStats::rowAnyNAs(mat)
    signal_data <- signal_data[!has_na, ]
  }
  else
  {
    core_log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y") ," Invalid inpute value. Please provide one of the following: median, mean, knn.")
    stop()
  }

  rm(mat)
  gc()

  # drop rows with all NA (vectorized)
  all_na <- matrixStats::rowAlls(is.na(as.matrix(signal_data)))
  signal_data <- signal_data[!all_na, ]
  nrows_ex_post <- nrow(signal_data)
  if (nrows_ex_post < nrow_ex_ante)
  {
    core_log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y") ," Dropping rows with all missing values. Number of rows dropped: ", nrow_ex_ante - nrows_ex_post)
    core_log_event("JOURNAL: Dropped rows with all missing values. Number of rows dropped: ", nrow_ex_ante - nrows_ex_post)
  }
  n_na <- sum(is.na(signal_data))
  core_log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y") ," Imputation completed. Number of missing values: ", n_na)
  return(signal_data)
}
