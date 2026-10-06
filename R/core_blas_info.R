# BLAS detection + one-line WARNING at core_init_env() time.
#
# R ships with the single-thread reference BLAS by default. The limma/voom batch families
# limma/voom families spend most of their wall-clock in BLAS calls
# (solve / crossprod inside lmFit + eBayes) - they lose ~5-10x against
# a multi-threaded BLAS (Accelerate on macOS, OpenBLAS or MKL on Linux,
# OpenBLAS Rblas.dll on Windows). The user can't tell the difference
# without being told: the code runs, returns the right numbers, just
# on one core out of N.
#
# This file only DETECTS and WARNS - no user-facing helper, no system
# mutation. The fix command is printed in the warning itself,
# OS-specific.

.core_blas_info <- function() {
  blas_path <- tryCatch(unname(extSoftVersion()["BLAS"]),
                         error = function(e) "")
  if (!is.character(blas_path) || length(blas_path) == 0L) blas_path <- ""

  is_accelerate <- grepl("vecLib|Accelerate", blas_path, ignore.case = TRUE)
  is_openblas   <- grepl("openblas",          blas_path, ignore.case = TRUE)
  is_mkl        <- grepl("mkl",               blas_path, ignore.case = TRUE)
  is_reference  <- grepl("libRblas(\\.[0-9])?\\.(dylib|so|dll)$|Rblas\\.dll$",
                          blas_path)

  flavor <- if (is_accelerate) "Accelerate (vecLib)"
            else if (is_openblas)  "OpenBLAS"
            else if (is_mkl)       "MKL"
            else if (is_reference) "reference (single-thread)"
            else "unknown"

  list(path           = blas_path,
       multi_threaded = is_accelerate || is_openblas || is_mkl,
       flavor         = flavor,
       is_reference   = is_reference)
}

.core_warn_blas_single_thread <- function() {
  info <- .core_blas_info()
  if (info$multi_threaded) return(invisible(info))

  sysname <- Sys.info()[["sysname"]]
  rhome   <- R.home()
  fix_cmd <- switch(
    sysname,
    "Darwin" = sprintf(
      "cd %s/lib && sudo ln -sf libRblas.vecLib.dylib libRblas.dylib",
      rhome),
    "Linux"  = paste0(
      "sudo apt-get install -y libopenblas-dev && ",
      "sudo update-alternatives --config libblas.so.3-x86_64-linux-gnu  ",
      "# Debian/Ubuntu; on RHEL/Fedora use dnf install openblas-devel"),
    "Windows" = paste0(
      "Replace ", rhome, "\\bin\\x64\\Rblas.dll with an OpenBLAS build. ",
      "Prebuilt DLLs: https://github.com/xianyi/OpenBLAS/wiki"),
    paste0("Consult your R/BLAS docs for ", sysname, ".")
  )
  core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
            " BLAS is ", info$flavor, " (", info$path,
            "). limma/voom will run single-thread (~5-10x slower than",
            " multi-threaded BLAS). One-time fix: ", fix_cmd)
  invisible(info)
}
