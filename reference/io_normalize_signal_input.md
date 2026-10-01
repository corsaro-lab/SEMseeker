# Transparently normalise signal input for any technology.

If the input data frame has CHR/START columns (WGBS or long-read
format), it is converted to the probe-ID-indexed format expected by
SEMseeker. Illumina matrices (rownames = probe IDs) are returned
unchanged.

## Usage

``` r
io_normalize_signal_input(signal_data)
```

## Arguments

- signal_data:

  A data frame (probe-indexed or coordinate-based).

## Value

A data frame with rownames = probe IDs (real or synthetic).

## Details

Called inside `sem_analyze_batch()` before
[`core_get_meth_tech()`](https://corsaro-lab.github.io/SEMseeker/reference/core_get_meth_tech.md).
