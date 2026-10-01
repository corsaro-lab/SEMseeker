# Core SEMseeker pipeline (internal)

Internal entry point used by the public
[`semseeker`](https://corsaro-lab.github.io/SEMseeker/reference/semseeker.md)
dispatcher. Accepts already-normalised signal data (matrix or data frame
with probe-ID rownames) and runs the full SEM analysis. Users should
call
[`semseeker`](https://corsaro-lab.github.io/SEMseeker/reference/semseeker.md)
instead, which handles input normalisation, M-value conversion, and
tech/genome_build validation.

## Usage

``` r
sem_core(sample_sheet, signal_data, result_folder)
```

## Arguments

- sample_sheet:

  Data frame (or list of data frames) with a `Sample_ID` column.

- signal_data:

  Methylation matrix or data frame with probe-ID rownames (or list of
  such objects, one per batch).

- result_folder:

  Output directory.

## Value

Invisibly NULL; writes output files to `result_folder`.
