# Run SEMseeker on methylation data from any supported source

Public entry point. Accepts a wide range of inputs — bedmethyl files
(modkit / nanopolish), coordinate-based data frames (WGBS / long-read),
Illumina probe-indexed matrices, or already-loaded data frames —
normalises them to the internal format, validates the tech /
genome_build combination, and delegates to the core pipeline
[`sem_core`](https://corsaro-lab.github.io/SEMseeker/reference/sem_core.md).
Input values are passed through unchanged: if you need to convert
M-values to beta, call
[`sem_mvalue_to_beta`](https://corsaro-lab.github.io/SEMseeker/reference/sem_mvalue_to_beta.md)
explicitly before `semseeker()`.

## Usage

``` r
semseeker(
  input,
  sample_sheet,
  result_folder,
  input_type = c("auto", "bedmethyl", "coord_df", "matrix"),
  tech = NULL,
  genome_build = "hg19",
  strict_build_check = TRUE,
  ...
)
```

## Arguments

- input:

  Input signal data. See details for supported forms.

- sample_sheet:

  Data frame (or list of data frames) with a `Sample_ID` column
  identifying samples.

- result_folder:

  Output directory.

- input_type:

  One of `"auto"` (default), `"bedmethyl"`, `"coord_df"`, `"matrix"`.
  When `"auto"`, the type is inferred from the object class / file
  extension.

- tech:

  Optional technology label (`"K850"`, `"K450"`, `"K27"`, `"WGBS"`,
  `"LONGREAD"`). If `NULL` (default) it is auto-detected downstream by
  [`core_get_meth_tech()`](https://corsaro-lab.github.io/SEMseeker/reference/core_get_meth_tech.md).

- genome_build:

  Reference genome build. One of `"hg19"` (default), `"hg38"`, `"mm10"`,
  `"legacy"`.

- strict_build_check:

  If `TRUE` (default), impossible tech / genome_build combinations (e.g.
  `tech="LONGREAD"` with `genome_build="hg19"`) raise an error. If
  `FALSE`, a warning is emitted and the pipeline proceeds.

- ...:

  Additional arguments forwarded to
  [`sem_core`](https://corsaro-lab.github.io/SEMseeker/reference/sem_core.md)
  and
  [`core_init_env()`](https://corsaro-lab.github.io/SEMseeker/reference/core_init_env.md)
  (e.g. `parallel_strategy`, `alpha`, `LESIONS_BP`, `marker`, `areas`).
  `LESIONS_BP` (default 2000) is the maximum bp distance between two
  probes for them to be in the same LESIONS enrichment window — replaces
  the legacy `sliding_window_size` probe-count parameter (removed in
  AI-092).

## Value

Invisibly `NULL`; writes output files to `result_folder`.

## Details

Supported `input` forms:

- Character vector of bedmethyl file paths (`.bed`/`.tsv`/ `.bedmethyl`)
  — parsed via
  [`io_bedmethyl_read`](https://corsaro-lab.github.io/SEMseeker/reference/io_bedmethyl_read.md).

- Data frame with `CHR`/`START`\[`/END`\] columns (WGBS / long-read
  coordinate format) — normalised via
  [`io_normalize_signal_input`](https://corsaro-lab.github.io/SEMseeker/reference/io_normalize_signal_input.md).

- Matrix or data frame with probe-ID rownames (Illumina array) — passed
  through unchanged.

- List of any of the above, one element per batch.

## Examples

``` r
# Stub: see vignette('imprinting-disorders', package = 'SEMseeker') for a
# runnable Beckwith-Wiedemann workflow on the GSE133774 subset (AI-112b).
invisible(NULL)
if (FALSE) { # \dontrun{
# Bedmethyl (Nanopore / modkit):
semseeker(
  input = list.files("bedmethyl/", pattern = "\\.bed$", full.names = TRUE),
  sample_sheet = ss,
  result_folder = tempdir(),
  tech = "LONGREAD",
  genome_build = "hg38"
)

# Illumina beta matrix:
semseeker(
  input = beta_matrix,
  sample_sheet = ss,
  result_folder = tempdir()
)

# WGBS coordinate data frame:
semseeker(
  input = wgbs_df,              # columns: CHR, START, END, sample1, ...
  sample_sheet = ss,
  result_folder = tempdir(),
  tech = "WGBS"
)
} # }
```
