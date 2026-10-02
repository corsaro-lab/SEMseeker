# Build a canonical output file path

All SEMseeker output files have their filename UPPERCASED via
\[core_name_cleaning()\] (which calls \`toupper()\`). This is
intentional and load-bearing: it guarantees that semantic identifiers
(AREA / MARKER / FIGURE / Sample_ID) collapse to a stable case
regardless of how the caller spelled them, so pivot / per-sample-bed /
summary files always resolve to the same on-disk path.

## Usage

``` r
io_file_path_build(baseFolder, detailsFilename, extension, add_gz = FALSE)
```

## Arguments

- baseFolder:

  Directory the file lives in.

- detailsFilename:

  Character vector concatenated with "\_" and passed to
  \[core_name_cleaning()\] (uppercased + non-alnum → "\_").

- extension:

  File extension (no leading dot). Empty string skips.

- add_gz:

  If TRUE, append ".gz".

## Value

Full path as a character scalar.

## Details

Practical consequence: a literal name like \`"sample_sheet_result"\` is
written to disk as \`SAMPLE_SHEET_RESULT.csv\`. On case-INsensitive file
systems (macOS APFS/HFS, Windows NTFS) \`file.exists()\` finds either
spelling; on case-SENSITIVE ones (Linux ext4) only the uppercase form
resolves. Tests that hard-code an expected path must therefore use the
uppercase form, or - preferably - discover the file via \[list.files()\]
or by re-calling \`io_file_path_build()\` with the same arguments.
