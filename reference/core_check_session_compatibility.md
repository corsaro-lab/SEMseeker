# Check compatibility of SEMseeker sessions before meta-analysis

Reads `session_metadata.json` from each path in `session_list` and
enforces provenance rules before combining results across studies:

1.  **Stops** if `genome_build` differs - coordinates from different
    assemblies are physically incomparable and would produce silently
    wrong intersection results.

2.  **Warns** if `tech` differs - cross-array meta-analysis (e.g. K450 +
    K850) is statistically valid on the probe intersection but must be
    intentional.

A missing `session_metadata.json` raises a warning (legacy runs without
provenance data are tolerated but flagged).

## Usage

``` r
core_check_session_compatibility(session_list)
```

## Arguments

- session_list:

  character vector. Paths to SEMseeker result folders to compare (must
  be at least two for a meaningful comparison).

## Value

Invisibly returns a `data.frame` with one row per session containing the
`folder`, `genome_build`, and `tech` fields parsed from each metadata
file.

## Examples

``` r
# Two fake session folders with provenance metadata
d1 <- file.path(tempdir(), "study_A"); dir.create(d1, showWarnings = FALSE)
d2 <- file.path(tempdir(), "study_B"); dir.create(d2, showWarnings = FALSE)
writeLines(jsonlite::toJSON(list(genome_build = "hg19", tech = "K450"),
                            auto_unbox = TRUE),
           file.path(d1, "session_metadata.json"))
writeLines(jsonlite::toJSON(list(genome_build = "hg19", tech = "K450"),
                            auto_unbox = TRUE),
           file.path(d2, "session_metadata.json"))
core_check_session_compatibility(c(d1, d2))
```
