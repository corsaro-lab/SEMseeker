# init ssEnvonment

init ssEnvonment

## Usage

``` r
core_init_env(result_folder, maxResources = 90, ...)
```

## Arguments

- result_folder:

  where result of semseeker will be stored

- maxResources:

  percentage of how many available cores will be used (default 90
  percent, rounded to lowest integer)

- ...:

  additional session options, including:

  `parallel_strategy`

  :   parallelisation strategy for future: `"sequential"` (default),
      `"multisession"`, `"multicore"`, `"cluster"`

  `genome_build`

  :   reference genome assembly: `"hg19"` (default, matches Illumina
      array annotation), `"hg38"` (GRCh38, typical for long-read /
      Nanopore data), `"mm10"` (mouse - requires C-05). Stored in
      `ssEnv$genome_build` and written to session provenance metadata
      (C-06).

  `tech`

  :   override technology auto-detection: `"K850"`, `"K450"`, `"K27"`,
      `"WGBS"`, `"LONGREAD"`. Required for long-read data because
      LONGREAD cannot be distinguished from WGBS by probe-ID pattern
      alone. Example:
      `core_init_env(folder, tech = "LONGREAD", genome_build = "hg38")`

## Value

the working ssEnvonment
