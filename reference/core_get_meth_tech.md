# Detect the Illumina methylation array technology from a signal matrix

Identifies whether a methylation signal matrix originates from a 27k,
450k, or EPIC 850k array (or WGBS data). Detection runs in priority
order:

## Usage

``` r
core_get_meth_tech(signal_data)
```

## Arguments

- signal_data:

  A numeric matrix or `data.frame` with CpG probes as rows and samples
  as columns. Row names must be probe identifiers (e.g. `cg00000029`)
  unless a `PROBE` column is present.

## Value

The updated session environment (`ssEnv`), invisibly. The detected
technology is accessible via `core_get_session_info()$tech`.

## Details

1.  **Annotation-package overlap** - probe IDs are matched against each
    installed `IlluminaHumanMethylation*anno` package; the technology
    with the most matching probes wins. Accurate for any subset size
    (e.g. 20 k probes out of 866 k).

2.  **Row-count heuristics** - last resort for WGBS or unknown probe
    naming schemes.

The detected technology and beta/M-value flag are written to the session
environment.
