# Prepare a SIGNAL matrix for one sem_analyze_batch() invocation

Single source of truth for normalising a raw input matrix into a shape
that is consistent with the probe-level annotation used downstream. The
invariant guaranteed by this function is:

## Usage

``` r
sem_prepare_batch_signal(
  signal_data,
  tech = NULL,
  sex_chromosome_remove = TRUE
)
```

## Arguments

- signal_data:

  A data.frame whose rownames are probe identifiers (Illumina probe IDs
  like \`"cg00050873"\` for K27/K450/K850, or coordinate-encoded
  \`"\<CHR\>\_\<START\>"\` strings for WGBS/LONGREAD). Sample columns
  follow. Must be PROBE-keyed, i.e. already passed through
  \`io_normalize_signal_input()\`.

- tech:

  Character scalar. One of \`"K27"\`, \`"K450"\`, \`"K850"\`,
  \`"WGBS"\`, \`"LONGREAD"\`. If \`NULL\` (default) the function calls
  \`core_get_meth_tech()\` to detect it.

- sex_chromosome_remove:

  Logical. If \`TRUE\` (default), drop probes on \`CHR == "X"\` or \`CHR
  == "Y"\` from both \`probe_features\` and \`signal_data\`. Applied
  uniformly across all techs.

## Value

The input \`signal_data\` filtered to the autosomal manifest
intersection, with two attributes attached:

- \`probe_features\` - data.frame with columns \`PROBE, CHR, START,
  END\` (+ any extra columns produced by the tech-specific annotation
  builder), one row per surviving probe, in the same order as
  \`rownames(signal_data)\`.

- \`tech\` - the resolved tech string.

## Details

\`nrow(signal_data) == nrow(attr(signal_data, "probe_features"))\`

and \`rownames(signal_data)\` is exactly \`attr(.,
"probe_features")\$PROBE\` in the same order.

This replaces ~30 lines of scattered annotation/filter/align logic that
previously lived in \`sem_analyze_batch()\` fresh-path and was a source
of silent drift between the SIGNAL matrix and the probe_features used to
compute thresholds, write the POSITION pivot, and run downstream
analyses. The classic failure (visible v35-v43) was

\`Error in data.frame(probe_features, VALUE = values, row.names =
probe_features\$PROBE): arguments imply differing number of rows:
366948, 374705\`

which the alignment step in this function eliminates by construction.

## Invariants

Post-call:

- \`nrow(signal_data) == nrow(attr(signal_data, "probe_features"))\`

- \`identical(rownames(signal_data), attr(signal_data,
  "probe_features")\$PROBE)\`

- \`!any(attr(signal_data, "probe_features")\$CHR when
  \`sex_chromosome_remove = TRUE\`

- \`!anyDuplicated(attr(signal_data, "probe_features")\$PROBE)\`
