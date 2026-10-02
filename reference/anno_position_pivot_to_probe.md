# Convert a POSITION-keyed signal pivot into a PROBE-keyed lazy frame.

Returns a Polars LazyFrame rather than an R data.frame. Materialisation
downstream is opt-in only - sem_analyze_batch and
sem_analyze_population_bulk consume the lazy frame directly. This
eliminates the ~12 GB R-side peak that caused silent jetsam kills on
ewas-scale matrices (367k × 4013) in resume mode (v18, v21).

## Usage

``` r
anno_position_pivot_to_probe(signal_data)
```

## Arguments

- signal_data:

  polars LazyFrame, polars DataFrame, or R data.frame carrying
  CHR/START/END + sample columns. Lazy input keeps the entire pipeline
  in zero-copy mode; R data.frame is wrapped before the join.

## Value

A polars LazyFrame with the PROBE column and one column per sample.
CHR/START/END/tech annotation columns are dropped after the join, since
downstream consumers only need the probe-keyed shape.

## Details

Contract change: callers that need PROBE as a column see it as a regular
polars String column. Callers that previously read \`rownames\` must
switch to selecting / collecting the \`PROBE\` column.

Order: input must already be sorted (sort gate = io_signal_save chunked
per-chr by START). This function preserves order via
\`unique(maintain_order=TRUE)\`.
