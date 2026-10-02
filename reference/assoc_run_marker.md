# Run depth\>1 (area / region level) association for one marker

Extracted from association_analysis() (was inline at lines 294-418).
Iterates over the area/subarea keys of one marker, reads the
corresponding pivot parquet, chunks it (6e6 cells / ncol), transposes
and merges with sample_names, then applies the stat model per chunk.

## Usage

``` r
assoc_run_marker(
  prep,
  marker,
  family_test,
  fileNameResults,
  filter_p_value,
  ssEnv,
  selected_areas,
  results,
  start_time,
  processed_items,
  ...
)
```

## Arguments

- prep:

  list returned by sem_prepare_study_for_analysis().

- marker:

  character. The marker name (e.g. "MUTATIONS").

- family_test:

  character.

- fileNameResults:

  character. Path of the output CSV.

- filter_p_value:

  logical.

- ssEnv:

  list. Session environment.

- selected_areas:

  character vector or empty.

- results:

  data.frame. Accumulator carried over from depth=1.

- start_time:

  POSIXct. Job start, used by assoc_analysis_log().

- processed_items:

  integer. Counter carried over from depth=1.

- ...:

  forwarded to assoc_apply_stat_model().

## Value

list(results = data.frame, processed_items = integer). Side effect:
writes the CSV via assoc_analysis_save_results().
